#!/usr/bin/env python3
"""
Run a short Si Fig.8 runtime-compliance probe for modes a/b/c and check that
the runtime process/transport setup matches the frozen baseline assumptions.

Outputs:
  - CSV pass/fail table
  - Markdown summary report
  - Raw Geant4 logs and effective per-mode probe JSON configs
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


BASELINE_CONFIGS = {
    "a": Path("config/geant4/scan_si_fig8_paper_mode_a_sub5000nm_r100000nm.json"),
    "b": Path("config/geant4/scan_si_fig8_paper_mode_b_sub5000nm_r100000nm.json"),
    "c": Path("config/geant4/scan_si_fig8_paper_mode_c_sub5000nm_r100000nm.json"),
}

PROHIBITED_PROCESSES = {
    "StepLimiter",
    "eBrem",
    "ePairProd",
    "CoulombScat",
    "eCoulombScattering",
}


@dataclass
class CheckResult:
    check_id: str
    description: str
    passed: bool
    observed: str
    expected: str


def _run_command(cmd: List[str], log_path: Path) -> None:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_path.write_text(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Geant4 run failed with code {proc.returncode}: {' '.join(cmd)}\n"
            f"See log: {log_path}"
        )


def _extract_mode_flags(log_text: str) -> Dict[str, str]:
    out = {}
    mode_line = ""
    for line in log_text.splitlines():
        if "[MicroElec] mode=" in line:
            mode_line = line.strip()
            break
    out["mode_line"] = mode_line
    if not mode_line:
        return out

    # Example:
    # [MicroElec] mode=a | disable_initial_energy=false | surface=true | ... | strict=true | ...
    for token in mode_line.split("|"):
        token = token.strip()
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        out[k.strip().replace("[MicroElec] ", "")] = v.strip()
    return out


def _extract_process_list(log_text: str, particle_name: str = "e-") -> List[str]:
    lines = log_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if f"[ProcessList] particle: {particle_name}" in line:
            start = i
            break
    if start is None:
        return []

    processes: List[str] = []
    for line in lines[start + 1 :]:
        s = line.strip()
        if s.startswith("[ProcessList] particle:"):
            break
        m = re.search(r"\[\s*\d+\]\s+([^\|]+)\|", s)
        if m:
            processes.append(m.group(1).strip())
    return processes


def _extract_global_em_settings(log_text: str) -> Dict[str, str]:
    lines = log_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "--- EM step settings (global) ---" in line:
            start = i
            break
    if start is None:
        return {}

    settings: Dict[str, str] = {}
    for line in lines[start + 1 :]:
        s = line.strip()
        if s.startswith("---"):
            break
        if ":" not in s:
            continue
        k, v = s.split(":", 1)
        settings[k.strip()] = v.strip()
    return settings


def _build_probe_config(
    base_cfg_path: Path,
    mode: str,
    energy_ev: float,
    events: int,
    out_dir_tag: str,
) -> Dict[str, object]:
    cfg = json.loads(base_cfg_path.read_text())
    cfg["primary_energy_MeV"] = [energy_ev * 1.0e-6]
    cfg["events"] = int(events)
    cfg["output_dir"] = (
        f"results/probe_si_fig8_runtime_compliance_mode{mode}_"
        f"{int(round(energy_ev))}eV_events{events}_{out_dir_tag}"
    )
    return cfg


def _mode_expectations(mode: str) -> Dict[str, str]:
    if mode == "a":
        return {"disable_initial_energy": "false", "surface": "true"}
    if mode == "b":
        return {"disable_initial_energy": "true", "surface": "true"}
    if mode == "c":
        return {"disable_initial_energy": "true", "surface": "false"}
    raise ValueError(f"Unsupported mode '{mode}'")


def _evaluate_checks(mode: str, log_text: str) -> Tuple[List[CheckResult], Dict[str, object]]:
    flags = _extract_mode_flags(log_text)
    em = _extract_global_em_settings(log_text)
    plist = _extract_process_list(log_text, "e-")
    plist_set = set(plist)

    mode_exp = _mode_expectations(mode)
    checks: List[CheckResult] = []

    def add(check_id: str, desc: str, ok: bool, observed: object, expected: object) -> None:
        checks.append(
            CheckResult(
                check_id=check_id,
                description=desc,
                passed=bool(ok),
                observed=str(observed),
                expected=str(expected),
            )
        )

    add(
        "C01",
        "Runtime strict mode flag is enabled",
        flags.get("fig8_strict", "").lower() == "true",
        flags.get("fig8_strict", ""),
        "true",
    )
    add(
        "C02",
        "LO phonon is disabled in baseline",
        flags.get("lo_phonon", "").lower() == "false",
        flags.get("lo_phonon", ""),
        "false",
    )
    add(
        "C03",
        "Capture is disabled in baseline",
        flags.get("capture", "").lower().startswith("false"),
        flags.get("capture", ""),
        "false",
    )
    add(
        "C04",
        "Mode-specific disable_initial_energy mapping",
        flags.get("disable_initial_energy", "").lower() == mode_exp["disable_initial_energy"],
        flags.get("disable_initial_energy", ""),
        mode_exp["disable_initial_energy"],
    )
    add(
        "C05",
        "Mode-specific surface mapping",
        flags.get("surface", "").lower() == mode_exp["surface"],
        flags.get("surface", ""),
        mode_exp["surface"],
    )
    add(
        "C06",
        "EM step limit type is fUseSafety (strict profile)",
        em.get("MscStepLimitType (e-/e+)", "") == "fUseSafety",
        em.get("MscStepLimitType (e-/e+)", ""),
        "fUseSafety",
    )
    add(
        "C07",
        "Lowest electron energy is 0 eV in benchmark regime",
        em.get("LowestElectronEnergy", "").startswith("0"),
        em.get("LowestElectronEnergy", ""),
        "0 eV",
    )
    add(
        "C08",
        "MicroElec inelastic process is present",
        "e-_G4Dielectrics" in plist_set,
        int("e-_G4Dielectrics" in plist_set),
        1,
    )
    add(
        "C09",
        "MicroElec elastic process is present (Si target)",
        "e-_G4MicroElecElastic" in plist_set,
        int("e-_G4MicroElecElastic" in plist_set),
        1,
    )

    # Surface process required for a/b and forbidden for c.
    has_surface_proc = "e-_G4MicroElecSurface" in plist_set
    expected_surface_proc = (mode != "c")
    add(
        "C10",
        "Surface process presence matches mode",
        has_surface_proc == expected_surface_proc,
        int(has_surface_proc),
        int(expected_surface_proc),
    )

    # Weakly bound override line expected in b/c only.
    has_override_line = "WeaklyBoundInitialEnergy override active" in log_text
    expected_override_line = (mode in {"b", "c"})
    add(
        "C11",
        "WeaklyBoundInitialEnergy override line appears only for b/c",
        has_override_line == expected_override_line,
        int(has_override_line),
        int(expected_override_line),
    )

    prohibited_present = sorted([p for p in PROHIBITED_PROCESSES if p in plist_set])
    add(
        "C12",
        "Strict mode removed prohibited standard processes from e- stack",
        len(prohibited_present) == 0,
        ", ".join(prohibited_present) if prohibited_present else "none",
        "none",
    )

    parsed = {
        "mode_flags": flags,
        "em_settings": em,
        "processes_electron": plist,
    }
    return checks, parsed


def main() -> None:
    parser = argparse.ArgumentParser(description="Probe Si Fig.8 runtime compliance for baseline modes a/b/c.")
    parser.add_argument("--exe", default="build/SEE_in_vacuum")
    parser.add_argument("--conda-env", default="geant4")
    parser.add_argument("--events", type=int, default=2000)
    parser.add_argument("--energy-ev", type=float, default=300.0)
    parser.add_argument(
        "--out-dir",
        default="results/si_fig8_runtime_compliance_probe_300eV_events2000",
    )
    args = parser.parse_args()

    project_root = Path.cwd()
    out_dir = project_root / args.out_dir
    cfg_dir = out_dir / "configs"
    log_dir = out_dir / "logs"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Tuple[str, CheckResult]] = []
    mode_summary: Dict[str, Dict[str, object]] = {}

    for mode in ("a", "b", "c"):
        base_cfg = project_root / BASELINE_CONFIGS[mode]
        if not base_cfg.exists():
            raise FileNotFoundError(f"Missing baseline config: {base_cfg}")

        probe_cfg = _build_probe_config(
            base_cfg_path=base_cfg,
            mode=mode,
            energy_ev=args.energy_ev,
            events=args.events,
            out_dir_tag=Path(args.out_dir).name,
        )
        probe_cfg_path = cfg_dir / f"probe_mode_{mode}_{int(round(args.energy_ev))}eV_events{args.events}.json"
        probe_cfg_path.write_text(json.dumps(probe_cfg, indent=2) + "\n")

        log_path = log_dir / f"probe_mode_{mode}_{int(round(args.energy_ev))}eV_events{args.events}.log"
        cmd = [
            "conda",
            "run",
            "-n",
            args.conda_env,
            str(project_root / args.exe),
            str(probe_cfg_path),
        ]
        _run_command(cmd, log_path)
        log_text = log_path.read_text(errors="ignore")

        checks, parsed = _evaluate_checks(mode, log_text)
        for check in checks:
            all_rows.append((mode, check))

        mode_summary[mode] = {
            "config": str(probe_cfg_path),
            "log": str(log_path),
            "n_checks": len(checks),
            "n_pass": sum(1 for c in checks if c.passed),
            "flags": parsed["mode_flags"],
        }

    # CSV output
    csv_path = out_dir / "si_fig8_runtime_compliance_checks.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mode", "check_id", "description", "passed", "observed", "expected"])
        for mode, check in all_rows:
            w.writerow(
                [
                    mode,
                    check.check_id,
                    check.description,
                    int(check.passed),
                    check.observed,
                    check.expected,
                ]
            )

    # Markdown summary
    md_path = out_dir / "si_fig8_runtime_compliance_checks.md"
    with md_path.open("w") as f:
        f.write("# Si Fig.8 Runtime Compliance Probe\n\n")
        f.write(f"- Probe energy: `{args.energy_ev:.0f} eV`\n")
        f.write(f"- Events per mode: `{args.events}`\n")
        f.write(f"- Executable: `{project_root / args.exe}`\n")
        f.write(f"- Conda env: `{args.conda_env}`\n\n")

        f.write("## Mode summaries\n\n")
        for mode in ("a", "b", "c"):
            s = mode_summary[mode]
            f.write(
                f"- mode `{mode}`: `{s['n_pass']}/{s['n_checks']}` checks passed"
                f" | config `{s['config']}` | log `{s['log']}`\n"
            )
        f.write("\n")

        f.write("## Detailed checks\n\n")
        f.write("| mode | check | description | pass | observed | expected |\n")
        f.write("|---|---|---|---:|---|---|\n")
        for mode, check in all_rows:
            f.write(
                f"| {mode} | {check.check_id} | {check.description} | "
                f"{'1' if check.passed else '0'} | {check.observed} | {check.expected} |\n"
            )

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
