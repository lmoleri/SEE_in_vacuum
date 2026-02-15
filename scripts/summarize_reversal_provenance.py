#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter, defaultdict

import ROOT


CLASS_LABELS = {
    1: "entrance_exit",
    2: "opposite_exit",
    3: "lateral_exit",
    4: "stop_no_valid_exit",
}


def _find_root_files(results_dir):
    files = []
    for root, _, names in os.walk(results_dir):
        for name in names:
            if not name.endswith(".root"):
                continue
            if "_SEY_MonteCarlo" in name or name.startswith("summary"):
                continue
            files.append(os.path.join(root, name))
    return sorted(files)


def _meta_string(value):
    if value is None:
        return ""
    try:
        if isinstance(value, (bytes, bytearray)):
            raw = bytes(value)
        elif value.__class__.__name__ == "LowLevelView":
            raw = bytes(value)
        else:
            return str(value).split("\x00", 1)[0].strip()
        return raw.split(b"\x00", 1)[0].decode("utf-8", errors="ignore").strip()
    except Exception:
        return str(value).split("\x00", 1)[0].strip()


def _step_status_name(code):
    # G4StepStatus enum values
    mapping = {
        0: "fWorldBoundary",
        1: "fGeomBoundary",
        2: "fAtRestDoItProc",
        3: "fAlongStepDoItProc",
        4: "fPostStepDoItProc",
        5: "fUserDefinedLimit",
        6: "fExclusivelyForcedProc",
        7: "fUndefined",
    }
    return mapping.get(int(code), f"status_{int(code)}")


def _safe_div(a, b):
    if b == 0:
        return 0.0
    return float(a) / float(b)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize first-reversal provenance and step-limitation audit from EventDiagnostics."
    )
    parser.add_argument("--results-dir", required=True, help="Directory with ROOT outputs.")
    parser.add_argument("--output-csv", required=True, help="Output CSV path.")
    parser.add_argument(
        "--classes",
        default="2,4",
        help="Comma-separated primaryExitClass values to summarize (default: 2,4).",
    )
    args = parser.parse_args()

    classes = [int(x.strip()) for x in args.classes.split(",") if x.strip()]
    files = _find_root_files(args.results_dir)
    if not files:
        raise SystemExit("No ROOT files found.")

    rows = []
    for path in files:
        f = ROOT.TFile.Open(path)
        if not f or not f.IsOpen():
            continue
        meta = f.Get("RunMeta")
        tree = f.Get("EventDiagnostics")
        if not meta or not tree:
            f.Close()
            continue
        meta.GetEntry(0)
        step_nm = float(meta.maxStepNm)
        model = _meta_string(meta.emModel)
        e0_eV = float(meta.primaryEnergyMeV) * 1.0e6
        n_events = int(meta.primaryElectrons)

        per_class = {
            cls: {
                "N": 0,
                "first_rev_count": 0,
                "proc_counter": Counter(),
                "status_counter": Counter(),
                "delta_theta_sum": 0.0,
                "steplen_sum": 0.0,
                "preE_sum": 0.0,
                "postE_sum": 0.0,
                "audit_sum": defaultdict(float),
            }
            for cls in classes
        }

        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            cls = int(tree.primaryExitClass)
            if cls not in per_class:
                continue
            d = per_class[cls]
            d["N"] += 1

            # Step-audit counters
            d["audit_sum"]["nStepStatusGeomBoundary"] += float(tree.nStepStatusGeomBoundary)
            d["audit_sum"]["nStepStatusPostStepProc"] += float(tree.nStepStatusPostStepProc)
            d["audit_sum"]["nStepStatusAlongStepProc"] += float(tree.nStepStatusAlongStepProc)
            d["audit_sum"]["nStepStatusUserLimit"] += float(tree.nStepStatusUserLimit)
            d["audit_sum"]["nStepStatusOther"] += float(tree.nStepStatusOther)
            d["audit_sum"]["nProcMsc"] += float(tree.nProcMsc)
            d["audit_sum"]["nProcStepLimiter"] += float(tree.nProcStepLimiter)
            d["audit_sum"]["nProcTransportation"] += float(tree.nProcTransportation)
            d["audit_sum"]["nProcEIoni"] += float(tree.nProcEIoni)
            d["audit_sum"]["nProcOther"] += float(tree.nProcOther)

            if int(tree.firstDirectionReversalStep) > 0:
                d["first_rev_count"] += 1
                proc = _meta_string(tree.firstDirectionReversalProcess)
                status_name = _step_status_name(int(tree.firstDirectionReversalStepStatus))
                d["proc_counter"][proc] += 1
                d["status_counter"][status_name] += 1
                d["delta_theta_sum"] += float(tree.firstDirectionReversalDeltaThetaDeg)
                d["steplen_sum"] += float(tree.firstDirectionReversalStepLenNm)
                d["preE_sum"] += float(tree.firstDirectionReversalPreEnergyEv)
                d["postE_sum"] += float(tree.firstDirectionReversalPostEnergyEv)

        for cls in classes:
            d = per_class[cls]
            n = d["N"]
            nrev = d["first_rev_count"]
            top_proc, top_proc_n = ("", 0)
            if d["proc_counter"]:
                top_proc, top_proc_n = d["proc_counter"].most_common(1)[0]
            top_status, top_status_n = ("", 0)
            if d["status_counter"]:
                top_status, top_status_n = d["status_counter"].most_common(1)[0]

            row = {
                "root_file": os.path.basename(path),
                "step_nm": step_nm,
                "em_model": model,
                "E0_eV": e0_eV,
                "events_total": n_events,
                "exit_class": cls,
                "exit_class_label": CLASS_LABELS.get(cls, str(cls)),
                "N_class": n,
                "f_class": _safe_div(n, n_events),
                "N_first_reversal": nrev,
                "f_first_reversal_in_class": _safe_div(nrev, n),
                "top_first_reversal_process": top_proc,
                "top_first_reversal_process_frac": _safe_div(top_proc_n, nrev),
                "top_first_reversal_step_status": top_status,
                "top_first_reversal_step_status_frac": _safe_div(top_status_n, nrev),
                "mean_first_reversal_deltaTheta_deg": _safe_div(d["delta_theta_sum"], nrev),
                "mean_first_reversal_stepLen_nm": _safe_div(d["steplen_sum"], nrev),
                "mean_first_reversal_preE_eV": _safe_div(d["preE_sum"], nrev),
                "mean_first_reversal_postE_eV": _safe_div(d["postE_sum"], nrev),
            }

            for key in [
                "nStepStatusGeomBoundary",
                "nStepStatusPostStepProc",
                "nStepStatusAlongStepProc",
                "nStepStatusUserLimit",
                "nStepStatusOther",
                "nProcMsc",
                "nProcStepLimiter",
                "nProcTransportation",
                "nProcEIoni",
                "nProcOther",
            ]:
                row[f"mean_{key}"] = _safe_div(d["audit_sum"][key], n)

            rows.append(row)
        f.Close()

    if not rows:
        raise SystemExit("No rows produced. Check that EventDiagnostics has the new branches.")

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(args.output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["step_nm"], r["exit_class"])):
            writer.writerow(row)

    print(f"Saved CSV: {args.output_csv}")


if __name__ == "__main__":
    main()
