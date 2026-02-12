#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
from collections import defaultdict

import numpy as np
import uproot


QOI_COLUMNS = [
    "step_nm",
    "E0_eV",
    "f_stop",
    "f_ent",
    "f_opp",
    "mean_edep",
    "mode_edep",
    "q50_eexit_opp",
    "q90_eexit_opp",
    "sey",
]

METRIC_ORDER = [
    "f_stop",
    "f_ent",
    "f_opp",
    "mean_edep",
    "mode_edep",
    "q50_eexit_opp",
    "q90_eexit_opp",
    "sey",
]

FRACTION_METRICS = {"f_stop", "f_ent", "f_opp"}
ENERGY_METRICS = {"mean_edep", "mode_edep", "q50_eexit_opp", "q90_eexit_opp"}


def _scalar(value):
    if isinstance(value, np.ndarray):
        if len(value) == 0:
            return np.nan
        value = value[0]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _string_value(value):
    value = _scalar(value)
    if isinstance(value, bytes):
        return value.decode(errors="ignore").strip()
    return str(value).strip()


def _find_root_files(inputs):
    paths = []
    for token in inputs:
        expanded = glob.glob(token)
        if not expanded:
            expanded = [token]
        for path in expanded:
            if os.path.isdir(path):
                for name in sorted(os.listdir(path)):
                    if not name.endswith(".root"):
                        continue
                    if "_SEY_MonteCarlo" in name or name.startswith("summary"):
                        continue
                    paths.append(os.path.join(path, name))
            elif os.path.isfile(path) and path.endswith(".root"):
                paths.append(path)
            else:
                raise SystemExit(f"Input not found or unsupported: {path}")

    uniq = sorted(set(paths))
    if not uniq:
        raise SystemExit("No ROOT files found from the provided inputs.")
    return uniq


def _hist_stats(values, edges):
    total = float(np.sum(values))
    if total <= 0.0:
        return np.nan, np.nan
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean = float(np.sum(centers * values) / total)
    mode = float(centers[int(np.argmax(values))])
    return mean, mode


def _hist_quantile(values, edges, quantile):
    total = float(np.sum(values))
    if total <= 0.0:
        return np.nan

    target = quantile * total
    cum = 0.0
    for i, count in enumerate(values):
        next_cum = cum + float(count)
        if next_cum >= target:
            low = float(edges[i])
            high = float(edges[i + 1])
            width = high - low
            if count <= 0.0 or width <= 0.0:
                return low
            frac = (target - cum) / float(count)
            frac = min(max(frac, 0.0), 1.0)
            return low + frac * width
        cum = next_cum
    return float(edges[-1])


def _safe_get_hist(root_file, name):
    obj = root_file.get(name)
    if obj is None:
        raise SystemExit(f"Missing histogram '{name}' in {root_file.file_path}")
    return obj


def _load_record(path):
    f = uproot.open(path)
    run_meta = f.get("RunMeta")
    if run_meta is None:
        raise SystemExit(f"Missing RunMeta in {path}")
    meta = run_meta.arrays(library="np")

    energy_eV = float(_scalar(meta.get("primaryEnergyMeV", np.array([np.nan])))) * 1.0e6
    step_nm = float(_scalar(meta.get("maxStepNm", np.array([np.nan]))))
    n_events = float(_scalar(meta.get("primaryElectrons", np.array([np.nan]))))
    sey = float(_scalar(meta.get("sey", np.array([np.nan]))))

    if not np.isfinite(energy_eV) or not np.isfinite(step_nm) or not np.isfinite(n_events):
        raise SystemExit(f"Missing required RunMeta fields in {path}")
    if n_events <= 0.0:
        raise SystemExit(f"Invalid primaryElectrons={n_events} in {path}")

    h_exit = _safe_get_hist(f, "PrimaryExitClass")
    exit_values, _ = h_exit.to_numpy(flow=False)
    if len(exit_values) < 4:
        raise SystemExit(
            f"PrimaryExitClass in {path} has {len(exit_values)} bins; expected at least 4."
        )
    n_ent = float(exit_values[1])  # class=1 in ROOT bin 2
    n_opp = float(exit_values[2])  # class=2 in ROOT bin 3
    n_lat = float(exit_values[3])  # class=3 in ROOT bin 4
    n_stop = max(0.0, n_events - (n_ent + n_opp + n_lat))

    f_ent = n_ent / n_events
    f_opp = n_opp / n_events
    f_stop = n_stop / n_events

    h_edep = _safe_get_hist(f, "EdepPrimary")
    edep_values, edep_edges = h_edep.to_numpy(flow=False)
    mean_edep, mode_edep = _hist_stats(edep_values, edep_edges)

    h_opp_e = _safe_get_hist(f, "PrimaryExitEnergyOpposite")
    opp_values, opp_edges = h_opp_e.to_numpy(flow=False)
    q50 = _hist_quantile(opp_values, opp_edges, 0.50)
    q90 = _hist_quantile(opp_values, opp_edges, 0.90)

    return {
        "file": path,
        "step_nm": step_nm,
        "E0_eV": energy_eV,
        "f_stop": f_stop,
        "f_ent": f_ent,
        "f_opp": f_opp,
        "mean_edep": mean_edep,
        "mode_edep": mode_edep,
        "q50_eexit_opp": q50,
        "q90_eexit_opp": q90,
        "sey": sey,
        "primary_particle": _string_value(meta.get("primaryParticle", np.array([""]))),
        "em_model": _string_value(meta.get("emModel", np.array([""]))),
        "thickness_nm": float(_scalar(meta.get("sampleThicknessNm", np.array([np.nan])))),
        "substrate_nm": float(_scalar(meta.get("substrateThicknessNm", np.array([np.nan])))),
        "radius_nm": float(_scalar(meta.get("sampleRadiusNm", np.array([np.nan])))),
    }


def _rel_delta(coarse, fine, eps):
    if not np.isfinite(coarse) or not np.isfinite(fine):
        return np.nan
    denom = max(abs(fine), eps)
    return abs(coarse - fine) / denom


def _abs_delta(coarse, fine):
    if not np.isfinite(coarse) or not np.isfinite(fine):
        return np.nan
    return abs(coarse - fine)


def _step_pairs(steps):
    pairs = []
    for coarse in steps:
        target = coarse / 2.0
        fine = None
        for candidate in steps:
            if math.isclose(candidate, target, rel_tol=1.0e-8, abs_tol=1.0e-10):
                fine = candidate
                break
        if fine is not None:
            pairs.append((coarse, fine))
    return pairs


def _to_csv_value(value):
    if isinstance(value, float) and (not np.isfinite(value)):
        return ""
    return value


def _write_csv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: _to_csv_value(row.get(col, "")) for col in columns})


def main():
    parser = argparse.ArgumentParser(
        description="Extract step-convergence QoIs and pairwise pass/fail matrix from ROOT outputs."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input ROOT files, result directories, and/or glob patterns.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where CSV outputs are written.",
    )
    parser.add_argument(
        "--label",
        default="step_convergence",
        help="Filename prefix for output CSV files.",
    )
    parser.add_argument(
        "--energy-min-ev",
        type=float,
        default=350.0,
        help="Minimum primary energy (eV) included in the convergence table.",
    )
    parser.add_argument(
        "--energy-max-ev",
        type=float,
        default=900.0,
        help="Maximum primary energy (eV) included in the convergence table.",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1.0e-6,
        help="Epsilon used in relative-delta denominators.",
    )
    parser.add_argument(
        "--fraction-threshold",
        type=float,
        default=0.02,
        help="Relative threshold for f_stop/f_ent/f_opp pass-fail.",
    )
    parser.add_argument(
        "--sey-threshold",
        type=float,
        default=0.02,
        help="Relative threshold for SEY pass-fail.",
    )
    parser.add_argument(
        "--energy-threshold-ev",
        type=float,
        default=10.0,
        help="Absolute threshold for mean/mode/quantile energy metrics.",
    )
    args = parser.parse_args()

    root_files = _find_root_files(args.inputs)
    all_records = [_load_record(path) for path in root_files]

    records = [
        rec for rec in all_records
        if rec["E0_eV"] >= args.energy_min_ev and rec["E0_eV"] <= args.energy_max_ev
    ]
    if not records:
        raise SystemExit(
            "No records in the requested energy window "
            f"[{args.energy_min_ev}, {args.energy_max_ev}] eV."
        )

    em_models = sorted({rec["em_model"] for rec in records if rec["em_model"]})
    particles = sorted({rec["primary_particle"] for rec in records if rec["primary_particle"]})
    if len(em_models) > 1:
        raise SystemExit(
            f"Found mixed EM models in inputs ({', '.join(em_models)}). "
            "Run one campaign/model at a time."
        )
    if len(particles) > 1:
        raise SystemExit(
            f"Found mixed primary particles in inputs ({', '.join(particles)}). "
            "Run one campaign/particle at a time."
        )

    key_map = {}
    for rec in records:
        key = (round(rec["step_nm"], 12), round(rec["E0_eV"], 6))
        if key in key_map:
            raise SystemExit(
                "Duplicate (step_nm, E0_eV) point found:\n"
                f"  {key_map[key]['file']}\n  {rec['file']}"
            )
        key_map[key] = rec

    qoi_rows = []
    for rec in sorted(records, key=lambda r: (-r["step_nm"], r["E0_eV"])):
        qoi_rows.append({col: rec[col] for col in QOI_COLUMNS})

    # Build pairwise rows and pass/fail matrix.
    per_step = defaultdict(dict)
    for rec in records:
        per_step[round(rec["step_nm"], 12)][round(rec["E0_eV"], 6)] = rec
    unique_steps = sorted((rec["step_nm"] for rec in records), reverse=True)
    step_pairs = _step_pairs(unique_steps)

    pairwise_rows = []
    matrix_rows = []

    for coarse, fine in step_pairs:
        coarse_key = round(coarse, 12)
        fine_key = round(fine, 12)
        common_energies = sorted(set(per_step[coarse_key]).intersection(per_step[fine_key]))
        if not common_energies:
            continue

        metric_status = {metric: True for metric in METRIC_ORDER}
        metric_counts = {metric: 0 for metric in METRIC_ORDER}
        metric_pass_counts = {metric: 0 for metric in METRIC_ORDER}

        for e_key in common_energies:
            coarse_rec = per_step[coarse_key][e_key]
            fine_rec = per_step[fine_key][e_key]
            for metric in METRIC_ORDER:
                x_coarse = float(coarse_rec[metric])
                x_fine = float(fine_rec[metric])

                if metric in FRACTION_METRICS:
                    delta_type = "relative"
                    threshold = args.fraction_threshold
                    delta = _rel_delta(x_coarse, x_fine, args.epsilon)
                elif metric in ENERGY_METRICS:
                    delta_type = "absolute"
                    threshold = args.energy_threshold_ev
                    delta = _abs_delta(x_coarse, x_fine)
                else:
                    delta_type = "relative"
                    threshold = args.sey_threshold
                    delta = _rel_delta(x_coarse, x_fine, args.epsilon)

                passed = bool(np.isfinite(delta) and delta < threshold)

                pairwise_rows.append(
                    {
                        "step_coarse_nm": coarse,
                        "step_fine_nm": fine,
                        "E0_eV": coarse_rec["E0_eV"],
                        "metric": metric,
                        "coarse_value": x_coarse,
                        "fine_value": x_fine,
                        "delta": delta,
                        "delta_type": delta_type,
                        "threshold": threshold,
                        "pass": int(passed),
                    }
                )

                metric_counts[metric] += 1
                if passed:
                    metric_pass_counts[metric] += 1
                else:
                    metric_status[metric] = False

        matrix_row = {
            "step_coarse_nm": coarse,
            "step_fine_nm": fine,
            "n_energies": len(common_energies),
        }
        for metric in METRIC_ORDER:
            matrix_row[f"{metric}_pass"] = int(metric_status[metric])
            matrix_row[f"{metric}_n_pass"] = metric_pass_counts[metric]
            matrix_row[f"{metric}_n_total"] = metric_counts[metric]
        matrix_row["all_qois_pass"] = int(all(matrix_row[f"{m}_pass"] == 1 for m in METRIC_ORDER))
        matrix_rows.append(matrix_row)

    os.makedirs(args.output_dir, exist_ok=True)
    prefix = os.path.join(args.output_dir, args.label)

    qoi_csv = f"{prefix}_qoi_table.csv"
    pairwise_csv = f"{prefix}_pairwise_deltas_passfail.csv"
    matrix_csv = f"{prefix}_pairwise_passfail_matrix.csv"

    _write_csv(qoi_csv, qoi_rows, QOI_COLUMNS)
    _write_csv(
        pairwise_csv,
        sorted(
            pairwise_rows,
            key=lambda r: (-r["step_coarse_nm"], r["E0_eV"], METRIC_ORDER.index(r["metric"])),
        ),
        [
            "step_coarse_nm",
            "step_fine_nm",
            "E0_eV",
            "metric",
            "coarse_value",
            "fine_value",
            "delta",
            "delta_type",
            "threshold",
            "pass",
        ],
    )
    _write_csv(
        matrix_csv,
        sorted(matrix_rows, key=lambda r: -r["step_coarse_nm"]),
        [
            "step_coarse_nm",
            "step_fine_nm",
            "n_energies",
            "f_stop_pass",
            "f_stop_n_pass",
            "f_stop_n_total",
            "f_ent_pass",
            "f_ent_n_pass",
            "f_ent_n_total",
            "f_opp_pass",
            "f_opp_n_pass",
            "f_opp_n_total",
            "mean_edep_pass",
            "mean_edep_n_pass",
            "mean_edep_n_total",
            "mode_edep_pass",
            "mode_edep_n_pass",
            "mode_edep_n_total",
            "q50_eexit_opp_pass",
            "q50_eexit_opp_n_pass",
            "q50_eexit_opp_n_total",
            "q90_eexit_opp_pass",
            "q90_eexit_opp_n_pass",
            "q90_eexit_opp_n_total",
            "sey_pass",
            "sey_n_pass",
            "sey_n_total",
            "all_qois_pass",
        ],
    )

    print(f"Wrote QoI table: {qoi_csv}")
    print(f"Wrote pairwise deltas/pass-fail: {pairwise_csv}")
    print(f"Wrote pairwise pass/fail matrix: {matrix_csv}")


if __name__ == "__main__":
    main()
