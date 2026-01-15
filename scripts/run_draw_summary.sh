#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <scan_output_dir> [summary_root] [root_macro]"
  echo "Example: $0 results/scan_thick5-10-15-20-25nm_energy1MeV_events100000"
  exit 1
fi

scan_dir="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
out_file="${2:-$scan_dir/summary.root}"
macro="${3:-$repo_dir/draw_summary.C}"

if [[ ! -d "$scan_dir" ]]; then
  echo "Error: directory not found: $scan_dir"
  exit 1
fi

if [[ ! -f "$macro" ]]; then
  echo "Error: macro not found: $macro"
  exit 1
fi

root -l -b -q "$macro(\"$scan_dir\", \"$out_file\")"
