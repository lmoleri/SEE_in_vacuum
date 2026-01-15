#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <scan_output_dir> [root_macro]"
  echo "Example: $0 scan_thick5-10-15-20-25nm_energy1MeV_events100000"
  exit 1
fi

scan_dir="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
macro="${2:-$repo_dir/draw_histo.C}"

if [[ ! -d "$scan_dir" ]]; then
  echo "Error: directory not found: $scan_dir"
  exit 1
fi

if [[ ! -f "$macro" ]]; then
  echo "Error: macro not found: $macro"
  exit 1
fi

shopt -s nullglob
roots=("$scan_dir"/*.root)
if [[ ${#roots[@]} -eq 0 ]]; then
  echo "No .root files found in $scan_dir"
  exit 0
fi

for f in "${roots[@]}"; do
  echo "Processing $f"
  root -l -b -q "$macro(\"$f\")"
done
