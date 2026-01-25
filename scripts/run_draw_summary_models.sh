#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <pai_dir> <livermore_dir> <penelope_dir> [summary_root] [root_macro]"
  echo "Example: $0 results/scan_*_modelPAI results/scan_*_modelLivermore results/scan_*_modelPenelope"
  exit 1
fi

pai_dir="$1"
liv_dir="$2"
pen_dir="$3"
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
out_file="${4:-$repo_dir/results/summary_models.root}"
macro="${5:-$repo_dir/draw_summary_models.C}"

if [[ ! -d "$pai_dir" ]]; then
  echo "Error: directory not found: $pai_dir"
  exit 1
fi
if [[ ! -d "$liv_dir" ]]; then
  echo "Error: directory not found: $liv_dir"
  exit 1
fi
if [[ ! -d "$pen_dir" ]]; then
  echo "Error: directory not found: $pen_dir"
  exit 1
fi

if [[ ! -f "$macro" ]]; then
  echo "Error: macro not found: $macro"
  exit 1
fi

root -l -b -q "$macro(\"$pai_dir\", \"$liv_dir\", \"$pen_dir\", \"$out_file\")"
