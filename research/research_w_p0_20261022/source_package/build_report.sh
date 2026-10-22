#!/usr/bin/env bash
set -eu
report_root="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
cd "$report_root/report_parts"
for report_pass in 1 2 3; do
    xelatex -interaction=nonstopmode -halt-on-error main.tex > "compile_${report_pass}.txt"
done
cp main.pdf "$report_root/Research_W_Report_ZH.pdf"
printf '%s\n' 'Built Research_W_Report_ZH.pdf'
