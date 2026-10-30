#!/usr/bin/env bash
# Build both documents with the cluster TeX Live module.
set -eo pipefail
source /usr/local/etc/modules.sh
module load texlive/20260310-GCC-13.2.0
# Prefer system zlib to avoid CPU-specific instructions in the module dependency.
if [[ -f /lib64/libz.so.1 ]]; then
    export LD_PRELOAD="/lib64/libz.so.1${LD_PRELOAD:+:$LD_PRELOAD}"
fi
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
for document in paper_pre entropy_theory; do
    for pass in 1 2 3; do
        echo "Building ${document}, pass ${pass}"
        pdflatex -interaction=nonstopmode -halt-on-error "${document}.tex" > "${document}.build.log" 2>&1
        # Stop once references and PDF outlines have stabilized.
        if ! grep -Eq 'undefined references|undefined citations|Citation.*undefined|Label\(s\) may have changed|Rerun to get' "${document}.log"; then
            break
        fi
    done
    if grep -Eq 'undefined references|undefined citations' "${document}.log"; then
        echo "Unresolved references in ${document}.log" >&2
        exit 1
    fi
    echo "Built ${document}.pdf"
    if grep -En 'undefined references|undefined citations|Rerun to get|Overfull' "${document}.log"; then
        echo "Review layout or reference warnings above."
    fi
done
