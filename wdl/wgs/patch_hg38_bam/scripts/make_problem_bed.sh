#!/bin/bash
#
# Downloads the FixItFelix duplicated/collapsed hg38 problem-region BEDs,
# pads each interval by 100kb on each side, and writes the merged result to
# data/wgs/giab_fix_intervals_100k_padded.bed.

set -euo pipefail

pad=100000
outdir="data/wgs"
out="${outdir}/giab_fix_intervals_100k_padded.bed"

mkdir -p "${outdir}"

tmpdir=$(mktemp -d)
trap '/bin/rm -rf "${tmpdir}"' EXIT

curl -sSL https://github.com/srbehera/FixItFelix/raw/refs/heads/main/duplicated.bed \
    | tr -s ' ' '\t' | cut -f1-3 | grep -v '_random' > "${tmpdir}/duplicated.bed"
curl -sSL https://github.com/srbehera/FixItFelix/raw/refs/heads/main/collapsed.bed \
    | tr -s ' ' '\t' | cut -f1-3 | grep -v '_random' > "${tmpdir}/collapsed.bed"
curl -sSL -o "${tmpdir}/hg38.chrom.sizes" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

bedtools slop -i "${tmpdir}/duplicated.bed" -g "${tmpdir}/hg38.chrom.sizes" -b "${pad}" \
    > "${tmpdir}/duplicated.padded.bed"
bedtools slop -i "${tmpdir}/collapsed.bed" -g "${tmpdir}/hg38.chrom.sizes" -b "${pad}" \
    > "${tmpdir}/collapsed.padded.bed"

cat "${tmpdir}/duplicated.padded.bed" "${tmpdir}/collapsed.padded.bed" \
    | bedtools sort -g "${tmpdir}/hg38.chrom.sizes" -i - > "${out}"
