#!/bin/bash
#
# Downloads the FixItFelix duplicated/collapsed hg38 problem-region BEDs,
# pads each interval by 100kb on each side, merges overlapping intervals,
# and annotates each resulting interval with Hugo gene symbols.
#
# Usage:
#   ./make_giab_fix_intervals.sh /path/to/gencode.annotation.gtf.gz

set -euo pipefail

GTF="$1"
OUT_BED="$2"
pad=100000

tmpdir=$(mktemp -d)
trap '/bin/rm -rf "${tmpdir}"' EXIT

# Only retain canonical chromosomes.
chrom_filter='
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/
'
curl -sSL https://github.com/srbehera/FixItFelix/raw/refs/heads/main/duplicated.bed |
    tr -s ' ' '\t' |
    cut -f1-3 |
    awk "${chrom_filter}" \
    > "${tmpdir}/duplicated.bed"

curl -sSL https://github.com/srbehera/FixItFelix/raw/refs/heads/main/collapsed.bed |
    tr -s ' ' '\t' |
    cut -f1-3 |
    awk "${chrom_filter}" \
    > "${tmpdir}/collapsed.bed"

curl -sSL -o "${tmpdir}/hg38.chrom.sizes" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

awk -v OFS='\t' '
    BEGIN {
        for (i = 1; i <= 22; i++)
            wanted["chr" i] = i
        wanted["chrX"] = 23
        wanted["chrY"] = 24
    }

    $1 in wanted {
        size[$1] = $2
    }

    END {
        for (i = 1; i <= 22; i++)
            print "chr" i, size["chr" i]

        print "chrX", size["chrX"]
        print "chrY", size["chrY"]
    }
' "${tmpdir}/hg38.chrom.sizes" \
    > "${tmpdir}/canonical.chrom.sizes"

# Pad and merge the FixItFelix intervals.
cat "${tmpdir}/duplicated.bed" "${tmpdir}/collapsed.bed" |
    bedtools slop \
        -g "${tmpdir}/canonical.chrom.sizes" \
        -b "${pad}" \
        -i - |
    bedtools sort \
        -g "${tmpdir}/canonical.chrom.sizes" \
        -i - |
    bedtools merge \
        -i - \
    > "${tmpdir}/regions.bed"

# Extract gene features and Hugo symbols from the GTF.
zcat -f "${GTF}" |
    awk -F'\t' '
        $0 !~ /^#/ &&
        $3 == "gene" &&
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ &&
        $9 ~ /gene_type "protein_coding"/ {
            gene_name = ""

            if (match($9, /gene_name "[^"]+"/))
                gene_name = substr($9, RSTART + 11, RLENGTH - 12)

            if (gene_name != "")
                print $1 "\t" $4 - 1 "\t" $5 "\t" gene_name
        }
    ' |
    bedtools sort \
        -g "${tmpdir}/canonical.chrom.sizes" \
        -i - \
    > "${tmpdir}/genes.bed"

# Find overlapping genes and collect their symbols by interval.
bedtools intersect \
    -a "${tmpdir}/regions.bed" \
    -b "${tmpdir}/genes.bed" \
    -wa -wb |
    awk -v OFS='\t' '
        {
            key = $1 OFS $2 OFS $3
            genes[key SUBSEP $7] = 1
        }

        END {
            for (item in genes) {
                split(item, parts, SUBSEP)
                key = parts[1]
                gene = parts[2]

                if (!(key in region_genes))
                    region_genes[key] = gene
                else
                    region_genes[key] = region_genes[key] "," gene
            }

            for (key in region_genes) {
                split(key, region, OFS)
                print region[1], region[2], region[3], region_genes[key]
            }
        }
    ' |
    bedtools sort \
        -g "${tmpdir}/canonical.chrom.sizes" \
        -i - \
    > "${OUT_BED}"

echo "Wrote ${OUT_BED}"

echo "Markdown for README:"
{
    printf '| Chromosome | Start | End | Genes |\n'
    printf '|---|---:|---:|---|\n'
    awk -F'\t' -v OFS='|' '
        {
            gsub(/,/, "<br>", $4)
            print "|"$1, $2, $3, $4"|"
        }
    ' "${OUT_BED}"
}
