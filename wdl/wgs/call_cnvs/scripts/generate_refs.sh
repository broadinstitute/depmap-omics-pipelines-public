#!/bin/zsh

wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz

bedtools intersect -v -a hg38.1000.bed -b hg38-blacklist.v2.bed.gz \
        > hg38.1000.bl_excluded.2.bed
bedtools sort -i hg38.1000.bl_excluded.2.bed -g hg38.chrom.sizes \
        > hg38.1000.sorted.bed

# GC wig
bedtools nuc -fi Homo_sapiens_assembly38.fasta -bed hg38.1000.sorted.bed \
        > hg38.1000.gc.txt

cut -f1,2,3,5 hg38.1000.gc.txt \
        | sed '1d' \
                > hg38.1000.gc.bedgraph

bedGraphToBigWig hg38.1000.gc.bedgraph hg38.chrom.sizes hg38.1000.gc.bw
bigWigToWig hg38.1000.gc.bw hg38.1000.gc.wig

sed 's/#bedGraph section /fixedStep chrom=/g' hg38.1000.gc.wig \
        | sed 's/:/ start=/g' \
        | sed 's/\-[0-9]*/ step=1000 span=1000/g' \
        | awk '{if(index($1,"chr") == 1) {print $4 } else print $0}' \
                > hg38.1000.gc.formatted.wig

# Mappability wig
bwtool/bwtool summary -fill=0 hg38.1000.sorted.bed k100.Umap.MultiTrackMappability.bw stdout \
        | awk '{ print $1"\t"$2"\t"$3"\t"$6 }' \
                > hg38.1000.with_mean_mapscore.bed

bedGraphToBigWig hg38.1000.with_mean_mapscore.bed hg38.chrom.sizes hg38.1000.mapscore.bw

bigWigToWig hg38.1000.mapscore.bw hg38.1000.mapscore.wig

sed 's/#bedGraph section /fixedStep chrom=/g' hg38.1000.mapscore.wig \
        | sed 's/:/ start=/g' \
        | sed 's/\-[0-9]*/ step=1000 span=1000/g' \
        | awk '{if(index($1,"chr") == 1) {print $4 } else print $0}' \
                > hg38.1000.mapscore.formatted.wig

# Bin-level mappability with gene annotation
printf 'chr\tstart\tend\tmap\tgene\n' > hg38.1000.with_mean_mapscore.annotated_gene.tsv

bedtools intersect -a hg38.1000.with_mean_mapscore.bed \
  -b gencode.v38.primary_assembly.annotation.protein_coding_genes.bed -wao \
  | awk 'BEGIN{FS=OFS="\t"} $5!="." {print $1,$2,$3,$4,$8}' \
  >> hg38.1000.with_mean_mapscore.annotated_gene.tsv

# Gene-level, % of bins mappable
awk -F'\t' '
        NR == 1 {
                for (i = 1; i <= NF; i++) col[$i] = i
                next
        }
        {
                g = $col["gene"]
                n[g]++
                if ($col["map"] + 0 < 0.9) low[g]++
        }
        END {
                for (g in n) {
                        pct = (low[g] / n[g]) * 100
                        printf "%s,%s\n", g, (pct == int(pct) ? sprintf("%.1f", pct) : sprintf("%.10g", pct))
                }
        }
' hg38.1000.with_mean_mapscore.annotated_gene.tsv \
        | sort -t, -k1,1 \
        | cat <(printf 'gene,pct_bins_map_lt_90pct\n') - \
                > cn_pct_bins_mappable_lt_90pct.csv
