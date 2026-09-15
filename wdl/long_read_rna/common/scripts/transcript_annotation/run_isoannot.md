# Running IsoAnnot

<https://github.com/ConesaLab/IsoAnnot>

IsoAnnot is run twice: once over the Ensembl 104 reference transcriptome (`--database ensembl`) and
once over the novel IsoQuant transcripts from the release (`--database mytranscripts`). See
[README.md](README.md) for where this sits in the overall pipeline.

The container is built from [Dockerfile](Dockerfile), which clones IsoAnnot `releases/0.9.0b1` to
`/opt/IsoAnnot` and symlinks it to `/home/ubuntu/IsoAnnot`.

## Installing the config files

`isoannot.sh` does not take a config path by default — it looks for a config at a fixed location
derived from `--database` and `--species`:

```
$ISOANNOT_DIR/config/${DATABASE}/${SPECIES}/config.yaml
```

where `$ISOANNOT_DIR` is the directory containing `isoannot.sh` (`/home/ubuntu/IsoAnnot` in the
container). IsoAnnot ships its own `config.yaml` at each of those paths, so the two configs in this
directory must **overwrite** the shipped ones before running:

```
cp Isoannot_config_annotated_transcripts.yaml \
  /home/ubuntu/IsoAnnot/config/ensembl/hsapiens/config.yaml

cp Isoannot_config_novel_transcripts.yaml \
  /home/ubuntu/IsoAnnot/config/mytranscripts/hsapiens/config.yaml
```

Alternatively, pass `--configfile <path>` to `isoannot.sh` to point at a config in place without
overwriting anything. `isoannot.sh` exits with "The snakefile or configfile requested do not exist"
if the resolved config is missing.

### How these configs differ from the ones IsoAnnot ships

Both files are derived from upstream's `config/{ensembl,mytranscripts}/hsapiens/config.yaml` with
two deliberate DepMap changes:

- **Ensembl release is pinned to 104** (upstream `releases/0.9.0b1` ships release 115). The
  `ensembl_cdna`, `ensembl_proteins`, `ensembl_gtf`, and `ensembl_reference` URLs are all pinned.
  This must match the GTF release used in `make_transcript_ann_for_portal_upload_depmapomics.ipynb`.
- **`transcript_versioned: False`** in the novel config (upstream ships `True`), because the
  IsoQuant transcript IDs we feed in are unversioned.

Note that the `uniprot_fasta` / `uniprot_dat` URLs use `ftp://` in the annotated config and
`https://` in the novel config. This is inherited verbatim from upstream's two config files and is
not a DepMap change.

## Running IsoAnnot for annotated transcripts

```
./isoannot.sh --database ensembl --species hsapiens --outputdir ./results_ensembl104
```

## Running IsoAnnot for novel transcripts

### Before running

Extract the novel (non-Ensembl) transcripts and their sequences from the release's combined
transcriptome. Substitute the current release's `combine_gtfs` output paths:

```
TRANSCRIPTOME_FA='gs://fc-secure-d068ecea-bdcf-4e38-a449-edcd8933a71c/submissions/final-outputs/9859061a-3979-468a-9a37-a55c853a25bd/combine_gtfs/1a26f9d7-6ee3-4d2b-8c24-0da75c499ded/call-filter_gtf_and_tracking/26q3_transcriptome.fa'
PROCESSED_GTF='gs://fc-secure-d068ecea-bdcf-4e38-a449-edcd8933a71c/submissions/final-outputs/9859061a-3979-468a-9a37-a55c853a25bd/combine_gtfs/1a26f9d7-6ee3-4d2b-8c24-0da75c499ded/call-process_gtf/26q3_processed.gtf'

gcloud storage cp "$TRANSCRIPTOME_FA" 26q3_transcriptome.fa
gcloud storage cp "$PROCESSED_GTF" 26q3_processed.gtf

# Keep only records whose ID does not start with ENST, i.e. the novel transcripts.
awk -v RS=">" -v ORS="" 'NR>1 && $1 !~ /^ENST/ {print ">" $0}' \
  26q3_transcriptome.fa > novel_transcripts.fa

# Novel transcript IDs, for reference / QC.
awk '$2=="IsoQuant"' 26q3_processed.gtf > isoquant_entries.gtf
awk -F'transcript_id "' '{print $2}' isoquant_entries.gtf \
  | awk -F'"' '!seen[$1]++{print $1}' > unique_transcript_ids.txt
```

### Running IsoAnnot

`fasta_cdna` is a Snakemake `--config` override pointing at the novel transcript FASTA produced
above; give it an absolute path.

```
./isoannot.sh --database mytranscripts --species hsapiens \
  --outputdir results_novel_transcripts_v2 \
  --config fasta_cdna=/home/ubuntu/my_novel_transcriptome_26Q3/novel_transcripts.fa
```

## Post-processing

The `--outputdir` values used above are the `ENSEMBL_OUTDIR` / `NOVEL_OUTDIR` settings in the
config cell of [Isoannot_merge_domain_annotations_gitv.ipynb](Isoannot_merge_domain_annotations_gitv.ipynb).
