# Transcript annotation

Produces the portal `Transcript` table: per-transcript PFAM domain annotations, gene symbols,
UniProt protein isoform IDs, and coding/NMD/structural-category flags, for both the annotated
(Ensembl 104) and novel (IsoQuant) transcripts in a long-read RNA release.

This is a manual, run on a VM plus a local machine — it is not a WDL
workflow.

## Pipeline

```
                    combine_gtfs release outputs (GCS)
                                 |
                  novel_transcripts.fa (non-ENST records)
                                 |
   [1] IsoAnnot x2 on VM  --------+--------  Ensembl 104 reference
        ensembl          |                  mytranscripts
                         v
   [2] Isoannot_merge_domain_annotations_gitv.ipynb   (on VM)
                         |
              merged_annotations_all_git.csv
                         |  (copy down off the VM)
                         v
   [3] make_transcript_ann_for_portal_upload.ipynb    (local)
         + Ensembl 104 GTF
         + OmicsLongReadGTF        (Taiga)
         + hgnc_complete_set       (Taiga)
         + uniprot_sprot_varsplic  (UniProt FTP)
                         |
              Transcript_forportal_ptc.csv --> Taiga
```

### 1. Build the image and run IsoAnnot

On the VM. Full instructions, including where the config files must be installed, are in
[run_isoannot.md](run_isoannot.md).

```
docker build -t isoannot .
```

[Dockerfile](Dockerfile) builds an sshd image (for use with `hermit`) with IsoAnnot
`releases/0.9.0b1` at `/opt/IsoAnnot`, symlinked to `/home/ubuntu/IsoAnnot`, and a conda env
`isoannot_env` created from IsoAnnot's `snakemake.yaml`. Exec in over ssh as `ubuntu`; the env is
activated by `.bashrc`.

Both IsoAnnot runs are long (many hours) and download large references — run them under `screen`.

### 2. Merge the domain annotations

[Isoannot_merge_domain_annotations_gitv_depmapomics.ipynb](Isoannot_merge_domain_annotations_gitv_depmapomics.ipynb), run on
the VM (Jupyter is installed into `isoannot_env` by the Dockerfile).

Set `ENSEMBL_OUTDIR` / `NOVEL_OUTDIR` in the config cell to the `--outputdir` values used in step 1.
The cell asserts every expected IsoAnnot output exists before doing any work.

For each side it filters transcripts, then attaches NMD calls and the set of PFAM domain names per
transcript:

| | annotated (`ensembl`) | novel (`mytranscripts`) |
| --- | --- | --- |
| filter | `gene_biotype == protein_coding` and a translated peptide exists | SQANTI `ORF_seq` is not null |
| domains from | `human_tappas_ensembl_annotation_file.gff3_mod` | `layers/layer_interproscan.gtf` |
| NMD calls from | `sqanti_NMD.txt` | `nmd_data.txt` |
| `sequence` is | the translated peptide (`sequence_y`, from `*.pep.all.fa`) | SQANTI's predicted `ORF_seq` |

`sequence` is what the portal notebook joins UniProt on, so both sides rename their peptide column
to that name. Novel entries are *predicted* ORFs, so few match Swiss-Prot exactly — most novel rows
end up with a null `ProteinIsoformID`/`UniprotID`, which is expected.

The two sides are then reduced to their common columns, `chrom` is prefixed with `chr`, and they are
concatenated. Output: `merged_annotations_all_git.csv`.


### 3. Build the portal table and upload

[make_transcript_ann_for_portal_upload.ipynb](make_transcript_ann_for_portal_upload.ipynb), run
locally. Copy `merged_annotations_all_git.csv` down off the VM first.

Two files must be downloaded before running (commands are in the notebook's config section):

- `Homo_sapiens.GRCh38.104.chr.gtf` — Ensembl 104 transcript GTF. Must be the same release the
  annotated IsoAnnot run used.
- `uniprot_sprot_varsplic.fasta` — UniProt Swiss-Prot varsplic sequences.

Update the Taiga pins (`RELEASE_PERMANAME`, `RELEASE_VERSION`, `HGNC_NAME`, `HGNC_VERSION`) in the
config cell for the current quarter.

The notebook keeps only transcripts expressed in the release (present in `OmicsLongReadGTF`),
resolves gene symbols against HGNC — dropping Ensembl gene IDs that map to multiple symbols, and
retrying failures on `gene_name` — joins UniProt isoform IDs on exact peptide sequence plus symbol,
and writes `Transcript_forportal_ptc.csv`.

**The last cell uploads to the release dataset on Taiga.** 

## Files

| file | purpose |
| --- | --- |
| [Dockerfile](Dockerfile) | IsoAnnot 0.9.0b1 + conda env + Jupyter, sshd image |
| [run_isoannot.md](run_isoannot.md) | how to install the configs and run both IsoAnnot passes |
| [Isoannot_config_annotated_transcripts.yaml](Isoannot_config_annotated_transcripts.yaml) | config for the `ensembl` run; installed to `config/ensembl/hsapiens/config.yaml` |
| [Isoannot_config_novel_transcripts.yaml](Isoannot_config_novel_transcripts.yaml) | config for the `mytranscripts` run; installed to `config/mytranscripts/hsapiens/config.yaml` |
| [Isoannot_merge_domain_annotations_gitv.ipynb](Isoannot_merge_domain_annotations_gitv.ipynb) | step 2, on the VM |
| [make_transcript_ann_for_portal_upload.ipynb](make_transcript_ann_for_portal_upload.ipynb) | step 3, local |

