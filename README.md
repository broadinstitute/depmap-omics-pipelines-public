depmap-omics-pipelines-public
---

This repo contains WDL used to process raw omics data (WGS, RNA, and long read RNA under [DMC](https://depmap.org/portal/home/#/depmap-consortium) embargo) for [depmap.org](https://depmap.org/). Each workflow consists of the following:

- [WDL 1.0](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) file (e.g. `call_cnvs.wdl`)
- JSON file of workflow inputs (`call_cnvs_inputs.json`)
- README (`README.md`)

Also available are:

- Dockerfiles needed to build images used in the WDL tasks
- other supporting files/scripts used in tasks
- custom Python modules
- scripts for (re)generating input reference files

See also the `common` folders, which contain supporting files shared between workflows.

## Running workflows

Workflows can be run on Terra or via [miniwdl](https://github.com/chanzuckerberg/miniwdl). They can be imported from Dockstore for URLs like this:

https://dockstore.org/workflows/github.com/broadinstitute/depmap-omics-pipelines-public/align_wgs_sample:main?tab=info

Note that tasks reference files stored on our Google Cloud Storage buckets, and Docker images that are stored on our private Artifact Registry repo. Some supporting files require software licenses for external databases, and these files can't be shared publicly. You'll need to override those URLs with your own files and images, which should ideally be colocated with your computing infrastructure. 
