# Telebuilder

A toolkit for building TE annotations

## Installing

### Create new mamba/conda environment (recommended)

```bash
mamba env create -n telebuilder https://github.com/mlbendall/telebuilder/raw/main/environment.yml
````

### Install using pip

```bash
pip install git+https://github.com/mlbendall/telebuilder.git@main#egg=telebuilder
```

## Tools

### `buildERV`

`buildERV` assembles ERV proviral loci - ERV insertions with recognizable
internal regions - by combining genomic regions matching ERV internal regions
with flanking LTR sequences. The initial genomic regions, or "hits", along with 
the associated alignment information are retrieved from the UCSC Genome Browser
database (`rmsk` table). The main arguments include the internal model name
and the LTR model names. 

##### Example: HERVH

HERVH elements have internal regions matching the "HERVH-int" model and
flanking LTRs matching either "LTR7", "LTR7A", "LTR7B", "LTR7C", or "LTR7Y".
An example command to build HERVH would be:

```bash
buildERV\
  --genome_build hg38\
  --chrom_sizes chrom.sizes\
  --cytoband cytoband.gtf\
  HERVH HERVH-int LTR7,LTR7A,LTR7B,LTR7C,LTR7Y
```

where "hg38" is the genome build, "chrom.sizes" is a file with the chromosome 
names and sizes, and "cytoband.gtf" is a GTF file with cytogenetic band
regions.
