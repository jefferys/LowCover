# LowCover
Identifies low coverage in selected genomic regions (e.g. genes) based on mosdepth-like coverage info (e.g. genome wide coverage annotation).

## Installing the LowCover package

LowCover can be installed from GitHub within R if you have `BiocManager` and `remotes` installed

```R
> # install.packages("BiocManager")
> # install.packages("remotes")
> BiodManager::install("Jefferys/LowCover")
```

## Running LowCover

LowCover is an R package with a "CLI" wrapper function. This allows it to be used from within R, run  from a local script, or executed within a container.

LowCover uses the _future_ package for parallel processing. Most users who know nothing about the _future_ package will just want to set the option `--parallel guess`; this  makes LowCover run with (all available) cores using a forking process where that is supported or with one R session per core where not. Without this parameter or additional environmental setup, LowCover runs purely sequentially. If you know something about _future_, e.g. you want to use your own custom `future::plan()` to manage running R code, LowCover will let you; it will use whatever plan already exists when it runs, or in the normal way use the R options that future defines to choose a default. It also tries to be friendly if you do choose to use the `--parallel` option, restoring any pre-existing future plan contexts when it exits.

### Running from within R

The main function `LowCoverApp()` is designed to be called from a script, as it automatically loads the R script command line parameters internally and parses them. However, it can be run in an R session by passing command line options and values as a vector, e.g.:

```R
LowCoverApp( c( "--targetsFile", "targets.bed",
                "--coverageFile", "mozdepth.coverage.gz",
                "--parallel", "guess" ))
```

### Running locally from a script

Create a script (or copy it from the `Docker` directory)

**LowCover**

```sh
!/usr/bin/env Rscript

LowCover::LowCoverApp()
```

Running this is done from the shell command line like a normal program, generally passing `--targetsFile`, `--coverageFile`, and --parallel.

```sh
$ LowCover --targetFile "targets.bed" --coverageFile "mozdepth.coverage.gz" \
    --parallel "guess"
```

You can copy this script somewhere on your path to make it globally available.

### Running using a container

You can use a pre-existing docker container or build you own

#### Getting the docker container from quay.io

```shell
docker login quay.io
docker pull `quay.io/unclineberger/jefferys/lowcover:latest`
```

#### Building your own Docker container

Assumes you have docker, and that the contents of the Docker directory from the LowCover repository is available locally.

```shell
cd Docker
docker build  -t <local name>:<tag> .
```

This builds a container named `<local name>:<tag>`, e.g. `lowcover:latest`

#### Running with Docker

Running involves simply running LowCover in the Docker container. For example

```shell
docker run -it --rm -v `pwd`:`pwd`, -w `pwd` LowCover \
    --targetFile "targets.bed" \
    --coverageFile "mozdepth.coverage.gz" \
    --parallel "guess"
```

#### Running with Singularity

To run with singularity, convert to a singularity file and run that. You can either pull to a sif directly from quay.io

```shell
singularity login quay.io
singlualrity pull docker://quay.io/unclineberger/jefferys/lowcover:latest <local name>:<tag>
```

or convert a local container with

```shell
singlualrity pull docker-daemon://quay.io/unclineberger/jefferys/lowcover:latest <local name>:<tag>
```

The sif will be named `<local name>_<tag>.sif`

Then run it with singularity

```shell
singularity run --B `pwd`:`pwd`, --pwd `pwd` LowCover \
    --targetFile "targets.bed" \
    --coverageFile "mozdepth.coverage.gz" \
    --parallel "guess"
```

## LowCover options

Help is available, e.g. from the shell with

```shell
LowCover --help
```

or from R with

LowCover::LowCoverApp( "--help" )

Which currently returns

```
usage: LowCover [--] [--help] [--version] [--force] [--targetsFile TARGETSFILE]
       [--coverageFile COVERAGEFILE] [--keep KEEP] [--chrY CHRY] [--regionsFile
       REGIONSFILE] [--badGenesFile BADGENESFILE] [--goodGenesFile GOODGENESFILE]
       [--summaryFile SUMMARYFILE] [--summaryFileNoY SUMMARYFILENOY] [--parallel
       PARALLEL] [--workers WORKERS]

Low coverage reporting

Reports low coverage genomic regions based on mosdepth.

flags:
  -h, --help            show this help message and exit
  -V, --version         Print version info to stderr and exit.
  -f, --force           Create any non-existing directories and overwrite any
                        existing files.

optional arguments:
  -t, --targetsFile     Bed4 file of target regions to limit coverage reporting
                        to (column #4 = gene name), gzipped ok. [default:
                        targets.bed]
  -c, --coverageFile    Bed4 file reporting regions and their coverage class
                        (column #4 = coverage tag), gzipped ok. [default:
                        coverage.bed]
  -k, --keep            Comma separated string of tags identifying low coverage
                        regions. [default: NO_COVERAGE,LOW_COVERAGE]
  -y, --chrY            Name of 'Y' chromosome, if want to filter by sex.
                        [default: NA]
  -r, --regionsFile     Bed4 file of regions with low coverage (column #4 = gene
                        name). [default: LowCover.regions.bed]
  -b, --badGenesFile    Text file listing genes with at least 1 base with low
                        coverage (one per line). [default: LowCover.badgenes.txt]
  -g, --goodGenesFile   Text file listing genes with no low covage (one per
                        line). [default: LowCover.goodgenes.txt]
  -s, --summaryFile     Stats table as a tab-delimited file. [default:
                        LowCover.summary.tsv]
  -S, --summaryFileNoY  Stats table ignoring chrY. Only if --chrY is set.
                        [default: LowCover.summaryNoY.tsv]
  -p, --parallel        How to parallelize. Uses the 'future' R package. One of
                        'default', 'guess', 'multicore', 'multisession', or
                        'sequential' where: 'default' = Uses the plan set before
                        running LowCover, or as defined by the R options used by
                        `future`.  'guess' = Chooses 'multicore' if supported,
                        otherwise choose 'multisession'.  'multicore' = Fork
                        separately to each core/worker. May not not be supported.
                        'multisession' = One R session per core/worker. 
                        'sequential' = Do not run in parallel; uses just one
                        core/worker. [default: default]
  -w, --workers         In parallel mode, set this to limit workers/cpus. By
                        default all available cpus as determined by
                        future::availableCores() are used. Note: this is ignored
                        if --parallel is either unset, 'default', or
                        'sequential'.
```
