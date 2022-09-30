# LowCover
Identifies low coverage in selected genomic regions (e.g. genes) based on mosdepth-like coverage info (e.g. genome wide coverage annotation).

## Running locally

### Install package
Can be installed from GitHub within R if you have `BiocManager` and `remotes` installed

```R
> # install.packages("BiocManager")
> # install.packages("remotes")
> BiodManager::install("Jefferys/LowCover")
```

The main function `LowCoverApp()` is designed to be called from a script, as it automatically loads the R script command line parameters and parses them, but it can be run in an R session by passing command line options and values as a vector, e.g.:

```
> LowCoverApp( c( "--targetsFile", "targets.bed",
                  "--coverageFile", "mozdepth.coverage.gz" ))
```

### Creating a script

Designed to be called from a shim script *you create*, e.g.

**LowCover**

```sh
!/usr/bin/env Rscript

LowCover::LowCoverApp()
```

### Running the script

Running this is done from the shell command line like a normal program, generally passing `--targetsFile` and `--coverageFile`.

```sh
$ LowCover --targetFile "targets.bed" --coverageFile "mozdepth.coverage.gz"
```

You can copy this script somewhere on your path to make it globally available.

### Options

Help is available, e.g.

```sh
$ LowCover --help
```

Which currently returns

```
usage: LowCover [--] [--help] [--version] [--force] [--targetsFile
       TARGETSFILE] [--coverageFile COVERAGEFILE] [--keep KEEP] [--chrY
       CHRY] [--regionsFile REGIONSFILE] [--badGenesFile BADGENESFILE]
       [--goodGenesFile GOODGENESFILE] [--summaryFile SUMMARYFILE]
       [--summaryFileNoY SUMMARYFILENOY]

Low coverage reporting

Reports low coverage genomic regions based on mosdepth.

flags:
  -h, --help            show this help message and exit
  -V, --version         Print version info to stderr and exit.
  -f, --force           Create any non-existing directories and
                        overwrite any existing files.

optional arguments:
  -t, --targetsFile     [Required] Bed4 file of target regions to limit
                        coverage reporting to (colum #4 = gene name),
                        gzipped ok. [default: targets.bed]
  -c, --coverageFile    [Required] Bed4 file reporting regions and
                        their coverage class (column #4 = coverage
                        tag), gzipped ok. [default: coverage.bed]
  -k, --keep            Comma separated string of tags identifying low
                        coverage regions. [default:
                        NO_COVERAGE,LOW_COVERAGE]
  -y, --chrY            Name of 'Y' chromosome, if want to filter by
                        sex. [default: NA]
  -r, --regionsFile     Bed4 file of regions with low coverage (column
                        #4 = gene name) [default: LowCover.regions.bed]
  -b, --badGenesFile    Text file listing genes with at least 1 base
                        with low coverage (one per line) [default:
                        LowCover.badgenes.txt]
  -g, --goodGenesFile   Text file listing genes with no low covage (one
                        per line) [default: LowCover.goodgenes.txt]
  -s, --summaryFile     Stats table [default: LowCover.summary.tsv]
  -S, --summaryFileNoY  Stats table ignoring chrY. Ignored if --chrY
                        not set [default: LowCover.summaryNoY.tsv]
```

## Running containerized

### Docker

Change to the Docker directory and use `docker build ...` in the normal way to build the container. You can copy the Docker directory anywhere to build, it needs only the contents of this directory.

Run using `docker run ... `, mounting the directory with input data and the working directory where output files will be written. The normal `LowCover <options>` command can be run as described above, or a shell can be run and then LowCover run within that.
