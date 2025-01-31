# GIndex
A multidollar-BWT based graph index.

## Compile
```bash
git clone --recursive https://github.com/dlcgold/graph_index.git
cd graph_index/gfatools
make -j2
cd ../mfmi
make -j2
cd ..
make -j2
```

## Index graph in GFA file format
``` bash
./gindex index -h
Usage: index [--help] [--version] [--input VAR] [--thread VAR] [--cache VAR] [--output VAR]

Index a given GFA file to perform exact pattern-matching

Optional arguments:
  -h, --help     shows help message and exits 
  -v, --version  prints version information and exits 
  -i, --input    input graph in GFA format 
  -t, --thread   number of threads (default 1) [nargs=0..1] [default: 1]
  -c, --cache    lenght of kmer to preprocess (default 1) [nargs=0..1] [default: 0]
  -o, --output   output prefix 
```

## Query index using reads in FASTA/FASTQ format

``` bash
./gindex query  -h
Usage: query [--help] [--version] [--input VAR] [--query VAR] [--thread VAR] [--full] [--nocache]

Query a given index, queries in FASTA/FASTQ format

Optional arguments:
  -h, --help     shows help message and exits 
  -v, --version  prints version information and exits 
  -i, --input    input prefix for the index 
  -q, --query    query FASTA/FASTQ file 
  -f, --full     compute all nodes path for a match (warning: cache will not be used, very slow) 
  -n, --nocache  ignore cache 
```

## Example

``` bash
./gindex index -i example/test.gfa -o example/test-index -t 5 - c 1
./gindex query -i example/test-index -q example/test.fq
# Matches (as nodes in .gfa):
# @read1    6
# @read2    4
# @read2    4
```

## Experiments
In `/snakemake`you will find a little `snakemake` pipeline to test `gindex`
against [GCSA2](https://github.com/jltsiren/gcsa2.git).
