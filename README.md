# iLSLS
A Novel Scaffolding Algorithm Based on Contig Error Correction and Path Extension

## pre-install

perl, python(2.7 and up), bowtie2

if you're not sure about whether your computer satisfy the environment, please run "./src/compile.sh" to test, which would remind you the lacking configuration.

## install

get into ./src/ fold, and input "make"

## usage

please follw the commandline:

```
./bin/iLSLS your_contig_file.fasta pair_end1.fastq pair_end2.fastq
```

The results file will be generated in your working fold, named by "Scaffold_output.fasta"