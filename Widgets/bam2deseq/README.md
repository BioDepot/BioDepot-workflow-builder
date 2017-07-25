# BAM to counts container


## Installation/Usage:
**Step 1**: Pull the docker image
```bash
$ docker pull biodepot/bam2counts
```

**Step 2**: Run the docker image
```bash
$ docker run -it --rm biodepot/bam2counts
```

**Step 3**: Run the script 
```bash
    Rscript bam2counts.R [GTF File] [Sample Table File] [BAM file Directory] [Workers number]

    Arguments:
    [GTF File] a file specifying the genomic features.
    [Sample Table File] is a table with detailed information for each of samples that links samples to the associated FASTQ and BAM files.
    [BAM file Directory] is a directory contains BAM files.
    [Workers number] is a parameter defined number of threads, default = 1
```
