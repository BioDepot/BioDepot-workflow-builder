# samtools container


## Installation/Usage:
**Step 1**: Pull the docker image
```bash
$ docker pull biodepot/samtools
```

**Step 2**: Run the docker image
```bash
$ docker run -it --rm biodepot/samtools
```

**Step 3**: Run the script 
```bash
$ samtools view -bS SRR1039508.Aligned.out.sam  -o SRR1039508.bam
```