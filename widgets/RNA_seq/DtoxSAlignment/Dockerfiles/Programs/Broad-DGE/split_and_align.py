from __future__ import print_function
import os
import os.path
import sys
import gzip
import itertools

## Parameters


bwa_aln_flags = "-l 24"
batch_size = 5000000
r1_length = 16
#r2_length = 33

## Helpers

# Open .fastq or .fastq.gz files for reading
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rb")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rb")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError, "Unknown file: " + filename

def bwa_cmd(fastq_path, reference_prefix, alignment_dir):  ## separated bwa from bsub
    return " ".join(["bwa aln", bwa_aln_flags, reference_prefix, fastq_path, "|", "bwa samse -n 20", \
                    reference_prefix, "-", fastq_path, ">", ".".join([fastq_path, "sam"])])

def write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, name_seq_qual_list, reference_prefix):
    fastq_path = os.path.join(alignment_dir, ".".join([sample_id, subsample_id, str(read_count), "fastq"]))
    with open(fastq_path, "w") as out:
        for name, seq, qual in name_seq_qual_list:
            print("\n".join(["@" + name, seq, "+", qual]), file=out)
    #os.system(bsub_bwa_cmd(fastq_path, reference_prefix, alignment_dir))   
    return bwa_cmd(fastq_path, reference_prefix, alignment_dir)

# Mask sequence by quality score
def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.izip(seq, qual))

## Main

# Parse command line parameters
if len(sys.argv) != 7:
    print("Usage: python " + sys.argv[0] + " <sample id> <subsample id> <read1 .fastq or .fastq.gz> " +\
            "<read2 fastq or .fastq.gz> <alignment dir> <reference prefix>",file=sys.stderr)
    sys.exit()

sample_id, subsample_id, r1_fastq, r2_fastq, alignment_dir, reference_prefix = sys.argv[1:]

# Read through R1 and R2 fastq files in parallel, add R1 barcode to R2 name, and launch batched
# BWA alignments via LSF with output written to alignment_dir
with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:
    read_count = 0
    buf = list()
    cmd_list = list() ### added

    r1_r2 = itertools.izip(r1_file, r2_file)
    for header1, header2 in r1_r2:
        seq1, seq2 = r1_r2.next()
        plus1, plus2 = r1_r2.next()
        qual1, qual2 = r1_r2.next()

        read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
        assert read_name1 == read_name2
        seq2, qual2 = seq2.rstrip(), qual2.rstrip()
        barcode, seq, qual = mask(seq1[0:6], qual1[0:6], min_qual=10) + mask(seq1[6:r1_length], qual1[6:r1_length]), seq2, qual2
        barcoded_name = ":".join([read_name2, barcode])

        buf.append((barcoded_name, seq, qual))
        read_count += 1
        if read_count % batch_size == 0:
            cmd_to_add = write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, buf, reference_prefix)
            cmd_list.append(cmd_to_add)
            buf = list()

    if len(buf) > 0:
    	cmd_to_add = write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, buf, reference_prefix)
    	cmd_list.append(cmd_to_add)
	
	print('2) run individual align commands:')
	for cmd in cmd_list:
		print(cmd)
