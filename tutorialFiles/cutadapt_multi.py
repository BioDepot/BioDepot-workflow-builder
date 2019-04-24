import os
from glob import glob
from optparse import OptionParser


def runCutAdapt(file1, file2, flags):
    print(
        "cutadapt {} -o tmp/{} -p tmp/{} {} {}".format(
            flags, file1, file2, file1, file2
        )
    )
    os.system(
        "cutadapt {} -o tmp/{} -p tmp/{} {} {}".format(
            flags, file1, file2, file1, file2
        )
    )


# get options - we are interested in -d for directory -q for quality -m for minimum length
parser = OptionParser()
parser.add_option("-d")
parser.add_option("-q")
parser.add_option("-m")
(options, args) = parser.parse_args()

flags = "-q {} -m {}".format(options.q, options.m)
# change directory to output directory
os.chdir(options.d)

# we use the fact that for our naming convention the paired end files will be nicely grouped in pairs
files = sorted(glob("SRR*.gz"))

# make a the temporary directory
if not os.path.exists("tmp"):
    os.makedirs("tmp")

# run cutadapt on pairs
for i in range(0, len(files) / 2):
    runCutAdapt(files[2 * i], files[2 * i + 1], flags)

# copy the filtered files and remove the temp directory
os.system("cp -r tmp/* . && rm -r tmp")
