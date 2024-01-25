import os
import argparse

def list_to_array(lis):
    """
    Converts a python string list into a string bash array format.
    ("<str1>","<str2>","<str3>")
    """
    stringbuilt = ''
    for item in lis[:-1]:
        stringbuilt += '\"'
        stringbuilt += item
        stringbuilt += '\" '
    stringbuilt += '\"'
    stringbuilt += lis[-1]
    stringbuilt += '\"'

    return stringbuilt


def string_output(var, output):
    """
    Transfers a variable and its content (as a string) to an output widget.
    Use a printf statement to pipe into /tmp/output/<output> to properly load content into output.

    var: Python variable to transfer, either string or list
    output: Output parameter name specified in the original widget (in "Outputs" tab)
    """
    
    # Variable is a list, then need to convert it to bash list string
    if type(var) == list:
        var = list_to_array(var)
    
    os.system('printf "{}"  > "/tmp/output/{}"'.format(var, output))



def main():
    parser = argparse.ArgumentParser()

    # Work directory arg
    parser.add_argument('-workdir',
                        type = str,
                        nargs = '+',
                        help = "Work directory where fastq, sam, and count output files are stored.",
                        required = True)

    # sample IDs arg
    parser.add_argument('-ids',
                        type = str,
                        nargs = '+',
                        help = "IDs used in alignment.",
                        required = True)
    
    # Option to gzip downloaded FASTQ files
    parser.add_argument('-gzip',
                        help = "Make fastq files gzipped when downloaded",
                        required = False,
                        action = 'store_true')
    
    args = parser.parse_args()

    # Get flag contents
    work_dir = args.workdir[0]
    ids = args.ids
    gzip = args.gzip
    print("Work directory: {}".format(work_dir))
    print("IDs: {}".format(ids))

    # Make list of expected downloaded fastq files from IDs for 1st and 2nd mates, as well as unpaired fastq
    fastq1 = ["{}/{}_1.fastq".format(work_dir, id) for id in ids]
    fastq2 = ["{}/{}_2.fastq".format(work_dir, id) for id in ids]
    fastq_unpaired = ["{}/{}.fastq".format(work_dir, id) for id in ids]

    # Make fastq files gzipped if gzip=True
    if gzip:
        fastq1 = [fastq + '.gz' for fastq in fastq1]
        fastq2 = [fastq + '.gz' for fastq in fastq2]
        fastq_unpaired = [fastq + '.gz' for fastq in fastq_unpaired]

    # Make list of expected output directories (outputs after Salmon Align step)
    outputs = ["{}/{}_quant".format(work_dir, id) for id in ids]


    # Load variables as outputs (found in Outputs tab)
    string_output(fastq1, 'mate_1')
    string_output(fastq2, 'mate_2')
    string_output(fastq_unpaired, 'unpaired')
    string_output(outputs, 'output')

if __name__ == '__main__':
    main()