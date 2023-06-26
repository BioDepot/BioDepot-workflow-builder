#!/bin/bash

# temporary batch command file, this is currently only used for autoDetermineRegions option
temp_batch=/tmp/batch.txt

# parse command line options, this is done in a function to preserve original arguments
parseOptions() {
	while [ $# -ne 0 ]; do
		case $1 in
			-g|--genome)
				shift
				genome_id=$1
				;;
			-b|--batch)
				shift
				batch_file=$1
				;;
			*)
				input_file+=($1)
				;;
		esac
		shift
	done
}

parseOptions $@

# check if batch file was passed along with auto-determine regions option
# this is not supported behavior
if [[ -n "$autoDetermineRegions" && -n "$batch_file" ]]; then
	echo "ERROR: auto-determine region(s) of interest does not work with batch command file input"
	exit 1
# if autoDetermineRegions is set, construct a batch script
elif [ -n "$autoDetermineRegions" ]; then
	# set to default state
	echo new > $temp_batch
	# load a reference genome file if set
	[ -n "$genome_id" ] && echo "genome $genome_id" >> $temp_batch
	for f in "${input_file[@]}"; do
		# verify that file is .maf type
		if ! echo $f | grep -q "maf$\|vcf$"; then
			echo "ERROR: $f is not .maf file type, cannot auto-determine regions of interest"
			exit 1
		fi
		# load .maf or .vcf file
		echo "load $f" >> $temp_batch
	done
	# Load regions of interest, these are parsed directly from the .maf file
	for f in "${input_file[@]}"; do
		if echo $f | grep -q "maf$"; then
			tail -n +3 $f | awk '{print "region",$5,$6,$7}' >> $temp_batch
		elif echo $f | grep -q "vcf$"; then
			awk '!/^#/ {print "region",$1,$2,$2}' $f >> $temp_batch
		fi
	done
	igv.sh -b $temp_batch
else
	igv.sh $@
fi
# exit with return status from igv.sh
exit $?
