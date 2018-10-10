from __future__ import print_function
import os
import sys
import re
import numpy as np

## Functions

def get_mismatches(seq1,seq2):
	mismatches = 0	
	for i in range(0,len(seq1)):
		if seq1[i] != seq2[i] or seq1[i] == "N" or seq2[i] == "N":
			mismatches += 1
	return mismatches
	
def find_well_barcode(seq,barcode_list):
	for barcode in barcode_list:
		mismatches = get_mismatches(seq,barcode)
		if mismatches <= 1:
			return barcode
	return None

## Main

if len(sys.argv) != 8:
    print("Usage: python " + sys.argv[0] + " <sample id> <sym2ref file> <ERCC fasta> <barcodes file> <alignment dir> <dge dir> <loose_barcodes: True/False>",
            file=sys.stderr)
    sys.exit()

sample_id, sym2ref, ercc_fasta, barcodes, aligned_dir, dge_dir, loose_barcodes = sys.argv[1:]

if loose_barcodes == "True":
	loose_barcodes = True
else:
	loose_barcodes = False

# Read gene symbol->RefSeq ID mapping
gene_list = list()
refseq_to_gene = dict()
with open(sym2ref, "rU") as f:
    for line in f:
        items = line.split()
        gene, refseqs = items[0], items[1].split(",")
        gene_list.append(gene)
        refseq_to_gene.update([(rs, gene) for rs in refseqs])
gene_to_idx = {gene: idx for idx, gene in enumerate(gene_list)}

# Read ERCC control IDs
ercc_list = list()
with open(ercc_fasta, "rU") as f:
    for line in f:
        if line[0] == ">":
            items = line.split()
            ercc_list.append(items[0][1:])
ercc_to_idx = {ercc: idx for idx, ercc in enumerate(ercc_list)}

# Read barcode->well mapping
well_list = list()
barcode_list = list()
barcode_to_well = dict()
with open(barcodes, "rU") as f:
    for line in f:
        plate, well, barcode = line.split()
        plate_well = "_".join([plate, well])
        well_list.append(plate_well)
        barcode_list.append(barcode)
        barcode_to_well[barcode] = plate_well

# Create well to array index mapping
well_to_idx = {well: idx for idx, well in enumerate(well_list)}



polyA_tail_pattern = re.compile(r"A{20,}$")
amb_umi_pattern = re.compile(r"[A-Z]*N[A-Z]*")
max_edit_dist = 1
max_best = 20
multimappings_list = [True, False]
multi_file_suffix = ".all"

for multimappings in multimappings_list:
	if not multimappings:
		multi_file_suffix = ".unq"
		
	# Initialize counters
	count_total_reads = 0
	count_assigned_reads = 0
	count_assigned_aligned_reads = 0

	count_spike_total = np.zeros((len(ercc_list), len(well_to_idx)), np.int)
	count_spike_umi = np.zeros((len(ercc_list), len(well_to_idx)), np.int)

	count_total = np.zeros((len(gene_list), len(well_to_idx)), np.int)
	count_umi = np.zeros((len(gene_list), len(well_to_idx)), np.int)

	count_assigned_mito_reads = 0
	count_assigned_mito_umi = 0

	count_assigned_unknown_reads = 0
	count_assigned_unknown_umi = 0

	seen_umi = set()
	unknown = set()

	# Read each partial alignment file and add to counts

	for aligned_filename in [f for f in os.listdir(aligned_dir)
							if f.startswith(sample_id + ".") and f.endswith(".sam")]:
		with open(os.path.join(aligned_dir, aligned_filename), "rU") as f:
			for line in f:
				if line.startswith("@"):
					continue
				count_total_reads += 1
				items = line.split()
				read_id, aligned_id = items[0], items[2]
				barcode = read_id.split(":")[-1]

				well = barcode_to_well.get(barcode[0:6], None)
 				if well == 'Bad':
 					continue
				if well is None:
					if loose_barcodes:
 						well_barcode = find_well_barcode(barcode[0:6], barcode_list)
 						if well_barcode is None:
 							barcode_to_well[barcode[0:6]] = 'Bad'
 							continue
 						well = barcode_to_well.get(well_barcode, None)
 						barcode_to_well[barcode[0:6]] = well
					else:
						continue
				if (re.search(amb_umi_pattern, barcode[6:16])):
					continue
				count_assigned_reads += 1
				if aligned_id == "*":
					continue
				read, edit_dist, best_hits = items[9], int(items[12].split(":")[-1]), int(items[13].split(":")[-1])
				if (edit_dist > max_edit_dist) or (best_hits > max_best) or (re.search(polyA_tail_pattern, read)):
					continue
				best_hits_list = list()
				try:
					items[19]
					best_hits_loc = items[19].split(":")[-1]
					split_loc = best_hits_loc.split(";")
					for loc in split_loc:
						new_loc = loc.split(",")[0]
						if new_loc:
							best_hits_list.append(new_loc)
				except:
					pass
				count_assigned_aligned_reads += 1
				if aligned_id.startswith("ERCC"):
					if (not multimappings) and best_hits_list:
						continue
					umi = "_".join([well, barcode[6:16], aligned_id])
					count_spike_total[ercc_to_idx[aligned_id], well_to_idx[well]] += 1
					if umi not in seen_umi:
						count_spike_umi[ercc_to_idx[aligned_id], well_to_idx[well]] += 1
						seen_umi.add(umi)
				elif aligned_id.startswith("chrM"):
					if (not multimappings) and best_hits_list:
						continue
					umi = "_".join([well, barcode[6:16], aligned_id])
					count_assigned_mito_reads += 1
					if umi not in seen_umi:
						count_assigned_mito_umi += 1
						seen_umi.add(umi)
				elif aligned_id in refseq_to_gene:
					gene = refseq_to_gene[aligned_id]
					multi_gene_hits = False
					if best_hits_list:
						for loc in best_hits_list:
							try: 
								refseq_to_gene[loc]
								if refseq_to_gene[loc] != gene:
									multi_gene_hits = True
							except:
								multi_gene_hits = True
					if (not multimappings) and multi_gene_hits:
						continue                
					umi = "_".join([well, barcode[6:16], gene])
					count_total[gene_to_idx[gene], well_to_idx[well]] += 1
					if umi not in seen_umi:
						count_umi[gene_to_idx[gene], well_to_idx[well]] += 1
						seen_umi.add(umi)
				else:
					if (not multimappings) and best_hits_list:
						continue
					unknown.add(aligned_id)
					umi = "_".join([well, barcode[6:16], aligned_id])
					count_assigned_unknown_reads += 1
					if umi not in seen_umi:
						count_assigned_unknown_umi += 1
						seen_umi.add(umi)

	# Write summary of sample contents
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".log.dat"), "w") as f:
		print("\t".join(["Sample_ID",
						 "Total",
						 "Assigned",
						 "Aligned",
						 "Spike_Total",
						 "Spike_UMI",
						 "Mito_Total",
						 "Mito_UMI",
						 "Refseq_Total",
						 "Refseq_UMI",
						 "Unknown_Total",
						 "Unknown_UMI"]),
				file=f)
		print("\t".join([sample_id] +
						 [str(x) for x in [count_total_reads,
											count_assigned_reads,
											count_assigned_aligned_reads,
											count_spike_total.sum(),
											count_spike_umi.sum(),
											count_assigned_mito_reads,
											count_assigned_mito_umi,
											count_total.sum(),
											count_umi.sum(),
											count_assigned_unknown_reads,
											count_assigned_unknown_umi]]),
				file=f)

	# Write list of RefSeq IDs that were not mapped to genes by UCSC
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".unknown_list"), "w") as f:
		for name in unknown:
			print(name, file=f)

	# Write gene x well total alignment counts
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".refseq.total.dat"), "w") as f:
		print("\t".join([""] + well_list), file=f)
		for gene in gene_list:
			print("\t".join([gene] +
							 [str(x) for x in count_total[gene_to_idx[gene], :]]),
					file=f)

	# Write gene x well UMI counts
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".refseq.umi.dat"), "w") as f:
		print("\t".join([""] + well_list), file=f)
		for gene in gene_list:
			print("\t".join([gene] + \
							 [str(x) for x in count_umi[gene_to_idx[gene], :]]), \
					file=f)

	# Write ERCC x well total alignment counts
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".spike.total.dat"), "w") as f:
		print("\t".join([""] + well_list), file=f)
		for ercc in ercc_list:
			print("\t".join([ercc] + \
							 [str(x) for x in count_spike_total[ercc_to_idx[ercc], :]]), \
					file=f)

	# Write ERCC x well UMI counts
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".spike.umi.dat"), "w") as f:
		print("\t".join([""] + well_list), file=f)
		for ercc in ercc_list:
			print("\t".join([ercc] + \
							 [str(x) for x in count_spike_umi[ercc_to_idx[ercc], :]]), \
					file=f)

	# Write by-well total and UMI counts
	with open(os.path.join(dge_dir, sample_id + multi_file_suffix + ".well_summary.dat"), "w") as f:
		print("\t".join([""] + well_list), file=f)
		print("\t".join(["Refseq_Total"] +[str(x) for x in count_total.sum(axis=0)]), file=f)
		print("\t".join(["Refseq_UMI"] +[str(x) for x in count_umi.sum(axis=0)]), file=f)
		print("\t".join(["Spike_Total"] +[str(x) for x in count_spike_total.sum(axis=0)]), file=f)
		print("\t".join(["Spike_UMI"] +[str(x) for x in count_spike_umi.sum(axis=0)]), file=f)
		

