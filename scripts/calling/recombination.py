import pysam
import sys
from collections import defaultdict
from utils_v6 import indel_mismatch_column_map, snp_call, build_fastq_single_read, fastq_dicts
bam_in = sys.argv[1]
fa_in =  sys.argv[2]
contig_name = sys.argv[3]
contig_length = sys.argv[4]
len_of_interval_around_mut = int(sys.argv[5])
min_base_qual = int(sys.argv[6])
min_map_qual = int(sys.argv[7])
initial_closest_indel_column = int(sys.argv[8])
minimum_coverage_to_call_snp = int(sys.argv[9])
length_repeat_block = int(sys.argv[10])
context_length_repeat_block = int(sys.argv[11])
unique_repeat_chunks = int(sys.argv[12])
threshold_closest_indel_column = int(sys.argv[13])
threshold_indel_previous_snp = int(sys.argv[14])
threshold_base_to_indel_ratio = float(sys.argv[15])
minimum_ratio_of_haplotype_one_bases = float(sys.argv[16])
minimum_ratio_of_haplotype_two_bases = float(sys.argv[17])
threshold_base_quality_at_gene_conversion = int(sys.argv[18]) 
threshold_base_quality_at_snp = int(sys.argv[19])
maximum_coverage_to_call_snp = int(sys.argv[20])

snp_filter_path = sys.argv[21]
gc_filter_path = sys.argv[22]
snp_distances_path = sys.argv[23]
runs_of_homozygosity_path = sys.argv[24]
single_reads_fastq_path = sys.argv[25]
multiple_reads_fastq_path = sys.argv[26]
control_path = sys.argv[27]

reference = pysam.FastaFile(fa_in)
bam_alignment = pysam.AlignmentFile(bam_in, "rb")

read_count_haplotype_one = defaultdict(int)
read_count_haplotype_two = defaultdict(int)
haplotype_one_reads_in_position = defaultdict(str)
haplotype_two_reads_in_position = defaultdict(str)
inter_snp_distance = defaultdict(int)
haplotype_at_snp = defaultdict(str)
bases_for_read_at_snps = defaultdict(str)
bases_for_hap1_at_snps = defaultdict(str)
bases_for_hap2_at_snps = defaultdict(str)
reads_in_snp = defaultdict(list)
snp_positions = []
indel_type = defaultdict(dict)
mismatch_type = defaultdict(dict)
indel_columns = set()
read_snp_pos = defaultdict(list)
recombinations_in_pos = defaultdict(int)
output_dict = defaultdict(list)


snp_filter = open(snp_filter_path, 'a')
gc_filter = open(gc_filter_path, 'a')
snp_distances = open(snp_distances_path, 'a')
runs_of_homozygosity = open(runs_of_homozygosity_path, 'a')
single_reads_fastq = open(single_reads_fastq_path, 'a')
#multiple_reads_fastq = open(multiple_reads_fastq_path, 'a')
#gc_control = open(control_path, 'a')

indel_mismatch = indel_mismatch_column_map(bam_alignment,
                 contig_name, 
                 min_base_qual, 
                 min_map_qual, 
                 indel_columns, 
                 threshold_base_to_indel_ratio,
                 minimum_coverage_to_call_snp, 
                 maximum_coverage_to_call_snp, 
                 indel_type,
                 mismatch_type,
                 reference,
                 minimum_ratio_of_haplotype_one_bases,
                 minimum_ratio_of_haplotype_two_bases)

indel_type, mismatch_type = indel_mismatch[0], indel_mismatch[1]

snp_call(bam_alignment, 
                              reference, 
                              contig_name, 
                              min_base_qual, 
                              min_map_qual, 
                              minimum_coverage_to_call_snp, 
                              maximum_coverage_to_call_snp, 
                              snp_positions, 
                              reads_in_snp, 
                              read_count_haplotype_one, 
                              read_count_haplotype_two, 
                              indel_columns, 
                              threshold_closest_indel_column, 
                              threshold_base_quality_at_snp, 
                              snp_filter, 
                              minimum_ratio_of_haplotype_one_bases, 
                              minimum_ratio_of_haplotype_two_bases, 
                              threshold_base_to_indel_ratio, 
                              length_repeat_block, 
                              unique_repeat_chunks, 
                              snp_distances, 
                              runs_of_homozygosity, 
                              bases_for_hap1_at_snps, 
                              bases_for_hap2_at_snps, 
                              bases_for_read_at_snps, 
                              haplotype_at_snp,
                              len_of_interval_around_mut,
                              indel_type,
                              mismatch_type,
                              read_snp_pos,
                              gc_filter)



print('contig', 
      'previous_snp', 
      'jump_snp', 
      'next_snp', 
      'read_name',
      'count_haplotype_one',
      'count_haplotype_two', 
      'ref_bases',
      'alt_bases',
      'read_bases',
      'haplotype_string', 
      'number_of_jumps',
      'jump_snp_pos_in_read',
      'snp_pos_in_read',
      'strand',
      'read_length',
      'cigar_string',
      sep='\t')

for read_name, haplotype_string in haplotype_at_snp.items():
    if 'X' in haplotype_string and 'O' in haplotype_string:
        i_index = [index for index, letter in enumerate(haplotype_string) if letter == 'x' or letter == 'X']
        o_index = [index for index, letter in enumerate(haplotype_string) if letter == 'o' or letter == 'O']
        min_x = min(i_index)
        min_o = min(o_index)
        if min_x < min_o:
            jump_snp = min_o
            
            candidate_previous_snps = [index for index in i_index if index < jump_snp]
            if len(candidate_previous_snps) == 0:
                previous_snp_pos_index = jump_snp
            else:
                previous_snp_pos_index = max(candidate_previous_snps)
            
            candidate_next_snps = [index for index in i_index if index > jump_snp]
            if len(candidate_next_snps) == 0:
                candidate_next_snps_index = jump_snp
            else:
                candidate_next_snps_index = min(candidate_next_snps)
        else:
            jump_snp = min_x

            candidate_previous_snps = [index for index in o_index if index < jump_snp]
            if len(candidate_previous_snps) == 0:
                previous_snp_pos_index = jump_snp
            else:
                previous_snp_pos_index = max(candidate_previous_snps)
            
            candidate_next_snps = [index for index in o_index if index > jump_snp]
            if len(candidate_next_snps) == 0:
                candidate_next_snps_index = jump_snp
            else:
                candidate_next_snps_index = min(candidate_next_snps)
        recombinations_in_pos[read_snp_pos[read_name][jump_snp]] += 1
        
        jumps = 0
        for i in range(len(haplotype_string)-1):
            if haplotype_string[i] == 'n':
                continue
            for j in range(i, len(haplotype_string)-1):
                if haplotype_string[i] != haplotype_string[j+1] and haplotype_string[i] == haplotype_string[i].upper() and haplotype_string[j+1] == haplotype_string[j+1].upper():
                    jumps +=1
                    break
                elif haplotype_string[i].upper() == haplotype_string[j+1]:
                   break
        
        for pileupcolumn in bam_alignment.pileup(contig_name, start = read_snp_pos[read_name][jump_snp]-1, 
                                                 end = read_snp_pos[read_name][jump_snp], 
                                                 min_base_quality = min_base_qual, 
                                                 min_mapping_quality = min_map_qual,
                                                 truncate=True):
            jump_snp_pos_in_read = pileupcolumn.get_query_positions()[pileupcolumn.get_query_names().index(read_name)]
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.query_name == read_name:
                    strand = pileupread.alignment.is_reverse
                    read_length = pileupread.alignment.query_alignment_end
                    cigar_string = pileupread.alignment.cigarstring
        
        strand_dict = {True: '-', False: '+' }
        output_dict[read_snp_pos[read_name][jump_snp]] = [contig_name, 
            read_snp_pos[read_name][previous_snp_pos_index], 
            read_snp_pos[read_name][jump_snp],
            read_snp_pos[read_name][candidate_next_snps_index],
            read_name,
            read_count_haplotype_one[read_name],
            read_count_haplotype_two[read_name],
            bases_for_hap2_at_snps[read_name],
            bases_for_hap1_at_snps[read_name],
            bases_for_read_at_snps[read_name],
            haplotype_string,
            jumps,
            jump_snp_pos_in_read,
            ','.join(str(x) for x in read_snp_pos[read_name]),
            strand_dict[strand],
            read_length,
            cigar_string]
                  
previous_pos = 0
for position, count in recombinations_in_pos.items():
    distance_to_prev_event = position - previous_pos
    if count == 1 and distance_to_prev_event > 5000:
        output_list = output_dict[position]
        print('\t'.join([str(ele) for ele in output_list]))
        for pileupcolumn in bam_alignment.pileup(contig_name, start = position-1,  end = position):
            build_fastq_single_read(pileupcolumn, output_list[4])
    previous_pos = position
        

single_read_dict = fastq_dicts()
for read in single_read_dict:
    single_reads_fastq.write('\n'.join(single_read_dict[read])+'\n')

# # for read in multi_read_dict:
# #     multiple_reads_fastq.write('\n'.join(multi_read_dict[read])+'\n')
