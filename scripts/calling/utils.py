def indel_mismatch_column_map(bam_in, 
                     contig_name, 
                     min_base_qual, 
                     min_map_qual, 
                     indel_columns, 
                     threshold_base_to_indel_ratio, 
                     minimum_coverage_to_call_snp, 
                     maximum_coverage_to_call_snp,
                     indel_type,
                     mismatch_type,
                     fa_in,
                     minimum_ratio_of_haplotype_one_bases,
                     minimum_ratio_of_haplotype_two_bases):
    
    """
    Runs through all columns and finds all indels in the mapping to the corresponding contig. It also determines 
    whether the current column should be classified as an indel column i.e. more than 10 % of the reads contain 
    an indel in the column. If the current position is an indel column, the position is added to the dictionary 
    called indel_column. This dictionary is used in the SNP calling.
    """
    
    ref = fa_in.fetch(contig_name)
    for pileupcolumn in bam_in.pileup(contig_name,
                                      min_base_quality = min_base_qual, 
                                      min_mapping_quality = min_map_qual, 
                                      truncate=True): # specify in which region we are looking for gene conversion
        
        list_of_bases_in_pileup = [x.upper() for x in pileupcolumn.get_query_sequences(add_indels = True)] # Extracts all bases in the current list_of_bases_in_pileup
        number_of_different_symbols = len(set(list_of_bases_in_pileup)) 
        length_of_nucleotide_set = len(set([i for i in list_of_bases_in_pileup if i.lower() in 'acgt']))
        number_of_bases = len([i for i in list_of_bases_in_pileup if i.lower() in 'acgt']) 
        coverage = pileupcolumn.n 

        column_index = 0
        for character in list_of_bases_in_pileup:
            if character.lower() not in 'acgt': # Checks if we see an indel
                all_read_names = get_info_for_all_reads(list_of_bases_in_pileup, pileupcolumn, ref)[0] # Extracts all read names
                indel_type[all_read_names[column_index]][pileupcolumn.pos+1] = character
            column_index += 1
        
        if number_of_different_symbols > 1 and \
            length_of_nucleotide_set >= 1 and \
            number_of_bases / coverage < threshold_base_to_indel_ratio and \
            coverage > minimum_coverage_to_call_snp and coverage < maximum_coverage_to_call_snp:
            indel_columns.add(pileupcolumn.pos+1)

        if length_of_nucleotide_set > 1:
            base_count_dict = {}
            info_all_reads = get_info_for_all_reads(list_of_bases_in_pileup, pileupcolumn, ref)
            all_read_names, haplotype_one_index, haplotype_two_index = info_all_reads[0], info_all_reads[1], info_all_reads[2]
            if len(haplotype_one_index) / coverage <= minimum_ratio_of_haplotype_one_bases or \
                  len(haplotype_two_index) / coverage <= minimum_ratio_of_haplotype_two_bases:
                for base in set(list_of_bases_in_pileup):
                     if base.lower() in 'acgt':
                        base_count_dict[base] =list_of_bases_in_pileup.count(base)
                mismatch_base = (min(base_count_dict, key=base_count_dict.get)).lower()
                column_index = 0
                for character in list_of_bases_in_pileup:
                    if character.lower() == mismatch_base: # Checks if we see an indel
                        mismatch_type[all_read_names[column_index]][pileupcolumn.pos+1] = character
                    column_index += 1
    return indel_type, mismatch_type      
                    
def build_fastq_single_read(pileupcolumn, 
                            recombination_read_name):
    
    """
    When a potential recombination event is called, we extract the read name, read sequence and read qualities
    to later build a fastq file from that dictionary of reads. We distinguish between events that are seen on
    a single read and events that are seen on multiple reads.
    """

    for pileupread in pileupcolumn.pileups:
        if pileupread.alignment.query_name == recombination_read_name:
            quality = ''.join([chr(x+33) for x in pileupread.alignment.get_forward_qualities()])
            single_read_dict[recombination_read_name] = ['@{read_name}'.format(read_name = pileupread.alignment.query_name), 
                                                         pileupread.alignment.get_forward_sequence(), '+', quality]

def build_fastq_multiple_reads(pileupcolumn, 
                               recombination_read_names):
    
    """
    When a potential recombination event is called, we extract the read name, read sequence and read qualities to later build a fastq file from that dictionary of reads. 
    We distinguish between events that are seen on a single read and events that are seen on multiple reads.
    """
    
    global multi_read_dict
    for pileupread in pileupcolumn.pileups:
        for read_name in recombination_read_names:
            if pileupread.alignment.query_name == read_name and read_name not in multi_read_dict:
                quality = ''.join([chr(x+33) for x in pileupread.alignment.get_forward_qualities()])
                multi_read_dict[read_name] = ['@{read_name}'.format(read_name = pileupread.alignment.query_name), pileupread.alignment.get_forward_sequence(), '+', quality]

def reverse_complement(seq):
    
    """
    Takes a nucleotide sequence (string) and translate it to its reverse complement
    """
    
    reverse_complement = ''
    watson_crick_pairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverse = seq[::-1]
    for base in reverse:
            reverse_complement += watson_crick_pairs[base]
    return reverse_complement

def get_info_for_all_reads(list_of_bases_in_pileup, 
                           pileupcolumn, 
                           ref):
    
    """
    Extracts the read names in the current pileupcolumn. It also separates the reads into the two haplotypes:
    If a base in the pileupcolumn does not agree with the base in the de novo assembly, the read is labled haplotype one
    If a base in the pileupcolumn does agree with the base in the de novo assembly, the read is labled haplotype two
    """
    
    read_names = pileupcolumn.get_query_names()
    haplotype_one_index = [i for i, base in enumerate(list_of_bases_in_pileup) if base != ref[pileupcolumn.pos] and base.lower() in 'acgt'] 
    haplotype_two_index = [i for i, base in enumerate(list_of_bases_in_pileup) if base == ref[pileupcolumn.pos]] 
    base_qualities = pileupcolumn.get_query_qualities()
    return read_names, haplotype_one_index, haplotype_two_index, base_qualities

def haplotype_counter(pileupcolumn,
                      contig_name,
                      threshold_base_quality_at_snp, 
                      list_of_bases_in_pileup, 
                      read_count_haplotype_one, 
                      read_count_haplotype_two, 
                      haplotype_one_reads, 
                      haplotype_two_reads, 
                      ref, 
                      bases_for_hap1_at_snps, 
                      bases_for_hap2_at_snps, 
                      bases_for_read_at_snps, 
                      haplotype_at_snp,
                      len_of_interval_around_mut,
                      indel_type,
                      mismatch_type,
                      read_snp_pos,
                      output_file):
    """
    This function is called when we are at a SNP position, and it assigns the haplotype to each
    read. If a read has a lower quality than the threshold for base qualities at SNPs, we do
    not assign haplotype to that read. For every read we create a haplotype string, to keep
    of how many times it change haplotypes. This information is later applied to distinguish
    between gene conversion, crossovers, boundary cases, and complex events. Furthermore,
    for each read we create a string that keeps track of which bases are seen in the SNPs
    that a read covers. This information is particularly interesting regarding GC-biased
    gene conversion. If a read contains an indel in the current SNP, we add an 'i' as the
    base
    """
    read_names = pileupcolumn.get_query_names()
    quality = pileupcolumn.get_query_qualities()
    hap1_base = [base for base in set(list_of_bases_in_pileup) if base!= ref[pileupcolumn.pos] and base.lower() in 'acgt' and len(base) == 1][0]
    hap2_base = ref[pileupcolumn.pos]
    
    for pileupread in pileupcolumn.pileups: # Iterates through all reads at a given list_of_bases_in_pileup
        read_name = pileupread.alignment.query_name 
        base_quality_gene_conversion_read = quality[read_names.index(read_name)]
        indel_filter_status = indel_filter(pileupcolumn, contig_name, len_of_interval_around_mut, read_name, indel_type, output_file)
        mismatch_filter_status = mismatch_filter(pileupcolumn, contig_name, len_of_interval_around_mut, read_name, mismatch_type, output_file)
        big_indel_status = big_indel_filter(pileupcolumn, contig_name, indel_type, read_name, output_file)

        read_snp_pos[read_name].append(pileupcolumn.pos+1)
 
        if read_name in haplotype_one_reads: # Checks if the read name is in the list of haplotype one reads
            bases_for_hap1_at_snps[read_name] += hap1_base
            bases_for_hap2_at_snps[read_name] += hap2_base
            bases_for_read_at_snps[read_name] += hap1_base
            if base_quality_gene_conversion_read < threshold_base_quality_at_snp or indel_filter_status == 'indel_filter_failed' or mismatch_filter_status == 'mismatch_filter_failed' or big_indel_status == 'big_indel_filter_failed':
                haplotype_at_snp[read_name] += 'o'
            else:
                haplotype_at_snp[read_name] += 'O'
            read_count_haplotype_one[read_name] +=1 # If this is the case, we add one
            continue
        
        elif read_name in haplotype_two_reads: # The same goes for haplotype two
            bases_for_hap2_at_snps[read_name] += hap2_base
            bases_for_hap1_at_snps[read_name] += hap1_base 
            bases_for_read_at_snps[read_name] += hap2_base
            if base_quality_gene_conversion_read < threshold_base_quality_at_snp or indel_filter_status == 'indel_filter_failed' or mismatch_filter_status == 'mismatch_filter_failed' or big_indel_status == 'big_indel_filter_failed':
                haplotype_at_snp[read_name] += 'x'
            else:
                haplotype_at_snp[read_name] += 'X'
            read_count_haplotype_two[read_name] +=1
            continue
        
        elif read_name not in haplotype_one_reads and read_name not in haplotype_two_reads:
            bases_for_hap1_at_snps[read_name] += hap1_base
            bases_for_hap2_at_snps[read_name] += hap2_base
            bases_for_read_at_snps[read_name] += 'n'
            haplotype_at_snp[read_name] +='n'
            continue

def snp_call(bam_in, 
             fa_in,
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
             output_file):
    
    """
    This function calls all SNVs in a contig. We define a position as a SNP if:
    1) The number of different nucleotides is equal to two
    2) The coverage is within a lower and upper threshold
    3) The surrounding region is not too repetitive based on a given threshold
    4) The distance to the closest indel column is less than a given threshold
    5) The SNV cannot occur at the boundary of two homopolymers
    """
    ref = fa_in.fetch(contig_name)
    global single_read_dict 
    single_read_dict = {}
    global multi_read_dict 
    multi_read_dict = {}

    for pileupcolumn in bam_in.pileup(contig_name,
                                      min_base_quality = min_base_qual, 
                                      min_mapping_quality = min_map_qual, 
                                      truncate=True): # specify in which region we are looking for gene conversion
        
        list_of_bases_in_pileup = [x.upper() for x in pileupcolumn.get_query_sequences(add_indels = True)] # Extracts all bases in the current list_of_bases_in_pileup
        coverage = pileupcolumn.n
        number_of_bases = len([i for i in list_of_bases_in_pileup if i.lower() in 'acgt']) 
        number_of_different_symbols = len(set(list_of_bases_in_pileup))
        number_of_different_nucleotides = len(set([i for i in list_of_bases_in_pileup if i.lower() in 'acgt']))
        if number_of_different_nucleotides > 2:
            snp_filter.write('no_snp_not_two_alleles\t{position}\tnumber_of_symbols:{symbols}\n'.format(
                position = pileupcolumn.pos+1, 
                symbols = number_of_different_nucleotides))
            continue

        if number_of_different_symbols > 1 and \
            number_of_different_nucleotides == 2 and \
                  number_of_bases / coverage >= threshold_base_to_indel_ratio and \
                      coverage > minimum_coverage_to_call_snp and \
                        coverage < maximum_coverage_to_call_snp:
          
            info_all_reads = get_info_for_all_reads(list_of_bases_in_pileup, pileupcolumn, ref)
            read_names, haplotype_one_index, haplotype_two_index = info_all_reads[0], info_all_reads[1], info_all_reads[2]
           
            
            if len(haplotype_one_index) / coverage > minimum_ratio_of_haplotype_one_bases and \
                len(haplotype_one_index) > 2 and \
                len(haplotype_two_index) > 2 and \
                len(haplotype_two_index) / coverage > minimum_ratio_of_haplotype_two_bases:
                repeat_block = ref[pileupcolumn.pos-length_repeat_block//2:pileupcolumn.pos+length_repeat_block//2]
                if len(set([repeat_block[i:i+4] for i in range(len(repeat_block))])) < unique_repeat_chunks:
                    snp_filter.write('no_snp_repeat\t{position}\tblock_length:{block_length}\n'.format(
                        position = pileupcolumn.pos+1, 
                        block_length = len(set([repeat_block[i:i+4] for i in range(len(repeat_block))]))))
                    continue
                
                info_homopolymer_filter = homopolymer_filter(ref, pileupcolumn, list_of_bases_in_pileup)
                homopolymer_filter_status, homopolymer_length = info_homopolymer_filter[0], info_homopolymer_filter[1]
                if homopolymer_filter_status == 'homopolymer_filter_failed':
                    snp_filter.write('no_snp_homopolymer\t{position}\thomopolymer_length:{length}\n'.format(
                        position = pileupcolumn.pos+1, 
                        length = homopolymer_length))
                    continue
                
                if len(indel_columns) == 0:
                    closest_indel_column = threshold_closest_indel_column + 1
                else:
                    closest_indel_column = min([abs(pileupcolumn.pos+1-indel_col) for indel_col in indel_columns])
            
                if len(indel_columns) > 0 and closest_indel_column < threshold_closest_indel_column:
                    snp_filter.write('no_snp_indel_column\t{position}\tdist_to_indel_col:{distance}\n'.format(
                        position = pileupcolumn.pos+1,
                        distance = abs(closest_indel_column)))
                    continue
                
                else:
                    haplotype_one_reads = [read_names[i] for i in haplotype_one_index] # Extracts the read names having haplotype one for the current position
                    haplotype_two_reads = [read_names[i] for i in haplotype_two_index] # Extracts the read names having haplotype two for the current position
                
                    if len(snp_positions) >= 1:
                        hap1_base = [base for base in set(list_of_bases_in_pileup) if base!= ref[pileupcolumn.pos] and base.lower() in 'acgt' and len(base) == 1][0]
                        hap2_base = ref[pileupcolumn.pos]
                        distance_between_snps = pileupcolumn.pos+1 - snp_positions[-1]
                        snp_distances.write('{contig}\t{position}\t{ref_base}\t{alt_base}\t{distance}\n'.format(
                            contig = contig_name, 
                            position = pileupcolumn.pos+1,
                            ref_base = hap1_base,
                            alt_base = hap2_base,
                            distance = distance_between_snps))
                        if distance_between_snps > 10**5:
                            runs_of_homozygosity.write('{contig}\t{previous_snp}\t{current_snp}\t{distance}\n'.format(
                                contig = contig_name, 
                                previous_snp = snp_positions[-1], 
                                current_snp = pileupcolumn.pos+1, 
                                distance = distance_between_snps))
                    snp_positions.append(pileupcolumn.pos+1)
                    reads_in_snp[pileupcolumn.pos+1] = pileupcolumn.get_query_names()
                    haplotype_counter(pileupcolumn,
                                      contig_name,  
                                      threshold_base_quality_at_snp, 
                                      list_of_bases_in_pileup, 
                                      read_count_haplotype_one, 
                                      read_count_haplotype_two, 
                                      haplotype_one_reads, 
                                      haplotype_two_reads, 
                                      ref, 
                                      bases_for_hap1_at_snps, 
                                      bases_for_hap2_at_snps, 
                                      bases_for_read_at_snps, 
                                      haplotype_at_snp,
                                      len_of_interval_around_mut,
                                      indel_type,
                                      mismatch_type,
                                      read_snp_pos,
                                      output_file)
    return snp_positions

def homopolymer_filter(ref,
                       pileupcolumn,
                       list_of_bases_in_pileup):
    """
    This function checks if the current position is at the boundary between two homopolymers. 
    """
    hap1_base = [base for base in set(list_of_bases_in_pileup) if base!= ref[pileupcolumn.pos] and base.lower() in 'acgt' and len(base) == 1][0]
    hap2_base = ref[pileupcolumn.pos]
    homopolymers_around_snp = [ref[pileupcolumn.pos-7:pileupcolumn.pos], ref[pileupcolumn.pos+1:pileupcolumn.pos+8]] # Lists of seven bases before the SNP and seven bases after the SNP
    for index, homopolymer in enumerate(homopolymers_around_snp):
        counter = 1
        if index == 0:
            for i in range(1, len(homopolymer)):
                if homopolymer[-i] == homopolymer[-i-1]:
                    counter +=1
                else:
                    break
            homopolymer_length_before_snp = counter
        else:
            for i in range(len(homopolymer)-1):
                if homopolymer[i] == homopolymer[i+1]:
                    counter +=1
                else:
                    break
            homopolymer_length_after_snp = counter
    if homopolymer_length_before_snp + homopolymer_length_after_snp > 6 and \
        set([hap1_base, hap2_base]) == set([homopolymers_around_snp[0][-1], homopolymers_around_snp[1][0]]):
        return 'homopolymer_filter_failed', homopolymer_length_before_snp + homopolymer_length_after_snp
    else:
        return 'homopolymer_filter_passed', homopolymer_length_before_snp + homopolymer_length_after_snp

def indel_filter(pileupcolumn,
                 contig_name,
                 len_of_interval_around_mut, 
                 read_name, 
                 indel_type,
                 output_file):
    
    for i in list(range(-len_of_interval_around_mut, 0))+list(range(1, len_of_interval_around_mut+1)):
        try:
            if read_name in indel_type and \
                    pileupcolumn.pos+1+i in indel_type[read_name].keys():
                output_file.write('{contig}\t{position}\t{read_name}\tclose_indel:{distance}\n'.format(
                    contig = contig_name, 
                    position = pileupcolumn.pos+1, 
                    read_name = read_name, 
                    distance = abs((pileupcolumn.pos+1)-(pileupcolumn.pos+1+i))))
                return 'indel_filter_failed'
        except IndexError:
            break
    return 'indel_filter_passed'

def mismatch_filter(pileupcolumn,
                    contig_name, 
                    len_of_interval_around_mut, 
                    read_name, 
                    mismatch_type,
                    output_file):
    for i in list(range(-len_of_interval_around_mut, 0)) + list(range(1, len_of_interval_around_mut + 1)): # This for loop checks if there is anything close to the SNP in the read that we are looking at. 
        try:
            if pileupcolumn.pos+1+i in mismatch_type[read_name].keys(): 
             # If there is a mutation/error close to the SNP we change the 'haplotype_one_status' to 'bad'. However, we accept the nearby mismatch if it is a SNP
                output_file.write('{contig}\t{position}\t{read_name}\tclose_mismatch:{distance}\n'.format(
                    contig = contig_name, 
                    position = pileupcolumn.pos+1, 
                    read_name = read_name,
                    distance = abs((pileupcolumn.pos+1)-(pileupcolumn.pos+1+i))))
                return 'mismatch_filter_failed'
        except IndexError:
                    break
    return 'mismatch_filter_passed'

def big_indel_filter(pileupcolumn, 
                     contig_name,
                     indel_type, 
                    read_name,
                    output_file):
    
    if read_name in indel_type:
        candidate_read_indels = list(indel_type[read_name].keys())
        candidate_read_indels_types = list(indel_type[read_name].values())
        distance_to_closest_indel = [abs(pileupcolumn.pos + 1-candidate_read_indel) for candidate_read_indel in candidate_read_indels]
        long_indels = [type for type in candidate_read_indels_types if len(type) - 4 > 10]

        for long_indel in long_indels:
            if distance_to_closest_indel[candidate_read_indels_types.index(long_indel)] < 250:
                output_file.write('{contig}\t{position}\t{read_name}\tbig_indel:{indel_size}\n'.format(
                contig = contig_name, 
                position = pileupcolumn.pos + 1, 
                read_name = read_name, 
                indel_size = len(long_indel) - 4))
                return 'big_indel_filter_failed'
                
def fastq_dicts():
    """
    Outputs the final dictionaries of the single and multi read recombinations. 
    """
    return single_read_dict