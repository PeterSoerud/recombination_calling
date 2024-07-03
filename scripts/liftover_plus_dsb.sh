#!/opt/homebrew/bin/bash

pantro6_to_hg38="../chimpanzee/liftover/panTro6ToHg38.over.chain"
gorgor4_to_hg38="../western_gorilla/liftover/gorGor4ToHg38.over.chain"
hg38_to_t2t="../human/liftover/hg38_to_chm13v1.chain"
dsb_coordinates="../dsb_map/dsb_t2t_map.bed"
declare -A samples

samples["human"]="hs25 hs35 hs50 ht20 ht55"
samples["chimpanzee"]="ct15 ct32" 
samples["western_gorilla"]="gt22 gt43"

for sample in ${samples["human"]}; do
        for type in co nco; do
            mkdir -p "../human/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "nco" ] && ( [ "${sample}" = "ht20" ] || [ "${sample}" = "ht55" ] ); then
                input_coordinates="../classified_reads/${type}/${sample}/${sample}_nco_germline_jump_snps.bed"
                bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
                printf "${sample} ${interval} $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            if [ "${type}" = "co" ] && ( [ "${sample}" = "ht20" ] || [ "${sample}" = "ht55" ] ); then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_${interval}.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
            if [ "${type}" = "nco" ] && ( [ "${sample}" = "hs25" ] || [ "${sample}" = "hs35" ] || [ "${sample}" = "hs50" ] ) ; then
                input_coordinates="../human/${sample}/output_files/${type}/${sample}_jump_snps_${type}.bed"
                bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
                printf "${sample} ${interval} $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            if [ "${type}" = "co" ] && ( [ "${sample}" = "hs25" ] || [ "${sample}" = "hs35" ] || [ "${sample}" = "hs50" ] ) ; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../human/${sample}/output_files/${type}/${sample}_${interval}_co.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done


for species in chimpanzee western_gorilla; do
    for sample in ${samples[$species]}; do
        for type in co nco; do
            mkdir -p "../${species}/${sample}/output_files/${type}/lifted_coordinates"
            mkdir -p "../${species}/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "nco" ]; then
                output_lifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_jump_snps_${type}_hg38.bed"
                output_unlifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_jump_snps_${type}_hg38_unlifted.bed"
                input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_jump_snps.bed"
                if [ "${species}" = "chimpanzee" ]; then
                    liftOver $input_coordinates $pantro6_to_hg38 $output_lifted $output_unlifted
                fi
                if [ "${species}" = "western_gorilla" ]; then
                    liftOver $input_coordinates $gorgor4_to_hg38 $output_lifted $output_unlifted
                fi
                new_input_coordinates=$output_lifted
                new_output_lifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_jump_snps_${type}_t2t.bed"
                new_output_unlifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_jump_snps_${type}_t2t_unlifted.bed"
                liftOver $new_input_coordinates $hg38_to_t2t $new_output_lifted $new_output_unlifted
                bedtools intersect -wa -a $new_output_lifted -b $dsb_coordinates > ../${species}/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
                printf "${sample} jump_snps $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../${species}/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${new_output_lifted})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed   
            fi
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    output_lifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_${interval}_${type}_hg38.bed"
                    output_unlifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_${interval}_${type}_hg38_unlifted.bed"
                    input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_${interval}.bed"
                    if [ "${species}" = "chimpanzee" ]; then
                        liftOver $input_coordinates $pantro6_to_hg38 $output_lifted $output_unlifted
                    fi
                    if [ "${species}" = "western_gorilla" ]; then
                        liftOver $input_coordinates $gorgor4_to_hg38 $output_lifted $output_unlifted
                    fi
                    new_input_coordinates=$output_lifted
                    new_output_lifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_${interval}_${type}_t2t.bed"
                    new_output_unlifted="../${species}/${sample}/output_files/${type}/lifted_coordinates/${sample}_${interval}_${type}_t2t_unlifted.bed"
                    liftOver $new_input_coordinates $hg38_to_t2t $new_output_lifted $new_output_unlifted
                    bedtools intersect -wa -a $new_output_lifted -b $dsb_coordinates > ../${species}/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk 'NR>1{a[$3]++} END{for(b in a) print b}' ../${species}/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${new_output_lifted})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
        done
    done
done


