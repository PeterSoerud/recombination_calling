#!/opt/homebrew/bin/bash

dsb_coordinates="../dsb_map/dsb_t2t_map.bed"
declare -A samples

samples["human"]="hs25 hs35 hs50 ht20 ht45 ht50 ht55 ht60 ht65"
samples["chimp"]="ct15 ct22 ct28 ct32"  
samples["gorilla"]=" gt21 gt22 gt43"
samples["gibbon"]="lt39"
samples["baboon"]="bt15" 
samples["macaque"]="mt5"

for sample in ${samples["human"]}; do
        for type in co gcv; do
            mkdir -p "../human/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "gcv" ] && ( [ "${sample}" = "ht20" ] || [ "${sample}" = "ht45" ] || [ "${sample}" = "ht50" ] || [ "${sample}" = "ht55" ] || [ "${sample}" = "ht60" ] || [ "${sample}" = "ht65" ]); then
                input_coordinates="../classified_reads/${type}/${sample}/${sample}_gcv_germline_jump_snps.bed"
                bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
                printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            if [ "${type}" = "co" ] && ( [ "${sample}" = "ht20" ] || [ "${sample}" = "ht45" ] || [ "${sample}" = "ht50" ] || [ "${sample}" = "ht55" ] || [ "${sample}" = "ht60" ] || [ "${sample}" = "ht65" ]); then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_${interval}.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
            if [ "${type}" = "gcv" ] && ( [ "${sample}" = "hs25" ] || [ "${sample}" = "hs35" ] || [ "${sample}" = "hs50" ] ) ; then
                input_coordinates="../sperm_reads/${type}/${sample}/${sample}_${type}_jump_snps.bed"
                bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
                printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            if [ "${type}" = "co" ] && ( [ "${sample}" = "hs25" ] || [ "${sample}" = "hs35" ] || [ "${sample}" = "hs50" ] ) ; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../sperm_reads/${type}/${sample}/${sample}_co_${interval}.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../human/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done

for sample in ${samples["chimp"]}; do
        for type in co gcv; do
            if [ "${type}" = "gcv" ]; then
            input_coordinates="../classified_reads/${type}/${sample}/${sample}_gcv_germline_jump_snps_human.bed"
            bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../chimp/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
            printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../chimp/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            mkdir -p "../chimp/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_${interval}_human.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../chimp/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../chimp/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done

for sample in ${samples["gorilla"]}; do
        for type in co gcv; do
            if [ "${type}" = "gcv" ]; then
            input_coordinates="../classified_reads/${type}/${sample}/${sample}_gcv_germline_jump_snps_human.bed"
            bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../gorilla/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
            printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../gorilla/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            mkdir -p "../gorilla/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../classified_reads/${type}/${sample}/${sample}_${type}_germline_${interval}_human.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../gorilla/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../gorilla/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done

for sample in ${samples["gibbon"]}; do
        for type in co gcv; do
            if [ "${type}" = "gcv" ]; then
            input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_gcv_jump_snps_human.bed"
            bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../gibbon/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
            printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../gibbon/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            mkdir -p "../gibbon/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_${type}_${interval}_human.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../gibbon/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../gibbon/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done

for sample in ${samples["baboon"]}; do
        for type in co gcv; do
            if [ "${type}" = "gcv" ]; then
            input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_gcv_jump_snps_human.bed"
            bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../baboon/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
            printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../baboon/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            mkdir -p "../baboon/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_${type}_${interval}_human.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../baboon/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../baboon/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done


for sample in ${samples["macaque"]}; do
        for type in co gcv; do
            if [ "${type}" = "gcv" ]; then
            input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_gcv_jump_snps_human.bed"
            bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../macaque/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed
            printf "${sample} jump_snps $(awk '{a[$3]++} END{for(b in a) print b}' ../macaque/${sample}/output_files/${type}/dsb_overlap/${sample}_jump_snps_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
            fi
            mkdir -p "../macaque/${sample}/output_files/${type}/dsb_overlap"
            if [ "${type}" = "co" ]; then
                for interval in jump_snps flanking_snps; do
                    input_coordinates="../unclassified_reads/${type}/${sample}/${sample}_${type}_${interval}_human.bed"
                    bedtools intersect -wa -a $input_coordinates -b $dsb_coordinates > ../macaque/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed
                    printf "${sample} ${interval} $(awk '{a[$3]++} END{for(b in a) print b}' ../macaque/${sample}/output_files/${type}/dsb_overlap/${sample}_${interval}_${type}_dsb_overlap.bed | wc -l) $(wc -l ${input_coordinates})\n" | awk -v OFS='\t' '{print $1, $2, $3, $4}' >> ../dsb_map/overlaps/${type}_overlaps_t2t.bed
                done
            fi
    done
done