#!/bin/bash

output_file="check.log"


nlayers=2
Nfield=4


echo "Code verification for ${nlayers} layers model" > $output_file

for nl in $(seq 1 $nlayers)
do
    index_layer=$((6*(nl-1) + 1))

    sed -n "${index_layer}p" mlswe_ref_FIN.txt >> $output_file

    for ifield in $(seq 2 $Nfield)
    do
        index_field=$((ifield+index_layer))
        # sed -n "${index_field}p" mlswe_ref_FIN.txt | awk '{printf("%-3s  = %10.5f %10.5f\n", $4, ($5-$6), ($6-$6))}' >> $output_file

        field_ref=$(sed -n "${index_field}p" mlswe_ref_FIN.txt | awk '{printf $4}')
        field_max_ref=$(sed -n "${index_field}p" mlswe_ref_FIN.txt | awk '{printf $5}')
        field_min_ref=$(sed -n "${index_field}p" mlswe_ref_FIN.txt | awk '{printf $6}')

        field_max=$(sed -n "${index_field}p" mlswe_FIN.txt | awk '{printf $5}')
        field_min=$(sed -n "${index_field}p" mlswe_FIN.txt | awk '{printf $6}')

        error_field_max=$(awk -v field_max_ref="$field_max_ref" -v field_max="$field_max" 'BEGIN {printf "%10.5f\n", (field_max_ref - field_max)}')
        error_field_min=$(awk -v field_max_ref="$field_min_ref" -v field_max="$field_min" 'BEGIN {printf "%10.5f\n", (field_max_ref - field_max)}')

        printf "%-3s  = %10.5e %10.5e\n" "$field_ref" "$error_field_max" "$error_field_min" >> $output_file
        
    done
done 

cat $output_file
