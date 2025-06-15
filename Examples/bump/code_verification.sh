#!/bin/bash

#!/bin/bash

# Configuration
output_file="Error.txt"
nlayers=2
Nfield=4

# Initialize output file
echo "Code verification for ${nlayers} layers model" > "$output_file"

# Loop through each layer
for nl in $(seq 1 "$nlayers"); do
    index_layer=$((6 * (nl - 1) + 1))

    # Write layer header to output file
    sed -n "${index_layer}p" ref_mlswe_FIN.txt >> "$output_file"

    # Loop through each field
    for ifield in $(seq 2 "$Nfield"); do
        index_field=$((ifield + index_layer))

        # Extract reference data
        field_ref=$(sed -n "${index_field}p" _ref_bump.CI | awk '{print $4}')
        field_max_ref=$(sed -n "${index_field}p" _ref_bump.CI | awk '{print $5}')
        field_min_ref=$(sed -n "${index_field}p" _ref_bump.CI | awk '{print $6}')

        # Extract computed data
        field_max=$(sed -n "${index_field}p" mlswe_FIN.txt | awk '{print $5}')
        field_min=$(sed -n "${index_field}p" mlswe_FIN.txt | awk '{print $6}')

        # Calculate errors

	error_field_max=$(awk -v ref="$field_max_ref" -v val="$field_max" 'BEGIN {printf "%10.5f", ((ref - val) < 0.0) ? -(ref - val) : (ref - val)}')
        error_field_min=$(awk -v ref="$field_min_ref" -v val="$field_min" 'BEGIN {printf "%10.5f", ((ref - val) < 0.0) ? -(ref - val) : (ref - val)}')

	abs_ref_max=$(awk -v ref="$field_max_ref" 'BEGIN {printf "%10.5f", (ref < 0) ? -ref : ref}')
	abs_ref_min=$(awk -v ref="$field_min_ref" 'BEGIN {printf "%10.5f", (ref < 0) ? -ref : ref}')

	res_error_max=$(awk -v n1="$error_field_max" -v n2="$abs_ref_min" 'BEGIN {print n1 / n2}')
	res_error_min=$(awk -v n1="$error_field_min" -v n2="$abs_ref_min" 'BEGIN {print n1 / n2}')

        # Write results to output file
        printf "%-3s  = %10.5e %10.5e\n" "$field_ref" "$error_field_max" "$error_field_min" >> "$output_file"
        #printf "%-3s  = %10.5e %10.5e\n" "$field_ref" "$res_error_max" "$res_error_min" >> "$output_file"
    done
done

cat $output_file
