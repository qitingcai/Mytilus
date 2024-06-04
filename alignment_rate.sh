#!/bin/bash

#extract alignment % from error files

# Set the working directory to where the .err files are located
cd /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map/G

# Create an output file for the alignment summary
output_file="alignment_summary.txt"

# Write the header for the summary file
echo -e "Sample\tOverall Alignment Rate" > "$output_file"

# Loop through all .err files in the directory
for err_file in *.err; do
  # Extract the sample name from the content of the .err file
  sample_name=$(grep ">>> Writing bisulfite mapping results to" "$err_file" | awk -F 'to ' '{print $2}' | awk -F '_R1_merged_bismark_bt2_pe.bam' '{print $1}')

  # Extract the overall alignment rate from the .err file
  alignment_rate=$(grep "overall alignment rate" "$err_file" | head -n 1 | awk '{print $1}')

  # If the alignment rate is not found, default to 0%
  if [ -z "$alignment_rate" ]; then
    alignment_rate="0%"
  fi

  # Print debugging information
  echo "Processing file: $err_file"
  echo "Extracted sample name: $sample_name"
  echo "Extracted alignment rate: $alignment_rate"

  # Write the sample name and alignment rate to the summary file
  echo -e "$sample_name\t$alignment_rate" >> "$output_file"
done


