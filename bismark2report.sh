#!/bin/bash

# Directory containing M-bias text files
input_dir="/hb/groups/kelley_lab/tina/mytilus/04_methylation_extractor/methylation"

# Directory to store HTML reports
#output_dir="/hb/groups/kelley_lab/tina/mytilus/04_methylation_extractor/mbias_plots_and_reports"
#mkdir -p "$output_dir"

# Loop through each M-bias text file in the input directory
for mbias_file in "$input_dir"/*_R1_merged_bismark_bt2_pe.M-bias.txt; do
    # Extract sample name from the file name
    sample_name=$(basename "$mbias_file" _R1_merged_bismark_bt2_pe.M-bias.txt)

    # Extract tissue type from sample name
    tissue_type=$(echo "$sample_name" | awk -F'[_-]' '{print $2}')

    # Determine alignment directory based on tissue type
    if [[ "$tissue_type" == "F" ]]; then
        alignment_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map/F"
    elif [[ "$tissue_type" == "G" ]]; then
        alignment_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map/G"
    else
        echo "Unknown tissue type for sample: $sample_name"
        continue
    fi

    # Find the alignment report file
    alignment_report_file=$(find "$alignment_dir/${sample_name}" -name "*_report.txt")

    # Generate report using bismark2report
    bismark2report --alignment_report "$alignment_report_file" --mbias_report "$mbias_file"

    # Extract tissue type from sample name
    echo "Sample name: $sample_name"
    tissue_type=$(echo "$sample_name" | awk -F'[_-]' '{print $2}')
    echo "Tissue type: $tissue_type"
done



