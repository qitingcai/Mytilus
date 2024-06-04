#!/bin/bash

# Set the working directory to where the .err and .out files are located
cd /hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603

# Directory to copy and rename the files
output_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map"
mkdir -p "$output_dir"

# Function to rename a file based on the sample name
rename_file() {
  file=$1
  sample_name=$(grep ">>> Writing bisulfite mapping results to" "$file" | awk -F 'to ' '{print $2}' | awk -F '_R1_merged_bismark_bt2_pe.bam' '{print $1}')

  if [ -n "$sample_name" ]; then
    extension="${file##*.}"
    new_name="${sample_name}.${extension}"
    cp "$file" "$output_dir/$new_name"
    echo "Renamed $file to $new_name"
  else
    echo "Sample name not found in $file"
  fi
}

# Loop through all .err and .out files in the directory
for file in *.err *.out; do
  if [ -f "$file" ]; then
    rename_file "$file"
  fi
done

#And then move F and G samples to the corresponding directory
#mv *-G_*.err G/
#mv *-F_*.err F/