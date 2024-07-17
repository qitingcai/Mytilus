#!/bin/bash

# Define the base directory containing the sample folders
base_dir="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/F"
output_file="/hb/groups/kelley_lab/tina/mytilus/03_mapping/02_alignment/final_map0603/final_map_min0.6/bismark_summary.txt"

# Write the header for the output file
echo -e "Sample\tSequence pairs analysed\tUnique best hit alignments\tMapping efficiency\tNo alignments\tDid not map uniquely\tDiscarded pairs\tCT/GA/CT\tGA/CT/CT\tGA/CT/GA\tCT/GA/GA\tComplementary strands rejected\tTotal Cs analysed\tMethylated Cs CpG\tMethylated Cs CHG\tMethylated Cs CHH\tMethylated Cs Unknown\tUnmethylated Cs CpG\tUnmethylated Cs CHG\tUnmethylated Cs CHH\tUnmethylated Cs Unknown\tCpG methylation\tCHG methylation\tCHH methylation\tUnknown context methylation" > $output_file

# Iterate through each sample directory
for sample_dir in $base_dir/*/; do
  sample=$(basename $sample_dir)
  report_file="${sample_dir}/${sample}_R1_merged_bismark_bt2_PE_report.txt"

  # Check if the report file exists
  if [[ -f $report_file ]]; then
    # Initialize variables
    sequence_pairs_analysed=""
    unique_best_hit=""
    mapping_efficiency=""
    no_alignments=""
    did_not_map_uniquely=""
    discarded_pairs=""
    ct_ga_ct=""
    ga_ct_ct=""
    ga_ct_ga=""
    ct_ga_ga=""
    complementary_rejected=""
    total_cs_analysed=""
    methylated_cs_cpg=""
    methylated_cs_chg=""
    methylated_cs_chh=""
    methylated_cs_unknown=""
    unmethylated_cs_cpg=""
    unmethylated_cs_chg=""
    unmethylated_cs_chh=""
    unmethylated_cs_unknown=""
    cpg_methylation=""
    chg_methylation=""
    chh_methylation=""
    unknown_methylation=""
    
    # Extract the required lines and format them
    while IFS=: read -r key value; do
      case "$key" in
        "Sequence pairs analysed in total") sequence_pairs_analysed=$(echo $value | xargs) ;;
        "Number of paired-end alignments with a unique best hit") unique_best_hit=$(echo $value | xargs) ;;
        "Mapping efficiency") mapping_efficiency=$(echo $value | xargs) ;;
        "Sequence pairs with no alignments under any condition") no_alignments=$(echo $value | xargs) ;;
        "Sequence pairs did not map uniquely") did_not_map_uniquely=$(echo $value | xargs) ;;
        "Sequence pairs which were discarded because genomic sequence could not be extracted") discarded_pairs=$(echo $value | xargs) ;;
        "CT/GA/CT") ct_ga_ct=$(echo $value | xargs) ;;
        "GA/CT/CT") ga_ct_ct=$(echo $value | xargs) ;;
        "GA/CT/GA") ga_ct_ga=$(echo $value | xargs) ;;
        "CT/GA/GA") ct_ga_ga=$(echo $value | xargs) ;;
        "Number of alignments to (merely theoretical) complementary strands being rejected in total") complementary_rejected=$(echo $value | xargs) ;;
        "Total number of C's analysed") total_cs_analysed=$(echo $value | xargs) ;;
        "Total methylated C's in CpG context") methylated_cs_cpg=$(echo $value | xargs) ;;
        "Total methylated C's in CHG context") methylated_cs_chg=$(echo $value | xargs) ;;
        "Total methylated C's in CHH context") methylated_cs_chh=$(echo $value | xargs) ;;
        "Total methylated C's in Unknown context") methylated_cs_unknown=$(echo $value | xargs) ;;
        "Total unmethylated C's in CpG context") unmethylated_cs_cpg=$(echo $value | xargs) ;;
        "Total unmethylated C's in CHG context") unmethylated_cs_chg=$(echo $value | xargs) ;;
        "Total unmethylated C's in CHH context") unmethylated_cs_chh=$(echo $value | xargs) ;;
        "Total unmethylated C's in Unknown context") unmethylated_cs_unknown=$(echo $value | xargs) ;;
        "C methylated in CpG context") cpg_methylation=$(echo $value | xargs) ;;
        "C methylated in CHG context") chg_methylation=$(echo $value | xargs) ;;
        "C methylated in CHH context") chh_methylation=$(echo $value | xargs) ;;
        "C methylated in Unknown context (CN or CHN)") unknown_methylation=$(echo $value | xargs) ;;
      esac
    done < $report_file
    
    # Write the extracted values to the output file
    echo -e "${sample}\t${sequence_pairs_analysed}\t${unique_best_hit}\t${mapping_efficiency}\t${no_alignments}\t${did_not_map_uniquely}\t${discarded_pairs}\t${ct_ga_ct}\t${ga_ct_ct}\t${ga_ct_ga}\t${ct_ga_ga}\t${complementary_rejected}\t${total_cs_analysed}\t${methylated_cs_cpg}\t${methylated_cs_chg}\t${methylated_cs_chh}\t${methylated_cs_unknown}\t${unmethylated_cs_cpg}\t${unmethylated_cs_chg}\t${unmethylated_cs_chh}\t${unmethylated_cs_unknown}\t${cpg_methylation}\t${chg_methylation}\t${chh_methylation}\t${unknown_methylation}" >> $output_file
  else
    echo "Report file for sample ${sample} not found."
  fi
done

echo "Summary file created at ${output_file}"
