# Commands to filter and format INDELs in the discovery cohort (Peifer+Berlin)

### Get the INDELs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -P "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_INDELS_mutect2_7.0_nofiltering_newstrandinfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/mutect2/bwa.mutect2."$i"-T1-DNA1-WGS1.full.vcf.gz "$i; done | sh

### Format and filter PASS INDELs from Peifer cohort 
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -P "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_indels_format_4.0.py INDELS_mutect_indels_patient_"$i"_nofilters_newinfostrands.txt "$i; done | sh

### Get the INDELs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -P "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_INDELS_mutect2_7.0_nofiltering_newstrandinfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/mutect2/bwa.mutect2."$i"-T1-DNA1-WGS1.full.vcf.gz "$i; done | sh

### Format and filter PASS INDELs from Berlin cohort 
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -P "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_indels_format_4.0.py INDELS_mutect_indels_patient_"$i"_nofilters_newinfostrands.txt "$i; done | sh

### Generate final file
cat *.FILTERPASS.svtype | sort -n > /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/RESULTS_all_patients_berlin+peifer_cohorts_INDELS_mutect2_PASS_newinfostrands_20.11.09.txt
