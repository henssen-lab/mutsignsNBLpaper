# Commands to filter and format SNVs in the discovery cohort (Peifer+Berlin)

### Get the SNVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -P "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_SNVs_mutect2_6.0_nofiltering_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/mutect2/bwa.mutect2."$i"-T1-DNA1-WGS1.full.vcf.gz "$i; done | sh

### Format and filter PASS SNVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -P "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_snvs_format_3.0.py SNV_mutect_snv_patient_"$i"_nofilters.txt "$i; done | sh

### Get the SNVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -P "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_SNVs_mutect2_6.0_nofiltering_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/mutect2/bwa.mutect2."$i"-T1-DNA1-WGS1.full.vcf.gz "$i; done | sh

### Format and filter PASS SNVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -P "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_snvs_format_3.0.py SNV_mutect_snv_patient_"$i"_nofilters.txt "$i; done | sh

### Generate final file
cat *.FILTERPASS.svtype | sort -n > /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/RESULTS_all_patients_berlin+peifer_cohorts_SNVS_mutect_PASS_20.11.09.txt

