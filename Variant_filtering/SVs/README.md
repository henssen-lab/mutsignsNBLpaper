# Commands to filter, merge, and format SVs in the discovery cohort (Peifer+Berlin)

Collapse SVs (>50bp) from different callers (Svaba, Delly and Novobreak)
Filter them at least 2 callers except insertions
Including SVs from indels from Svaba and insertions PASS from delly2

## INTRACHROMOSOMAL SVs

### Get the intrachromosomal SVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -e "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_merge_intraSVs_svaba_svabafromindels_delly2_novobreak_5.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/svaba/svaba-sv.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/svaba/svaba-sv-fromindels.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/delly2/delly2-svs.pre.ca-0.05.vcf /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/novobreak_results_hg19_joern/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt "$i; done | sh

### Format and filter the intrachromosomal SVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -e "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_filter_atleast2callers_4.0_withinsertions.py merge_intraSVs_diffvariantcallers_svaba_svaba_fromindels_delly2_novobreak_patient_"$i"_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt "$i; done | sh

### Get the intrachromosomal SVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -e "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_merge_intraSVs_svaba_svabafromindels_delly2_novobreak_5.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/svaba/svaba-sv.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/svaba/svaba-sv-fromindels.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/delly2/delly2-svs.pre.ca-0.05.vcf /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/novobreak_results_hg19_joern/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt "$i; done | sh

### Format and filter the intrachromosomal SVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -e "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_filter_atleast2callers_4.0_withinsertions.py merge_intraSVs_diffvariantcallers_svaba_svaba_fromindels_delly2_novobreak_patient_"$i"_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt "$i; done | sh

### Generate final file
cat *.FILTER2CALLERS.svtype | sort -n > /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/RESULTS_all_patients_berlin+peifer_cohorts_intraSVs_merge_svaba_svabafromindels_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt


## INTERCHROMOSOMAL SVs

### Get the interchromosomal SVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -e "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_merge_translocations_svaba_delly2_novobreak_4.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/svaba/svaba-sv.PASS.vcf.gz //fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC/"$i"/delly2/delly2-svs.pre.ca-0.05.vcf /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/novobreak_results_hg19_joern/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt "$i; done | sh

### Format and filter the interchromosomal SVs from Peifer cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Peifer_cohort_VC | grep -e "NB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_filter_atleast2callers_4.0_translocations.py merge_translocations_diffvariantcallers_svaba_delly2_novobreak_patient_"$i"_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt "$i; done | sh

### Get the interchromosomal SVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -e "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_merge_translocations_svaba_delly2_novobreak_4.0_collapsedups_nofiltering_newstrandinfo_precisioninfo_NEW.py /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/svaba/svaba-sv.PASS.vcf.gz /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC/"$i"/delly2/delly2-svs.pre.ca-0.05.vcf /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/novobreak_results_hg19_joern/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt "$i; done | sh

### Format and filter the interchromosomal SVs from Berlin cohort
for i in `ls /fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/results_Berlin_cohort_VC | grep -e "CB"`; do echo "python /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/script_filter_atleast2callers_4.0_translocations.py merge_translocations_diffvariantcallers_svaba_delly2_novobreak_patient_"$i"_noVCFdups_collapsedups_nofilters_newinfostrands_precisioninfo.txt "$i; done | sh

### Generate final file
cat *.FILTER2CALLERS.svtype | sort -n > /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/RESULTS_all_patients_berlin+peifer_cohorts_translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt
