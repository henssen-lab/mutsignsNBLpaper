# Commands to filter and format INDEL calls from Mutect in the discovery cohort

### Get the INDELs per patient
python script_INDELS_mutect2_7.0_nofiltering_newstrandinfo_NEW.py bwa.mutect2.tumorsample.vcf.gz sampleID

### Format and filter PASS per patient
python script_indels_format_4.0.py INDELS_mutect_indels_patient_sampleID_nofilters_newinfostrands.txt sampleID

### Generate sorted final file for ALL patients
cat *.FILTERPASS.svtype | sort -n > RESULTS_all_patients_berlin+peifer_cohorts_INDELS_mutect2_PASS_newinfostrands_20.11.09.txt
