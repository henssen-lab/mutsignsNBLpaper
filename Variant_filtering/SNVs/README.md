# Commands to filter and format SNVs from Mutect in the discovery cohort

### Get the SNVs per patient
python script_SNVs_mutect2_6.0_nofiltering_NEW.py bwa.mutect2.tumorsample.vcf.gz sampleID

### Format and filter PASS SNVs per patient
python script_snvs_format_3.0.py SNV_mutect_snv_patient_sampleID_nofilters.txt sampleID
