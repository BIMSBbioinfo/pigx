## workflow

| Steps | files/commands used | notes |
|:------|:-----------|:------|  
|prepare sample space for v-pipe | prep_sample_space_210320.py|makes a sample list from a directory full with sample.fasta.gz|  
|prepare sample.tsvs | either by hand or will be done automatically by v-pipe| vpipe-version should only be used if there are not multiple sample dirs in the sample dir and if the readlength is always >250bp|  
|run v-pipe | vpipe.config | no visualization, with `--until lofreq` specified to skip `rule snv` (the haplotyping/shorah step)|  
|make lofreq csv | vcfTocsv.py | to use it for the reports|  
|download sars_cov_2 database for vep| * `wget ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz`  * `mv sars_cov_2*.tar.gz ~/.vep && cd ~/.vep && tar xf sars_cov_2*.tar.gz`|only once|  
|run ensemble vep and parse the file | ensemble_vep.py | (so far I only did this manually using the online tool but that file should contain everything we need to do it properly)|  
|make the sample report | run_variant_report.R, variantreport_p_sample.rmd| command for the run-script:   `Rscript <path/to/run_variant_report.R> --reportFile=<path/to/variantreport_p_sample.rmd> --vep_txt_file=<path/to/vpipe/workingdir/samples/sample_parent_folder/sample/vep_reports/vep_sarscov2_[sample].txt> --snv_csv_file=<path/to/vpipe/workingdir/samples/sample_parent_folder/sample/variants/SNVs/snvs.csv> --location_sigmuts="<path/to/folder/with/outbreak.info/signature_mutations> --sample_dir=<path/to/vpipe/workingdir/samples/sample_parent_folder/sample> --sample_name="<sample>" ` ! The report works only with the output of the online Vep so far, I'll rework the parsing script for the CLI output|  
 
### planned
* additional small skript that runs various samtools commands within the sample/alignment folder and output a txt file with 
some stats information, this one should be included in the reports too
  