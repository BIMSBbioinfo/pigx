# This file was originally taken and modified from
# https://github.com/katwre/makeWGBSnake/blob/master/Rules/Meth_preprocessing_methyldackel_rules.py

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

# ==========================================================================================
# Merge methylation samples

def get_unite_tabixfiles(treatment, tool, context):
    samplelist =  get_sampleids_from_treatment(treatment)
    protocols = [ samples(sampleid,'Protocol') for sampleid in samplelist]
    dedup_tag = [dedupe_tag(prot) for prot in protocols]
    if tool.lower() == "methylkit":
        paired_flag = [len(samples(sampleid, 'files'))==2 for sampleid in samplelist]
        aligner_tag = ["_1_val_1_bt2.sorted"  if flag else  "_se_bt2.sorted" for flag in paired_flag ]
    else:
        aligner_tag = ["" for sample in samplelist]
    prefix = [x+y+z for x,y,z in zip(samplelist, aligner_tag, dedup_tag)]
    files = []
    for sample in prefix:
        file = os.path.join(DIR_methcall,tool, "tabix_"+context, 
                sample+"_"+context+".txt.bgz")
        files.append(file)
    return(files)

rule unite_meth_calls:
    input:
        samples = lambda wc: get_unite_tabixfiles(wc.treatment,wc.tool,wc.context)     
    output:
        tabixfile = DIR_diffmeth+"{treatment}/methylBase_{treatment}_{context}_{tool}.txt.bgz", 
        tabixindex = DIR_diffmeth+"{treatment}/methylBase_{treatment}_{context}_{tool}.txt.bgz.tbi" 
    params:
        inputfiles = lambda wc, input: ",".join(input.samples),
        samples    = lambda wc: ",".join(get_sampleids_from_treatment(wc.treatment)),
        treatments = lambda wc: ",".join([samples(sample,'Treatment') for
            sample in get_sampleids_from_treatment(wc.treatment)]),
        assembly   = ASSEMBLY,
        context    = "{context}",
        destrand   = lambda wc: destrand(wc.context),
        cores      = int(config['general']['differential-methylation']['cores']),
        outdir     = DIR_diffmeth+"{treatment}/",
        suffix     = "{treatment}_{context}_{tool}"
    log: 
        DIR_diffmeth+"{treatment}/{treatment}_{context}_{tool}_unite.log"
    shell:
        nice('Rscript',
                ["{DIR_scripts}/methUnite.R",
                    "--inputfiles={params.inputfiles}",
                    "--samples={params.samples}",
                    "--treatments={params.treatments}",
                    "--assembly={params.assembly}",
                    "--context={params.context}",
                    "--destrand={params.destrand}",
                    "--cores={params.cores}",
                    "--outdir={params.outdir}",
                    "--suffix={params.suffix}",
                    "--logFile={log}"],"{log}")



# ==========================================================================================
# Perform differential methylation analysis:

rule diffmeth:
    ## paths inside input and output should be relative
    input:
        inputfile = DIR_diffmeth+"{treatment}/methylBase_{treatment}_{context}_{tool}.txt.bgz"  
    output:
        methylDiff_tabix_file   = os.path.join(DIR_diffmeth,"{treatment}","methylDiff_{treatment}_{context}_{tool}_full.txt.bgz"),
        methylDiff_tabix_index   = os.path.join(DIR_diffmeth,"{treatment}","methylDiff_{treatment}_{context}_{tool}_full.txt.bgz.tbi"),
        results_file   = os.path.join(DIR_diffmeth,"{treatment}","methylDiff_{treatment}_{context}_{tool}_results.tsv"),
    params:
        sampleids   = lambda wc: ','.join(get_sampleids_from_treatment(wc.treatment)),
        treatments = lambda wc: ",".join([samples(sample,'Treatment') for 
            sample in get_sampleids_from_treatment(wc.treatment)]),
        assembly    = ASSEMBLY,
        context    = "{context}",
        destranded = lambda wc: destrand(wc.context),
        cores       = int(config['general']['differential-methylation']['cores']),
        methylDiff_results_suffix   = "full",
        outdir     = DIR_diffmeth+"{treatment}/"
    log:
        os.path.join(DIR_diffmeth+"{treatment}_{context}_{tool}_diffmeth.log")
    message: fmt("Calculating differential methylation.")
    shell:
        nice('Rscript', 
                ['{DIR_scripts}/methDiff.R',
                    '--inputfile={input.inputfile}',
                    '--sampleids={params.sampleids}',
                    '--treatments={params.treatments}',
                    '--assembly={params.assembly}',
                    '--context={params.context}',
                    '--destranded={params.destranded}',
                    '--cores={params.cores}',
                    '--methylDiff_results_suffix={params.methylDiff_results_suffix}',
                    '--resultsFile={output.results_file}',
                    "--outdir={params.outdir}",
                    '--logFile={log}'
                    ])

