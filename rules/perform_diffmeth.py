# This file was originally taken and modified from
# https://github.com/katwre/makeWGBSnake/blob/master/Rules/Meth_preprocessing_methyldackel_rules.py

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

# ==========================================================================================
# Merge methylation samples

def get_unite_tabixfiles(analysis, tool, context):
    
    samplelist = get_sampleids_from_analysis(analysis)
    protocols = [ samplesheet(sampleid,'Protocol') for sampleid in samplelist]
    dedup_tag = [dedupe_tag(prot) for prot in protocols]
    if tool.lower() == "methylkit":
        paired_flag = [len(samplesheet(sampleid, 'files'))==2 for sampleid in samplelist]
        aligner_tag = ["_1_val_1_bt2.sorted"  if flag else  "_se_bt2.sorted" for flag in paired_flag ]
    else:
        aligner_tag = ["" for sample in samplelist]
    prefix = [x+y+z for x,y,z in zip(samplelist, aligner_tag, dedup_tag)]
    context = context.replace("_destranded","")
    files = []
    for sample in prefix:
        file = os.path.join(DIR_methcall,tool, "tabix_"+context, 
                sample+"_"+context+".txt.bgz")
        files.append(file)
    return(files)

rule unite_meth_calls:
    input:
        samples = lambda wc: get_unite_tabixfiles(wc.analysis,wc.tool,wc.context)     
    output:
        tabixfile = DIR_diffmeth+"{analysis}/methylBase_{analysis}_{context}_{tool}.txt.bgz", 
        tabixindex = DIR_diffmeth+"{analysis}/methylBase_{analysis}_{context}_{tool}.txt.bgz.tbi" 
    params:
        inputfiles = lambda wc, input: ",".join(input.samples),
        sampleids    = lambda wc: ",".join(get_sampleids_from_analysis(wc.analysis)),
        treatments = lambda wc: ",".join([samplesheet(sample,'Treatment') for
            sample in get_sampleids_from_analysis(wc.analysis)]),
        assembly   = ASSEMBLY,
        context    = lambda wc: wc.context.replace("_destranded",""),
        destrand   = lambda wc: True if "destranded" in wc.context else False,
        cores      = int(config['general']['differential-methylation']['cores']),
        outdir     = DIR_diffmeth+"{analysis}/",
        suffix     = "{analysis}_{context}_{tool}"
    log: 
        DIR_diffmeth+"{analysis}/{analysis}_{context}_{tool}_unite.log"
    message: fmt("Uniting samples for differential analysis {wildcards.analysis}")
    shell:
        nice('Rscript',
                ["{DIR_scripts}/methUnite.R",
                    "--inputfiles={params.inputfiles}",
                    "--sampleids={params.sampleids}",
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
        inputfile = os.path.join(DIR_diffmeth,"{analysis}","methylBase_{analysis}_{context}_{tool}.txt.bgz")
    output:
        methylDiff_tabix_file   = os.path.join(DIR_diffmeth,"{analysis}","methylDiff_{analysis}_{context}_{tool}_full.txt.bgz"),
        methylDiff_tabix_index   = os.path.join(DIR_diffmeth,"{analysis}","methylDiff_{analysis}_{context}_{tool}_full.txt.bgz.tbi"),
        results_file   = os.path.join(DIR_diffmeth,"{analysis}","methylDiff_{analysis}_{context}_{tool}_results.tsv"),
    params:
        sampleids   = lambda wc: ','.join(get_sampleids_from_analysis(wc.analysis)),
        treatments = lambda wc: ",".join([samplesheet(sample,'Treatment') for
            sample in get_sampleids_from_analysis(wc.analysis)]),
        assembly    = ASSEMBLY,
        context    = lambda wc: wc.context.replace("_destranded",""),
        destranded = lambda wc: True if "destranded" in wc.context else False,
        treatment_group = lambda wc: config["DManalyses"][wc.analysis]["treatment_sample_groups"],
        control_group = lambda wc: config["DManalyses"][wc.analysis]["control_sample_groups"],
        cores       = int(config['general']['differential-methylation']['cores']),
        methylDiff_results_suffix   = "full",
        outdir     = DIR_diffmeth+"{analysis}/"
    log:
        os.path.join(DIR_diffmeth,"{analysis}","{analysis}_{context}_{tool}_diffmeth.log")
    message: fmt("Calculating differential methylation for analysis {wildcards.analysis}")
    shell:
        nice('Rscript', 
                ['{DIR_scripts}/methDiff.R',
                    '--inputfile={input.inputfile}',
                    '--sampleids={params.sampleids}',
                    '--treatments={params.treatments}',
                    '--assembly={params.assembly}',
                    '--context={params.context}',
                    '--destranded={params.destranded}',
                    '--treatment_group={params.treatment_group}',
                    '--control_group={params.control_group}',
                    '--cores={params.cores}',
                    '--methylDiff_results_suffix={params.methylDiff_results_suffix}',
                    '--resultsFile={output.results_file}',
                    "--outdir={params.outdir}",
                    '--logFile={log}'
                    ],"{log}")


# ==========================================================================================
# Generate the final report for differential methylation between pairs of analysis values:

rule diffmeth_report:
    input:
        methylBase_tabix_file   = os.path.join(DIR_diffmeth,"{analysis}","methylBase_{analysis}_{context}_{tool}.txt.bgz"),
        methylDiff_tabix_file   = os.path.join(DIR_diffmeth,"{analysis}","methylDiff_{analysis}_{context}_{tool}_full.txt.bgz"),
        methylDiff_results_file   = os.path.join(DIR_diffmeth,"{analysis}","methylDiff_{analysis}_{context}_{tool}_results.tsv"),
        template           = os.path.join(DIR_templates,"diffmeth.Rmd"),
        chrom_seqlengths   = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv")
    output:
        report        = os.path.join(DIR_final, "{analysis}", "{analysis}_{context}_{tool}.diffmeth-report.html")
    params:
        # specific for this report
        sampleids              = lambda wc: ','.join(get_sampleids_from_analysis(wc.analysis)),
        treatments             = lambda wc: ",".join([samplesheet(sample,'Treatment') for
            sample in get_sampleids_from_analysis(wc.analysis)]),
        assembly               = ASSEMBLY,
        context                = lambda wc: wc.context.replace("_destranded",""),
        destranded             = lambda wc: True if "destranded" in wc.context else False,
        treatment_group        = lambda wc: config["DManalyses"][wc.analysis]["treatment_sample_groups"],
        control_group          = lambda wc: config["DManalyses"][wc.analysis]["control_sample_groups"],
        scripts_dir            = DIR_scripts,
        cpgIsland_bedfile      = CPGISLAND_BEDFILE,
        refGenes_bedfile       = REFGENES_BEDFILE,
        chrom_seqlengths       = os.path.join(DIR_mapped,"Refgen_"+ASSEMBLY+"_chromlengths.csv"),
        qvalue                 = float(config['general']['differential-methylation']['qvalue']),
        difference             = float(config['general']['differential-methylation']['difference']),
        webfetch               = config['general']['differential-methylation']['annotation']['webfetch'],
        # required for any report
        bibTexFile             = BIBTEXPATH,
        prefix                 = "{analysis}_{context}_{tool}",
        workdir                = os.path.join(DIR_final,"{analysis}"),
        logo                   = LOGOPATH
    log:
        os.path.join(DIR_final,"{analysis}","{analysis}_{context}_{tool}_diffmeth-report.log")
    message: fmt("Compiling differential methylation report " + "for analysis " + "{wildcards.analysis}")
    # run:
    #     generateReport(input, output, params, log, "")
    shell:
        nice('Rscript', ["{DIR_scripts}/generate_report.R",
                           "--reportFile={input.template}",
                           "--outFile={output.report}",
                           "--workdir={params.workdir}",
                           "--logo={params.logo}",
                           "--bibTexFile={params.bibTexFile}",
                           "--prefix={params.prefix}",
                           "--report.params='{{"+
                           ",".join([
                               '"sampleids":"{params.sampleids}"',
                               '"treatments":"{params.treatments}"',
                               '"assembly":"{params.assembly}"',
                               '"context":"{params.context}"',
                               '"destranded":"{params.destranded}"',
                               
                               '"treatment_group":"{params.treatment_group}"',
                               '"control_group":"{params.control_group}"',
                               
                               '"methylBase_file":"{input.methylBase_tabix_file}"',
                               '"methylDiff_file":"{input.methylDiff_tabix_file}"',
                               '"methylDiff_results_file":"{input.methylDiff_results_file}"',
                               
                               '"scripts_dir":"{params.scripts_dir}"',
                               '"cpgIsland_bedfile":"{params.cpgIsland_bedfile}"',
                               '"refGenes_bedfile":"{params.refGenes_bedfile}"',
                               '"chrom_seqlengths":"{params.chrom_seqlengths}"',
                               '"qvalue":"{params.qvalue}"',
                               '"difference":"{params.difference}"',
                               '"webfetch":"{params.webfetch}"'
                           ])+"}}'",
                           "--logFile={log}"], "{log}", "echo 'something went wrong' ")
