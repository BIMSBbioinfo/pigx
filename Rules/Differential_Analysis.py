# ---------------------------------------------------------------------------- #
def get_differential_analysis_bedfile(wc):
    analysisDict = config['differential_analysis'][wc.analysis]
    if ('Peakset' in analysisDict) and ( analysisDict['Peakset']):
        infiles = os.path.join(PATH_RDS_FEATURE,
                wc.analysis + '_FeatureCombination.bed'),
    else:
        infiles = os.path.join(PATH_PEAK,  
                wc.analysis, 
                wc.analysis + "_peaks.narrowPeak" )
    return(infiles)

def get_differential_analysis_countfiles(wc):
    analysisDict = config['differential_analysis'][wc.analysis]
    infiles = expand(
            os.path.join(
                PATH_RDS_COUNTS,
                '{analysis}_{name}_FeatureCounts.rds'), 
            analysis = wc.analysis, 
            name = analysisDict['Case'] + analysisDict['Control']
            )
    print(infiles)
    return(infiles)

# ----------------------------------------------------------------------------- #
rule feature_counting:
    input:
        bedfile  = lambda wc: get_differential_analysis_bedfile(wc),
        bamfile = os.path.join(PATH_MAPPED, GENOME_TYPES['Main'], "{name}", "{name}.sorted.bam")
    output:
        outfile  = os.path.join(PATH_RDS_COUNTS,'{analysis}_{name}_FeatureCounts.rds')
    params:
        threads   = config['execution']['rules']['feature_counting']['threads'],
        name = "{analysis}", 
        params_tool  = PARAMS['feature_counting'],
        scriptdir = SCRIPT_PATH,
        Rscript   = SOFTWARE['Rscript']['executable']
    log:
        logfile = os.path.join(PATH_LOG, '{analysis}_{name}_feature_counting.log')
    message:
        """
            Running: feature_counting:
                features:   {input.bedfile}
                sample:     {input.bamfile}
                output:     {output.outfile}
            """
    run:
        if (params.name in set(CUSTOM_PARAMS.keys())):
            if CUSTOM_PARAMS[params.name]:
                if 'feature_counting' in CUSTOM_PARAMS[params.name]:
                    params.params_tool.update(CUSTOM_PARAMS[params.name]['feature_counting'])
        RunRscript(input, output, params, log.logfile, 'Feature_Counting.R')

# ----------------------------------------------------------------------------- #
rule merge_counts:
    input:
        countfiles = lambda wc: get_differential_analysis_countfiles(wc)
    output:
        outfile  = os.path.join(PATH_REPORTS,'{analysis}','{analysis}_FeatureCounts.tsv'),
        statfile  = os.path.join(PATH_REPORTS,'{analysis}','{analysis}_FeatureCounts_Stats.tsv')
    params:
        scriptdir = SCRIPT_PATH,
        Rscript   = SOFTWARE['Rscript']['executable']
    log:
        logfile = os.path.join(PATH_LOG, '{analysis}_merge_counts.log')
    message:
        """
            Running: merge_counts:
                features:   {input}
                output:     {output.outfile}
            """
    run:
        RunRscript(input, output, params, log.logfile, 'Merge_Counts.R')
