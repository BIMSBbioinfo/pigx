# PiGx ChIPseq Pipeline.
#
# Copyright © 2020 Alexander Blume <alexander.blume@mdc-berlin.de>
#
# This file is part of the PiGx ChIPseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Set of Rules to prepare input for differential analysis report.


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
    groups = flatten([analysisDict['Case']] + [analysisDict['Control']])
    samps = [lookup(
                STRUCTURE_VARIABLES['SAMPLE_SHEET_GROUP_NAME'],
                sample,
                ['SampleName']
                ) for sample in groups]
    samps = flatten(samps)
    infiles = expand(
            os.path.join(
                PATH_RDS_COUNTS,
                '{analysis}_{name}_FeatureCounts.rds'), 
            analysis = wc.analysis, 
            name = samps

            )
    return(infiles)

def get_differential_analysis_settings(wc, setting):
    analysisDict = config['differential_analysis'][wc.analysis]
    if setting in analysisDict.keys():
        value = ','.join(flatten([analysisDict[setting]]))
    else:
        value = ''
    return(value)

# ----------------------------------------------------------------------------- #
rule feature_counting:
    input:
        bedfile  = lambda wc: get_differential_analysis_bedfile(wc),
        bamfile = os.path.join(PATH_MAPPED, GENOME_TYPES['Main'], "{name}", "{name}.sorted.bam")
    output:
        outfile  = os.path.join(PATH_RDS_COUNTS,'{analysis}_{name}_FeatureCounts.rds')
    wildcard_constraints:
        analysis = DIFF_ANN_CONSTRAINT
    params:
        threads   = config['execution']['rules']['feature_counting']['threads'],
        name = "{analysis}", 
        params_tool  = PARAMS['feature_counting'],
        library_type = lambda wc: get_library_type(wc.name),
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
    wildcard_constraints:
        analysis = DIFF_ANN_CONSTRAINT
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

#----------------------------------------------------------------------------- #
rule translate_sample_sheet_for_report:
    input: 
        SAMPLE_SHEET_FILE
    output: 
        outfile = os.path.join(PATH_REPORTS, "colData.tsv")
    run:
        import csv
        with open(output.outfile, 'w', newline='') as file:
            writer = csv.writer(file, delimiter = "\t")
            header = list(set(SAMPLE_SHEET[0].keys()) - 
                    set(STRUCTURE_VARIABLES['SAMPLE_SHEET_COLUMN_NAMES']+
                        [STRUCTURE_VARIABLES['SAMPLE_SHEET_GROUP_NAME']]))  
            header =[STRUCTURE_VARIABLES['SAMPLE_SHEET_GROUP_NAME']] + header        
            writer.writerow([''] + header)
            for sample in SAMPLE_NAMES:
                writer.writerow(lookup('SampleName', sample, ['SampleName'] + header))

#----------------------------------------------------------------------------- #
rule knit_differential_analysis_report:
    input:
        countDataFile = os.path.join(PATH_REPORTS,'{analysis}','{analysis}_FeatureCounts.tsv'),
        colDataFile   = os.path.join(PATH_REPORTS, "colData.tsv"),
        gtfFile       = rules.prepare_annotation.output.outfile
    output:
        outfile    = os.path.join(PATH_REPORTS,'{analysis}','{analysis}_DeseqReport.html')
    wildcard_constraints:
        analysis = DIFF_ANN_CONSTRAINT
    params:
        report_template     = REPORT_DA_TEMPLATE,
        peakFile            = lambda wc: get_differential_analysis_bedfile(wc),
        caseSampleGroups    = lambda wc: get_differential_analysis_settings(wc,'Case'),
        controlSampleGroups = lambda wc: get_differential_analysis_settings(wc,'Control'),
        covariates          = lambda wc: get_differential_analysis_settings(wc,'Covariates'),
        groupCol            = STRUCTURE_VARIABLES['SAMPLE_SHEET_GROUP_NAME'],
        prefix              = '{analysis}', 
        workdir             = os.path.join(PATH_REPORTS, '{analysis}'),
        scriptdir           = SCRIPT_PATH,
        organism            = config['general']['organism'] if 'organism' in config['general'] else '',
        logo                = LOGO_PATH 
    log:
        logfile = os.path.join(PATH_LOG, '{analysis}_knit_deseq_report.log')
    message:"""
            Running: knit_differential_analysis_report:
                Analysis:   {wildcards.analysis}
                output:     {output.outfile}
        """
    run:
        RunRscript(input, output, params, log.logfile, 'Knit_Report.R')
