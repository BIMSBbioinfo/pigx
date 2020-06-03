# PiGx BSseq Pipeline.
#
# Copyright 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
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

# This file was originally taken and modified from
# https://github.com/katwre/makeWGBSnake/blob/master/Rules/Meth_preprocessing_methyldackel_rules.py

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

# ==========================================================================================
# Extract methylation counts with methylDackel
#


def protocol(wc):
    return config["SAMPLES"][wc.sample]['Protocol'].upper()


def keepDups(protocol):
    keepDupsFlag = str(config["general"]["methylation-calling"]["keep-Dups"]).lower()
    keepDups = ""
    if keepDupsFlag == 'auto':
        if (protocol == "RRBS"):
            keepDups = "--keepDupes"
    elif keepDupsFlag in {'true','yes'}:
        keepDups = "--keepDupes"
    
    return keepDups


## TODO: compress methylKit files  
## TODO: make extraction of chg and chh files optional?? 
rule methyldackel_extract_methylKit:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        cpgCallFile = DIR_methcall + "methylDackel/" +"{sample}_methyldackel_CpG.methylKit",
        chgCallFile = DIR_methcall + "methylDackel/" + "{sample}_methyldackel_CHG.methylKit",
        chhCallFile = DIR_methcall + "methylDackel/" + "{sample}_methyldackel_CHH.methylKit"
    wildcard_constraints:
        sample=".+(?<!deduped)"
    params:
        threads = config['execution']['rules']['methyldackel_extract']['threads'],
        prefix = DIR_methcall + "methylDackel/" + "{sample}_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = lambda wc: keepDups(protocol(wc)),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_calls.log"
    message: fmt("Extract methylation calls from bam file using MethylDackel for sample {{sample}} and protocol {params.protocol}")
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "-@ {params.threads}", "{params.keepDups}",
              "--methylKit", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


## TODO: compress methylKit files  
rule methyldackel_extract_methylKit_deduped:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        cpgCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CpG.methylKit",
        chgCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CHG.methylKit",
        chhCallFile = DIR_methcall + "methylDackel/" + \
            "{sample}.deduped_methyldackel_CHH.methylKit"
    params:
        threads = config['execution']['rules']['methyldackel_extract']['threads'],
        prefix = DIR_methcall + "methylDackel/" + "{sample}.deduped_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = lambda wc: keepDups(protocol(wc)),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.deduped.methyldackel_calls.log"
    message: fmt("Extract methylation calls from bam file using MethylDackel for sample {{sample}} and protocol {params.protocol}")
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "-@ {params.threads}", "{params.keepDups}",
              "--methylKit", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


# ==========================================================================================
# Extract methylation bias with methylDackel
#

rule methyldackel_mbias:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        txt = DIR_methcall + "methylDackel/" + "{sample}_mbias_methyldackel.txt",
        p1 = DIR_methcall + "methylDackel/"+ "{sample}_mbias_OB.svg",
        p2 = DIR_methcall + "methylDackel/" + "{sample}_mbias_OT.svg"
    params:
        threads = config['execution']['rules']['methyldackel_extract']['threads'],
        prefix = DIR_methcall + "methylDackel/" + "{sample}_mbias",
        protocol = lambda wc: protocol(wc),
        keepDups = keepDups("{params.protocol}"),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_mbias.log"
    message: fmt("Calculate methylation bias using MethylDackel for sample {{sample}}.")
    shell:
        nice("methyldackel",
             ["mbias", "{input.genome}", "{input.bamfile}",
              "{params.prefix}", "{params.keepDups}",
              "-@ {params.threads}", "--txt",
              "-q {params.minqual}", "> {output.txt}",
              "2> {log}"])


# ==========================================================================================
# Extract methylation bias with methylDackel
#

rule methyldackel_cytosine_report:
    input:
        bamfile = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome = GENOMEFILE
    output:
        DIR_methcall + "methylDackel/" + "{sample}_methyldackel.cytosine_report.txt"
    params:
        threads = config['execution']['rules']['methyldackel_extract']['threads'],
        prefix = DIR_methcall + "methylDackel/" + "{sample}_methyldackel",
        protocol = lambda wc: protocol(wc),
        keepDups = keepDups("{params.protocol}"),
        minqual = int(config['general']
                      ['methylation-calling']['minimum-quality'])
    log:
        DIR_methcall + "{sample}.methyldackel_cytosine_report.log"
    message: fmt("Extract cytosine report from bam file using MethylDackel for sample {{sample}} and context {params.context}")
    shell:
        nice("methyldackel",
             ["extract", "{input.genome}", "{input.bamfile}",
              "-o {params.prefix}", "{params.keepDups}", "-@ {params.threads}",
              "--cytosine_report", "--CHH", "--CHG", "-q {params.minqual}"],
             ("{log}"))


# ==========================================================================================
# Convert to tabix files with methylKit
#

rule tabix_methyldackelfile:
    input:
        DIR_methcall + "methylDackel/" + "{prefix}_methyldackel_{context}.methylKit",
    output:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz",
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz.tbi"
    params:
        sampleid = "{prefix}_{context}",
        assembly = ASSEMBLY,
        treatment = lambda wc: samplesheet(
            wc.prefix.replace(".deduped", ""), 'Treatment'),
        context = "{context}",
        dbdir = DIR_methcall + "methylDackel/" + "/tabix_{context}/",
        mincov = int(config['general']['methylation-calling']
                     ['minimum-coverage'])
    log:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.makeTabix.log"
    message: fmt("Create Tabix file from MethylDackel file for sample {{sample}} and context {params.context}")
    shell:
        nice('Rscript', ["{DIR_scripts}/makeTabix.R",
                         "--location={input}",
                         "--sample.id={params.sampleid}",
                         "--assembly={params.assembly}",
                         "--treatment={params.treatment}",
                         "--context={params.context}",
                         "--mincov={params.mincov}",
                         "--dbdir={params.dbdir}",
                         "--logFile={log}"], "{log}")
