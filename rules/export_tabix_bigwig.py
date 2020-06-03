
# PiGx BSseq Pipeline.
#
# Copyright © 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
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

# ==========================================================================================
# Export methylation from tabix to bigwig
#

def destrand(context):
    return config['general']['export-bigwig']['context'][context.lower()]['destrand']

rule export_tabix_bigwig:
    input:
        seqlengths = os.path.join(
            DIR_mapped, ASSEMBLY + "_chromlengths.csv"),
        filepath = os.path.join(
            DIR_methcall, "{tool}", "tabix_{context}/{prefix}_{context}.txt.bgz")
    output:
        bw = os.path.join(DIR_bigwig, "{prefix}.{context}_{tool}.bw")
    wildcard_constraints:
        context=".+(?<!destranded)"
    params:
        assembly = ASSEMBLY,
        destrand = lambda wc: destrand(wc.context)
    log:
        DIR_bigwig + "{prefix}.{context}.{tool}.export_tbx2bw.log"
    message: fmt("exporting methylation as bigwig for sample {wildcards.prefix}.")
    shell:
        nice('Rscript', ["{DIR_scripts}/export_tbx2bw.R",
                         "--filepath={input.filepath}",
                         "--seqlengths_path={input.seqlengths}",
                         "--assembly={params.assembly}",
                         "--destrand={params.destrand}",
                         "--out_path={output}",
                         "--logFile={log}"], "{log}")

rule export_tabix_bigwig_destrand:
    input:
        seqlengths = os.path.join(
            DIR_mapped, ASSEMBLY + "_chromlengths.csv"),
        filepath = os.path.join(
            DIR_methcall, "{tool}", "tabix_{context}/{prefix}_{context}.txt.bgz")
    output:
        bw = os.path.join(DIR_bigwig, "{prefix}.{context}_destranded_{tool}.bw")
    params:
        assembly = ASSEMBLY,
        destrand = lambda wc: destrand(wc.context)
    log:
        DIR_bigwig + "{prefix}.{context}.{tool}.export_tbx2bw.log"
    message: fmt("exporting methylation as bigwig for sample {wildcards.prefix}.")
    shell:
        nice('Rscript', ["{DIR_scripts}/export_tbx2bw.R",
                         "--filepath={input.filepath}",
                         "--seqlengths_path={input.seqlengths}",
                         "--assembly={params.assembly}",
                         "--destrand={params.destrand}",
                         "--out_path={output}",
                         "--logFile={log}"], "{log}")
