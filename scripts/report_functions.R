# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017, 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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
#
## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included
render2HTML <- function(reportFile,
                            outFile,
                            outDir,
                            finalReportDir,
                            report.params=NULL,
                            self.contained=TRUE,
                            logFile = NULL)
{
  
  output_format = rmarkdown::all_output_formats(reportFile, 'UTF-8')
  
  if(!dir.exists(finalReportDir)) dir.create(finalReportDir, recursive = TRUE)
  
  if(is.null(report.params)) report.params <- list()
  
  ## render single report
  message("render single report")
  
  ## make independent intermediate dirs
  interDir <- paste0(outDir,"/",outFile,"_tmp")
  
  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    intermediates_dir = interDir,
    output_file = outFile,
    knit_root_dir = outDir,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = FALSE,
      code_folding = "hide",
      self_contained = self.contained,
      includes = list(in_header = "pigx_bsseq_logo.html"),
      bibliography= "reports.bib"
    ),
    params=report.params,
    quiet = FALSE,
    clean = TRUE
  )

  
  on.exit(unlink(interDir,recursive = TRUE),add = TRUE)
  
  
}

