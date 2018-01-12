# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <b.osberg@tum.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
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

args <- commandArgs(trailingOnly = TRUE)
logFile     <- args[1]
refgenesLoc <- args[2]
assembly    <- args[3]
scripts_dir <- args[4]
genome_dir  <- args[5]

source(paste0(scripts_dir, "fetch_procedures.R"))

## catch output and messages into log file
out <- file(logFile, open = "wt")
sink(out, type = "output")
sink(out, type = "message")

lookupBedFile(type = "refGene", filename = refgenesLoc, dir = genome_dir, assembly = assembly)
