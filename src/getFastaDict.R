# crispr-DART pipeline
#
# Copyright © 2017-2020 Bora Uyar <bora.uyar@mdc-berlin.de>
#
# This file is part of the crispr-DART pipeline
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

library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

fastaFile <- args[1]
outFile <- args[2]

seq <- Biostrings::readDNAStringSet(fastaFile)

# print a header line in the style of sam file specification
# the line should contain the length of the amplicon and its name
# e.g. "@SQ	SN:dpy-10	LN:1320"
write(x = paste0("@SQ\t", 
                 "SN:", names(seq), "\t",
                 "LN:", width(seq)), 
      file = outFile)
