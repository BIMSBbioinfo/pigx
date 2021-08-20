#!/bin/bash

while IFS=, read -r col1 col2 col3
do
           newname=${col2//[[:blank:]]/}

     mv *${col1}*R1* ${newname%%[[:cntrl:]]}_R1.fastq.gz
     mv *${col1}*R2* ${newname%%[[:cntrl:]]}_R2.fastq.gz
     done < < $(tail -n +2 $1)