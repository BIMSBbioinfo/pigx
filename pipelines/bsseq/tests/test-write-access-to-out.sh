#!/bin/sh

set -e
set -u

test_name=`basename $0`
dir=$(mktemp -d "/tmp/pigx.XXXX")
mkdir -p ${dir}/${test_name}
cd ${dir}/${test_name}

# create test data
mkdir -p in
mkdir -p genome
touch genome/foo.fasta

cat <<EOF > settings.yaml
locations:
  input-dir: in/
  output-dir: out/
  genome-dir: genome/
general:
  differential-methylation:
    treatment-groups: []
EOF
cat <<EOF > samplesheet.csv
Read1,Read2,SampleID,ReadType,Treatment
EOF

# run the pipeline
export PIGX_UGLY=1
export PIGX_UNINSTALLED=1

# This should fail, because the genome directory does not contain
# bisulfite genome files, and we have no write access to the genome
# directory.
chmod -w genome

${builddir}/pigx-bsseq --target=genome-prep -s $PWD/settings.yaml $PWD/samplesheet.csv 2>&1 >/dev/null | grep "TODO"
rm -rf ${dir}
