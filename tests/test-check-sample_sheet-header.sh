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

cat <<EOF > settings.yaml
locations:
  input-dir: in/
  output-dir: out/
  genome-dir: genome/
EOF
cat <<EOF > sample_sheet.csv
foo,bar,baz
EOF

# run the pipeline
export PIGX_UGLY=1
export PIGX_UNINSTALLED=1

${builddir}/pigx-bsseq -s settings.yaml sample_sheet.csv 2>&1 | grep "First columns of the input table have to be"
rm -rf ${dir}
