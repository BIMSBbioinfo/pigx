report_scripts_dir=$1
sample_sheet_file=$2
cut_sites_file=$3
pigx_output_dir=$4
target_name=$5
site_dir=$6
RSCRIPT_PATH=$7

# copy and rename site files to target folder
#also change @PARAMS_TARGET_NAME@ line in yaml header to the target_name: {target_name}
cp ${report_scripts_dir}/index.Rmd ${site_dir}
ls ${report_scripts_dir}/*.Rmd | grep -v 'index.Rmd' | \\
    while read f; do bs=`basename ${f}`; cat $f | sed "s/@PARAMS_TARGET_NAME@:/target_name: ${target_name}/" > ${site_dir}/"${target_name}.${bs}"; done

# create a config.yml file that will be shared across all report scripts
config_file="${site_dir}/config.yml"

echo "sample_sheet: ${sample_sheet_file}" > ${config_file}
echo "cut_sites_file: ${cut_sites_file}" >> ${config_file}
echo "pigx_output_dir: ${pigx_output_dir}" >> ${config_file}


# render the report website
echo "Rendering website at ${site_dir}"
${RSCRIPT_PATH} -e "library(rmarkdown); rmarkdown::render_site(\"${site_dir}\")"
