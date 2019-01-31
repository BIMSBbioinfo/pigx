report_scripts_dir=$1
sample_sheet_file=$2
cut_sites_file=$3
pigx_output_dir=$4
site_dir=$5
RSCRIPT_PATH=$6

# copy site files to target folder
cp ${report_scripts_dir}/* ${site_dir}

# create a config.yml file that will be shared across all report scripts
config_file="${site_dir}/config.yml"

echo "sample_sheet: ${sample_sheet_file}" > ${config_file}
echo "cut_sites_file: ${cut_sites_file}" >> ${config_file}
echo "pigx_output_dir: ${pigx_output_dir}" >> ${config_file}


# render the report website
echo "Rendering website at ${site_dir}"
${RSCRIPT_PATH} -e "library(rmarkdown); rmarkdown::render_site(\"${site_dir}\")"
