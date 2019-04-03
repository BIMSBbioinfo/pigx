# ---------------------------------------------------------------------------- #
import subprocess
import re
import os
import sys
import itertools

# ---------------------------------------------------------------------------- #
# uses global variable SAMPLE_SHEET
def lookup(column, predicate, fields=[]):
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column]==predicate]
    return [record[field] for record in records for field in fields]

# ---------------------------------------------------------------------------- #
# given a sample name returns library type, depending on number of fastq files
# uses global variable SAMPLE_SHEET
def get_library_type(name):
    library_type = lookup('SampleName', name, ['library_type'])
    return(library_type[0])

# ---------------------------------------------------------------------------- #
# given a sample name returns whether the sample was spiked
# uses global variable SAMPLE_SHEET
# return No by default
def get_spikein_information(name):

    return_value = 'No'
    if 'Spike-in' in set(SAMPLE_SHEET[0].keys()):
        return_value = lookup('SampleName', name, ['Spike-in'])[0]

    return return_value


# ---------------------------------------------------------------------------- #
# function which reads the sample sheet from a csv/xlsx file
def read_SAMPLE_SHEET(config):

    # Checks whether the sample sheet file exists
    if config['locations']['sample-sheet'] == None:
        message = 'ERROR: sample-sheet file is not defined\n'
        sys.exit(message)
    elif not os.path.isfile(config['locations']['sample-sheet']):
        print(config['locations']['sample-sheet'])
        message = 'ERROR: sample-sheet file does not exist\n'
        sys.exit(message)

    SAMPLE_SHEET_FILE = config['locations']['sample-sheet']
    # Check for the allowed extensions
    if not SAMPLE_SHEET_FILE.endswith(('.csv','.xlsx','.xls')):
        sys.exit('ERROR: File format of the sample_sheet has to be: csv or excel table (xls,xlsx).\n')
    ## Load sample sheet
    # either from excel file
    if SAMPLE_SHEET_FILE.endswith(('.xlsx','.xls')):
        with xlrd.open_workbook(SAMPLE_SHEET_FILE) as book:
            # assume that the first book is the sample sheet
            sheet = book.sheet_by_index(0)
            rows = [sheet.row_values(r) for r in range(0, sheet.nrows)]
    # or from csv file
    elif SAMPLE_SHEET_FILE.endswith('.csv'):
        with open(SAMPLE_SHEET_FILE, 'r') as fp:
            rows = [row for row in csv.reader(fp, delimiter=',')]
    # in both cases we need to strip leading or trailing whitespaces
    rows = [list(map(str.strip, row)) for row in rows]
    header = rows[0]; rows = rows[1:]
    # we do not want duplicate rows in the sample sheet
    rows.sort()
    rows = list(row for row,_ in itertools.groupby(rows))
    # then we create a dictionary from header and rows
    SAMPLE_SHEET = [OrderedDict(zip(header, row)) for row in rows if row]

    # Goes through the sample sheet, and defines the library type based on
    # the existence of one or two Read input files
    for input_sample in SAMPLE_SHEET:
        files = [input_sample['Read'], input_sample['Read2']]
        files = [file for file in files if file]
        if len(files) == 2:
            library_type = "paired"
        else:
            library_type = "single"

        input_sample['library_type'] = library_type

    return(SAMPLE_SHEET)


# ---------------------------------------------------------------------------- #
# given a sample name returns fastq location(s)
# uses global variable SAMPLE_SHEET
def get_fastq_input(name,prefix):
    samps = TRIM_GALORE_DICT[name][prefix]['raw']
    infiles = [os.path.join(PATH_FASTQ, i) for i in flatten(samps) if i]
    return(infiles)

# given a fastq file strips of either of the following extensions:
# "fastq.gz","fastq","fq.gz","fq"
# this is required for parsing trim galore output
def replace_fastq_ext(fqfile,replacement):
    p = re.compile('.f(ast)*q(\.gz)*')
    repl = p.sub(replacement,fqfile)
    return(repl)

# ---------------------------------------------------------------------------- #
# The function generates the locations of the genome link, genome prefix, and bowtie2
# index files for the given genome
# it is used for the main and spike-in genomes
# CAUTION: GENOME_HASH is a global dict
def generate_genome_files(genome_location, index_path, genome_type, genome_name=None):

    GENOME_HASH[genome_type] = {}
    if genome_name == None:
        genome_name = os.path.basename(GENOME_LOCATION)

    GENOME_HASH[genome_type]['genome_location'] = genome_location
    GENOME_HASH[genome_type]['genome_name']   = genome_name
    GENOME_HASH[genome_type]['genome_prefix'] = os.path.join(index_path, genome_type, GENOME_HASH[genome_type]['genome_name'])
    GENOME_HASH[genome_type]['genome_link'] = GENOME_HASH[genome_type]['genome_prefix'] + '.fa'
    GENOME_HASH[genome_type]['bowtie_index'] = GENOME_HASH[genome_type]['genome_prefix'] + '.1.bt2'

# ---------------------------------------------------------------------------- #
# given a sample name returns fastq location(s)
def get_trimming_dict(name):
    fqfiles = [file for file in lookup('SampleName',name,['Read','Read2']) if file]
    lib_type = get_library_type(name)
    trimming_dict = dict()
    # defaults resemble single end case
    prefix = name
    num_reps = len(fqfiles)
    read_ext = ['R']
    # and get updated if library is paired end
    if ( (len(fqfiles) % 2) == 0 and lib_type == 'paired'):
        num_reps = int(len(fqfiles)/2)
        read_ext = ['R1','R2']
    for tec_rep in range(num_reps):
        if num_reps > 1:
            prefix = "{}_tr{}".format(name,int(tec_rep)+1)
        trimming_dict[prefix] = dict()
        trimming_dict[prefix]['trimmed']  = [os.path.join(PATH_TRIMMED,name, "{}_{}.fastq.gz".format(prefix,read)) for read in read_ext]
        if lib_type == 'paired':
            trimming_dict[prefix]['raw']      = [fqfiles[tec_rep*2:(tec_rep+1)*2]]
        else:
            trimming_dict[prefix]['raw']      = [fqfiles[tec_rep]]
    return(trimming_dict)

# ---------------------------------------------------------------------------- #
# given a config dictionary sets the default value
def set_default(param, default, dictionary):
	VALUE = ''
	if param in dictionary.keys() and dictionary[param] != None:
		VALUE = dictionary[param]
	else:
		VALUE = default

	return(VALUE)

# ---------------------------------------------------------------------------- #
# given a app name calls the help and parses the parameters
def get_app_params(app_name):
    app_return = subprocess.check_output(SOFTWARE[app_name]['executable'] +' '+ SOFTWARE[app_name]['help'], shell=True)
    app_return = str(app_return)
    vals = list(set(re.findall('(\-{1,2}[a-zA-Z][\w\-]*)\W' , app_return)))
    keys = [re.sub('^-+','',i) for i in vals]
    args = dict(zip(keys, vals))
    return(args)

# ---------------------------------------------------------------------------- #
# WRITE TESTS FOR THIS FUNCTION!!!
def join_params(app, app_params, params_set):
    app_name    = os.path.basename(app)
    params_all  = get_app_params(app_name)
    names       = set(params_all.keys())
    params_diff = set(params_set) - names
    if len(params_diff) > 0:
        print(app_name + 'contains unknown parameters:' + ", ".join(list(params_diff)))
        return(0)
    params_set_keys = set(params_set.keys())
    # check whether some of the arguments are not allowed
    if 'remove' in SOFTWARE[app_name].keys():
        params_set_keys = params_set_keys- set(SOFTWARE[app_name]['remove'])
    params_set_keys = list(params_set_keys)
    params = [params_all[i] +' '+ str(params_set[i]) for i in params_set_keys]
    params = " ".join(params)
    return(params)

# ---------------------------------------------------------------------------- #
# given a sample name extracts whether the sample is broadPeak or narrowPeak
def get_macs2_suffix(name, dictionary):
    sample = dictionary[name]
    suffix = 'narrowPeak'

    if not sample == None:
        if 'macs2' in set(sample.keys()):
            if 'broad' in set(sample['macs2'].keys()):
                suffix = 'broadPeak'
    return(suffix)

# ---------------------------------------------------------------------------- #
# given a list of lists, returns a flattened version
def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


# ---------------------------------------------------------------------------- #
def RunRscript(input, output, params, logfile, script):

    # if isinstance(input, list):
    #     input = dict(zip(input, input))

    params_dump = json.dumps(dict(params.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    input_dump  = json.dumps(dict(input.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    output_dump = json.dumps(dict(output.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)

    cmd = " ".join([SOFTWARE['nice']['executable'],SOFTWARE['nice']['args'],
                    str(params.Rscript),
                    SOFTWARE['Rscript']['args'],
                    os.path.join(SCRIPT_PATH, script),
                    "--basedir", SCRIPT_PATH,
                    "--input",  "{input_dump:q}",
                    "--output", "{output_dump:q}",
                    "--params", "{params_dump:q}",
                    " 2> ", logfile,
                    " 1>&2" ])

    shell(cmd)

# ---------------------------------------------------------------------------- #
# tries to make symbolic link
def trylink(infile, outfile):
    if os.path.exists(outfile):
        os.remove(outfile)
    try:
        os.symlink(infile, outfile)
    except:
        "Symbolic link exists"
