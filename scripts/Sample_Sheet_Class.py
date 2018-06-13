import csv
import inspect

## Experiment Class
class experiment:

    # ----------------------------------------------------------------------------- #
    def __init__(self, config = [], name="sc_rnaseq"):
        """Return an experiment object whose name is *name*"""
        self.config = config
        self.name = name
        self.SAMPLE_SHEET = []
    
    # ----------------------------------------------------------------------------- #
    def init_SAMPLE_SHEET(self, PATH_SAMPLE_SHEET):
        """Load the SAMPLE_SHEET as csv and set the *SAMPLE_SHEET* attribute"""
        self.name = PATH_SAMPLE_SHEET
        with open(PATH_SAMPLE_SHEET, 'r') as fp:
            rows =  [row for row in csv.reader(fp, delimiter=',')]
            header = rows[0]; rows = rows[1:]
            # removes empty rows from the sample_sheet
            SAMPLE_SHEET = [dict(zip(header, row)) for row in rows if len(row) > 0]
            self.SAMPLE_SHEET = SAMPLE_SHEET
            
        self.validate_sheet_init()
    
    # ----------------------------------------------------------------------------- #
    def lookup(self, column, predicate, fields=[]):
        """Function for fetching elements from the *SAMPLE_SHEET*"""
        if inspect.isfunction(predicate):
            records = [line for line in self.SAMPLE_SHEET if predicate(line[column])]
        else:
            records = [line for line in self.SAMPLE_SHEET if line[column]==predicate]
        return [record[field] for record in records for field in fields]

    # ----------------------------------------------------------------------------- #        
    def list_attr(self, attr):
        """Function for listing values for *atrr*"""
        return [line[attr] for line in self.SAMPLE_SHEET]

    # ----------------------------------------------------------------------------- #        
    def list_rows(self, column, value):
        records = [line for line in self.SAMPLE_SHEET if line[column]==value]
        return records
 
    # ----------------------------------------------------------------------------- #       
    def collapse_technical_replicates(self, PATH_TO_NEW):
        """Function for collapsing technical replicates"""
        sample_names = set(self.list_attr('sample_name'))
        to_collapse = {}
        collapsed = {}
        tmp = []
        for name in sample_names:
            barcode = set(self.lookup('sample_name', name, ['barcode']))
            reads = set(self.lookup('sample_name', name, ['reads']))
            to_collapse[name]= {'barcode': barcode, 'reads': reads}
            collapsed[name] = {'barcode': name + '_barcode.fastq.gz', 
                                'reads': name + '_reads.fastq.gz'}
            tmp.append(self.list_rows('sample_name', name)[0])
            tmp[len(tmp) - 1]['barcode'] = collapsed[name]['barcode']
            tmp[len(tmp) - 1]['reads'] = collapsed[name]['reads']
        self.SAMPLE_SHEET = tmp
        setattr(self, "to_collapse", to_collapse)
        setattr(self, "collapsed", collapsed)
        with open(PATH_TO_NEW, 'w') as new_file:
            data = self.SAMPLE_SHEET
            columns = list(data[0].keys())
            dict_writer = csv.DictWriter(new_file, columns)
            dict_writer.writeheader()
            dict_writer.writerows(data)

    # ----------------------------------------------------------------------------- #            
    def fetch_reads(self, sample_name):
        return self.to_collapse[sample_name]['reads']

    # ----------------------------------------------------------------------------- #        
    def fetch_barcodes(self, sample_name):
        return self.to_collapse[sample_name]['barcode']

    # ----------------------------------------------------------------------------- #        
    def get_fastq_files(self, s_name, slot):
        """Function to get the FASTQ files associated with a sample"""
        if slot == 'reads':
            reads = []
            for fq_file in self.fetch_reads(s_name):
                reads.append(os.path.join(PATH_FASTQ, fq_file))
            reads.sort()
            return reads
        elif slot == 'barcode':
            barcodes = []
            for fq_file in self.fetch_barcodes(s_name):
                barcodes.append(os.path.join(PATH_FASTQ, fq_file))
            barcodes.sort()
            return barcodes

    # ----------------------------------------------------------------------------- #            
    def validate_sheet_init(self):
        """Function to validate the sample sheet"""
        sample_sheet = self.SAMPLE_SHEET
        
        # Check if the required fields are found in the sample sheet
        required_fields = set(['sample_name', 'barcode', 'reads', 'method'])
        not_found = required_fields.difference(set(sample_sheet[0].keys()))
        if len(not_found) > 0:
            raise Exception("ERROR: Required field(s) {} could not be found in the sample sheet file '{}'".format(not_found, self.name))

        # Check if any of the input files are generated by unsupported single-cell rna-seq methods
        methods = set(self.config['adapter_parameters'].keys())
        for sample in sample_sheet:
            method = sample['method']
            if not method in methods:
                message = 'Sample sheet contains unknown method:' + sample + '\n'
                message = message + 'Supported methods are:' + " ".join(list(methods)) + '\n'
                sys.exit(message)
        
        samples = {}
        # Check that reads files exist
        for row in sample_sheet:
            filenames = [row['barcode'], row['reads']]
            for filename in filenames:
                fullpath = os.path.join(self.config['locations']['reads-dir'], filename)
                if not os.path.isfile(fullpath):
                    raise Exception('ERROR: missing reads file: {}'.format(fullpath))
