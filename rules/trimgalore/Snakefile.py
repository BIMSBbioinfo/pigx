


PATHIN = config['PATHIN']
DIR_trimmed = config['DIR_trimmed']


#
# #
# # ==========================================================================================
# # trim the reads
# 
def get_trimgalore_input(wc):

   print('---------------wc------get_trimgalore_input-----1111')
   print(list(wc))
   samps = config['SAMPLES'][wc.sample]['fastq']

   if type(samps) is str:
        samps = [samps]

   if(len(samps)==2):
     return( [os.path.join(PATHIN, samps[0]), os.path.join(PATHIN, samps[1])] )
   else:
     return( [os.path.join(PATHIN, samps[0])] )


rule trimgalore_pe:
    input:
        #infile = get_trimgalore_input
        #infile = lambda wc: [os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0]), os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][1])] if len(config['SAMPLES'][wc.sample]['fastq_name'])==2 else [os.path.join(PATHIN, config['SAMPLES'][wc.sample]['fastq_name'][0])]
        infile = lambda wc: [os.path.join(PATHIN, config['SAMPLES'][wc.name]['fastq_name'][0]+".fq.gz"), os.path.join(PATHIN, config['SAMPLES'][wc.name]['fastq_name'][1]+".fq.gz")]
    output:
        output1 = DIR_trimmed+"{name}_val_1.fq.gz",
        output2 = DIR_trimmed+"{name}_val_2.fq.gz",
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir " + DIR_trimmed,
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
    log:
        log=DIR_trimmed+"{name}.trimgalore.log"
    message:
        " ---------  Trimming raw paired-end read data using {TRIMGALORE} -------  "
    run:

        print("---------------------print trim galore rule")
        if(len(input.infile)==0):
          print("AAAAAAAAAAA")
          print(output)
          print(intput)

        cmd = " ".join(
        [TRIMGALORE,
         " --paired",
         input.infile[0], input.infile[1],
         params.extra, params.outdir, params.gz, params.cutadapt,
         "-o "+ DIR_trimmed,
         " 2> ", log.log
         ])

        # Here a hack to generate files named like samples and not like units, even though
        # trimgalore will produce files named using units
        cmd = cmd + "; touch "+output.output1 + "; touch "+output.output2
        shell(cmd)


rule trimgalore_se:
    input:
        infile=PATHIN+'{sample}'+'.fq.gz'
    output:
        output = DIR_trimmed+"{sample}_trimmed.fq.gz"
    params:
        extra          = config.get("trim_galore_args", ""),
        outdir         = "--output_dir " + DIR_trimmed,
        gz             = "--gzip",
        cutadapt       = "--path_to_cutadapt "+ CUTADAPT,
    log:
        log=DIR_trimmed+"{sample}.trimgalore.log"
    message:
        " ---------  Trimming raw single-end read data using {TRIMGALORE} -------  "
    run:
        args = " ".join(
        [params.extra,
         params.outdir,
         params.gz,
         params.cutadapt,
         "-o "+ DIR_trimmed
         ])
        cmd=TRIMGALORE+args + " "+input.infile+" 2> "+ log.log
        cmd = cmd + "; touch "+output.output
        shell(cmd)
