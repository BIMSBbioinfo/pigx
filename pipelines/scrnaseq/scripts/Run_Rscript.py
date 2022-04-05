# ---------------------------------------------------------------------------- #
# the rule is called in the following way:
"""
run:
    RunRscript(input, output, params, 'SCRIPT.R')
"""

def RunRscript(input, output, params, path_script, script):

    # if isinstance(input, list):
    #     input = dict(zip(input, input))

    print(input.items())
    params_dump = json.dumps(dict(params.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    input_dump  = json.dumps(dict(input.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)
    output_dump = json.dumps(dict(output.items()),sort_keys=True,
                       separators=(",",":"), ensure_ascii=True)

    cmd = " ".join(['nice -19',str(params.Rscript),
                    config['tools']['Rscript']['args'],
                    os.path.join(path_script, script),
                    "--basedir", path_script,
                    "--input",  "{input_dump:q}",
                    "--output", "{output_dump:q}",
                    "--params", "{params_dump:q}"])

    shell(cmd)
