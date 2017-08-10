#!/usr/bin/python

import json, io
from snakemake.utils import update_config 

# config_file = open("config.json")
# config = json.load(config_file)

min_config = { 
	"params": { 
		"params_bam_methCall": {
        	"mincov": 10,
        	"minqual": 20
    	}     
	}
}


config = json.load(open("config.json"))



# if "params_bam_methCall" not in config.get("params",{})

if "params" not in config:
	update_config(config, min_config)
else:
	if "params_bam_methCall" not in config.get("params"):
		update_config(config, min_config)

	# if "params_bam_methCall" not in config["params"]:
	# 	config["params"] = {}
	# 	config["params"]["params_bam_methCall"] = {"mincov": 10,"minqual": 20}
	# 	config["params"]["params_bam_methCall2"] =  {"mincov": 10,"minqual": 20}

dumps = json.dumps(config,
					indent=4, sort_keys=True,
                    separators=(',', ': '), ensure_ascii=False)

outfile = open("new.config.json","w")
outfile.write(dumps)
outfile.close()