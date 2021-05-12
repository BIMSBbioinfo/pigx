# Installation 

### Download/Install the pigx-sarscov2-ww

[...]

### Prepare databases

Before the pipeline can be run, 3 databases must be downloaded and their location will need to be provided in the settings file. Depending on the size of the databses this can take some time.
Be sure that the pigx-sarscov2-ww pipeline is downloaded and the tools are installed, bevore preparing the databases. The folder structure is suggested like: 
```
databases
│
│__Kraken_db
│
│__Krona_db
│
│__?
``` 

#### Kraken2 database

There are several libraries of genomes that can be used to classify the (unaligned) reads. It is up to you which one to use, but be sure that they fulfill the necessity stated by Kraken2 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases). We recommend to use the Plus-PFP library provided [here](https://benlangmead.github.io/aws-indexes/k2). 
It is also possible to have multiple Kraken2 databases, just be sure to provide the wanted one to the settings file.

First download and unpack the database in the ~databases/Kraken_db/:
```
 wget/rsync [...]
 tar -xvf 
```

Next move one folder up to databases/ and run:
```
DBNAME=$Kraken_DB	# must point towards the location/folder of the database
kraken2-build --download-taxonomy --db $DBNAME
kraken2-build --build --db $DBNAME
```

#### Krona database

Krona Tools needs a two files, which have to be installed in the ~databases/Krona_db folder:
```
FOLDER=$krona_db
~/.guix-profile/share/krona-tools/updateTaxonomy.sh $folder	#the scripts are stored a priori in that folder
~/.guix-profile/share/krona-tools/updateAccessions.sh $folder
```

Now the folder structure should look like:
```
databases
│
│__Kraken_db
│   │_hash.k2d
│   │_opts.k2d
│   │_taxo.k2d
│
│__Krona_db
│   │_all.accession2taxid.sorted
│   │_taxonomy.tab
│
│__?

```


#### settings sheet 

provide both locations to the settings sheet 
