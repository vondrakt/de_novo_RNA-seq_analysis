# Running the reference based RNA-seq analysis

Create a conda enviroment and install necessary dependecies

`mamba env create -n de-novo_RNA-seq_analysis -f DE_analysis.yaml`

Activate the conda enviroment

`conda activate de-novo_RNA-seq_analysis`

Running the master python script

`./master_RNA_seq.py  -p /path/to/directory/ -s sample1,sample2,sample3 -r repA,repB,repC -c config_DE.txt -x GENE/TRANSCRIPT`
