#!/usr/bin/env python3
from optparse import OptionParser
import os
import re

# THIS IS A MASTER PYTHON SCRIPT FOR RUNNING A DE-NOVO RNA-SEQ ANALYSIS

# defining the arguments that need to be passed to the scripts
arguments = OptionParser()

arguments.add_option('-p', '--path', dest='path', help='path to wroking directory')
arguments.add_option('-s', '--samples', dest='samples', help='sample IDs, separated by comma')
arguments.add_option('-r', '--replicates', dest='replicates', help='replicate IDs, separated by comma')
arguments.add_option('-c', '--config-file', dest='config_DE', help='configuration file for the DE analysis')
arguments.add_option('-x', '--count-option', dest='count_option', help='the DE analysis will be performed using either the gene or transcript count matrix,'
                                                                       ' the input is either GENE or TRANSCRIPT')

(options, args) = arguments.parse_args()
if options.path is None or options.samples is None or options.replicates is None\
        or options.config_DE is None or options.count_option is None:
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Due to incompatibility with rnaspades, this step is skipped

os.chdir(options.path)
print(os.getcwd())

# Processing names
samples = options.samples.split(',')
replicates = options.replicates.split(',')

for replicate in replicates:
    name1 = options.path+samples[0]+replicate
    name2 = options.path+samples[1]+replicate
    print(name1, name2)
    # Trimming the reads
    trimming_command = 'trimmomatic PE -threads 30 -phred33 %s.fastq.gz' \
                       ' %s.fastq.gz %s_1P.fastq.gz %s_1U.fastq.gz %s_2P.fastq.gz %s_2U.fastq.gz ' \
                       'ILLUMINACLIP:illumina_multiplex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' % (name1, name2, name1, name1, name2, name2)
    print(trimming_command)
    os.system(trimming_command)

    # Fastqc
    fastqc_command = 'fastqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    multiqc_command = 'multiqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    print(fastqc_command)
    print(multiqc_command)
    os.system(fastqc_command)
    os.system(multiqc_command)

# Transcriptome de novo assembly
# print(samples)
# print(replicates)
names_sample1 = []
names_sample2 = []
for replicate in replicates:
    names_sample1.append(samples[0]+replicate+'_1P.fastq.gz')
    names_sample2.append(samples[1] + replicate + '_2P.fastq.gz')
names_sample1 = ','.join(names_sample1)
names_sample2 = ','.join(names_sample2)

# RUNNING TRINITY FOR ALL SAMPLES AND REPLICATES TOGETHER
print(names_sample1, names_sample2)
# added the no_salmon parameter
Trinity_command = 'Trinity --seqType fq --single %s,%s --CPU 8 --max_memory 20G --min_contig_length 150 --no_salmon' % (names_sample1, names_sample2)
print(Trinity_command)
os.system(Trinity_command)

# Structural annotation of the reference transcriptome created by Trinity
long_orf = 'TransDecoder.LongOrfs -t trinity_out_dir.Trinity.fasta'
predict_coding_seq = 'TransDecoder.Predict -t trinity_out_dir.Trinity.fasta'
convert_from_gff_to_gtf = 'agat_convert_sp_gff2gtf.pl -gff trinity_out_dir.Trinity.fasta.transdecoder.gff3 -o trinity_out_dir.Trinity.fasta.transdecoder_agat.gtf'

os.system(long_orf)
os.system(predict_coding_seq)
os.system(convert_from_gff_to_gtf)

# The STAR alignment and read counting is done for each replicate of each sample separately
print(samples)

# The reference will be done just once as to save on time and memory
STAR_reference = 'STAR --runMode genomeGenerate --runThreadN 30 --genomeDir trinity_out_dir' \
                 '  --genomeFastaFiles trinity_out_dir.Trinity.fasta --genomeSAindexNbases 10' \
                 ' --sjdbGTFfile trinity_out_dir.Trinity.fasta.transdecoder_agat.gtf'
os.system(STAR_reference)

for replicate in replicates:
    print(replicate)

    # STAR alignment
    STAR_command_reads_sample1 = 'STAR --genomeDir trinity_out_dir --runThreadN 10 --readFilesIn %s_1P.fastq.gz ' \
                                 '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s.sorted.bam -' \
                                 '-outSAMtype BAM SortedByCoordinate' % (samples[0] + replicate, samples[0] + replicate)

    STAR_command_reads_sample2 = 'STAR --genomeDir trinity_out_dir --runThreadN 10 --readFilesIn %s_2P.fastq.gz ' \
                                 '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s.sorted.bam -' \
                                 '-outSAMtype BAM SortedByCoordinate' % (samples[1] + replicate, samples[1] + replicate)

    print(STAR_command_reads_sample1)
    print(STAR_command_reads_sample2)
    os.system(STAR_command_reads_sample1)
    os.system(STAR_command_reads_sample2)

    # count reads and expression matrix based star BAM
    # Create a stringtie samples file
    stringtie_command_sample1 = 'stringtie %s.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G trinity_out_dir.Trinity.fasta.transdecoder_agat.gtf -e -o %s.gtf -A %s.gene_abundances.tsv' % \
                                (samples[0] + replicate, samples[0] + replicate, samples[0] + replicate)

    stringtie_command_sample2 = 'stringtie %s.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G trinity_out_dir.Trinity.fasta.transdecoder_agat.gtf -e -o %s.gtf -A %s.gene_abundances.tsv' % \
                                (samples[1] + replicate, samples[1] + replicate, samples[1] + replicate)
    print(stringtie_command_sample1)
    print(stringtie_command_sample2)
    os.system(stringtie_command_sample1)
    os.system(stringtie_command_sample2)

# Create the samples.txt file
samples_txt = open('prepDE_samples.txt', 'w')

for sample in samples:
    for replicate in replicates:
        list = [sample+replicate, options.path+sample+replicate+'.gtf']
        list = '\t'.join(list)
        samples_txt.write(list+'\n')

samples_txt.close()
# Generating count matrices
command = 'prepDE.py -i prepDE_samples.txt'
print(command)
os.system(command)

# DE analysis
# The configuration file for the DE analysis must be in format
# condition1    condition1_R1
# condition1    condition1_R2
# condition1    condition1_R3
# condition2    condition2_R1
# condition2    condition2_R2
# condition2    condition2_R3

if options.count_option == 'GENE':
    # first convert the .csv file to a .tsv file
    out_tsv = open('./gene_count_matrix.tsv', 'w')
    with open('./gene_count_matrix.csv') as c:
        for line in c:
            line = re.sub(',', '\t', line)
            out_tsv.write(line)
    out_tsv.close()

    DE_analysis_command_edgeR = 'run_DE_analysis.pl --matrix gene_count_matrix.tsv --samples_file %s%s ' \
                        '--reference_sample condition1 --method edgeR --output edgeR_genes' % (options.path, options.config_DE)
    print(DE_analysis_command_edgeR)
    os.system(DE_analysis_command_edgeR)

    DE_analysis_command_DESeq2 = 'run_DE_analysis.pl --matrix gene_count_matrix.tsv --samples_file %s%s ' \
                          '--reference_sample condition1 --method DESeq2 --output DESeq2_genes' % (options.path, options.config_DE)
    print(DE_analysis_command_DESeq2)
    os.system(DE_analysis_command_DESeq2)

    # Extracting differentially expressed transcripts and generating heatmaps
    # Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed
    # at a significance of <= 0.05 in any of the pairwise sample comparisons:

    # The working directory needs to be changed
    new_dir = '%sedgeR_genes' % options.path
    print(new_dir)
    os.chdir(new_dir)
    analyze_DE = 'analyze_diff_expr.pl --matrix ../gene_count_matrix.tsv --samples %s%s -P 0.05 -C 2' % (options.path,options.config_DE)
    print(analyze_DE)
    os.system(analyze_DE)

    new_dir = '%sDESeq2_genes' % options.path
    print(new_dir)
    os.chdir(new_dir)
    analyze_DE = 'analyze_diff_expr.pl --matrix ../gene_count_matrix.tsv --samples %s%s -P 0.05 -C 2' % (options.path,options.config_DE)
    print(analyze_DE)
    os.system(analyze_DE)

else:
    out_tsv = open('./transcript_count_matrix.tsv', 'w')
    with open('./transcript_count_matrix.csv') as c:
        for line in c:
            line = re.sub(',', '\t', line)
            out_tsv.write(line)
    out_tsv.close()

    DE_analysis_command_edgeR = 'run_DE_analysis.pl --matrix transcript_count_matrix.tsv --samples_file %s%s ' \
                                '--reference_sample condition1 --method edgeR --output edgeR_genes' % (options.path, options.config_DE)
    print(DE_analysis_command_edgeR)
    os.system(DE_analysis_command_edgeR)

    DE_analysis_command_DESeq2 = 'run_DE_analysis.pl --matrix transcript_count_matrix.tsv --samples_file %s%s ' \
                                 '--reference_sample condition1 --method DESeq2 --output DESeq2_genes' % (options.path, options.config_DE)
    print(DE_analysis_command_DESeq2)
    os.system(DE_analysis_command_DESeq2)

    # Extracting differentially expressed transcripts and generating heatmaps
    # Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed
    # at a significance of <= 0.05 in any of the pairwise sample comparisons:

    # The working directory needs to be changed
    new_dir = '%sedgeR_genes' % options.path
    print(new_dir)
    os.chdir(new_dir)
    analyze_DE = 'analyze_diff_expr.pl --matrix ../transcript_count_matrix.tsv --samples %s%s -P 0.05 -C 2' % (options.path, options.config_DE)
    print(analyze_DE)
    os.system(analyze_DE)

    new_dir = '%sDESeq2_genes' % options.path
    print(new_dir)
    os.chdir(new_dir)
    analyze_DE = 'analyze_diff_expr.pl --matrix ../transcript_count_matrix.tsv --samples %s%s -P 0.05 -C 2' % (options.path, options.config_DE)
    print(analyze_DE)
    os.system(analyze_DE)