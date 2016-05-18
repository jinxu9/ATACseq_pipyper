#!/usr/bin/env python
"""
ATACseq  pipeline
"""
__author__="Jin Xu"
__email__="xujin937@gmail.com"


from argparse import ArgumentParser
import os
import os.path
import sys
from subprocess import call
import subprocess
import re
import pypiper
import yaml 
from datetime import datetime




# Argument Parsing
# #######################################################################################
parser = ArgumentParser(description='Pypiper arguments.')
parser = pypiper.add_pypiper_args(parser, all_args=True)

#Add any pipeline-specific arguments
#parser.add_argument('-e', '--ercc', default="ERCC92",
#parser.add_argument('-em', '--ercc-mix',
#parser.add_argument('-f', dest='filter', action='store_false', default=True)
# Core-seq as optional parameter
#parser.add_argument('-cs', '--core-seq', default=False, dest='coreseq', action='store_true', help='CORE-seq Mode')
args = parser.parse_args()

if args.single_or_paired == "paired":
	args.paired_end = True
else:
	args.paired_end = False

###
class Container:
	pass

paths = Container()
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))
pipelines_config_file = os.path.join(paths.scripts_dir, "Configure_ATACseq.yaml")
config = yaml.load(open(pipelines_config_file, 'r'))

# Resources
paths.bowtie_ref= config["resources"]["ref"]
paths.ref_size=config["resources"]["ref_size"]
paths.refGene_TSS=config["resources"]["refGene_TSS"]
paths.blacklist=config["resources"]["blacklist"]
paths.adaptor=config["resources"]["adaptor"]
# Tools
paths.java=config["tools"]["java"]
paths.bowtie2=config["tools"]["bowtie2"]
paths.samtools =config["tools"]["samtools"]
paths.bedtools=config["tools"]["bedtools"]
paths.rmdup=config["tools"]["MarkDuplicates"]
paths.norm_bedgraph=config["tools"]["norm_bedGraph"]
paths.MakeVplot=config["tools"]["pyMakeVplot"]
paths.fragment_pl=config["tools"]["fragment_length_dist_pl"]
paths.fragment_R=config["tools"]["fragment_length_dist_R"]
paths.bam2bed_shift=config["tools"]["bam2bed_shift"]
#paths.adapterTrim=config["tools"]["adapterTrim"]
paths.trimmo=config["tools"]["trimmo"]
# Output
paths.pipeline_outfolder = os.path.join(args.output_parent, args.sample_name + "/")
# Initialize
mypiper = pypiper.PipelineManager(name="ATACseq", outfolder=paths.pipeline_outfolder, args=args)
myngstk = pypiper.NGSTk(pm=mypiper)
#args.paired_end=T

################################################################################
#local_input_file = myngstk.create_local_input(args.output_parent, args.input, args.sample_name)
print("Local input file: " + args.input[0]) 
print("Local input file: " + args.input2[0]) 
#ATACseq pipeline

# Adaptor trimming
mypiper.timestamp("### Adapter trimming: ")

output = paths.pipeline_outfolder + args.sample_name
trimmed_fastq= output +"_R1_trimmed.fq"
trimmed_fastq_R2= output +"_R2_trimmed.fq"
cmd = paths.java +" -Xmx" + str(mypiper.mem) +" -jar " + paths.trimmo + " PE " + " -threads " + str(mypiper.cores) + " "
cmd += args.input[0] + " "
cmd += args.input2[0] + " " 
cmd += trimmed_fastq + " "
cmd += output + "_R1_unpaired.fq "
cmd += trimmed_fastq_R2 + " "
cmd += output + "_R2_unpaired.fq "
#cmd +=  "ILLUMINACLIP:"+ paths.adaptor + ":2:30:10:LEADING:3TRAILING:3SLIDINGWINDOW:4:15MINLEN:36" 
cmd +=  "ILLUMINACLIP:"+ paths.adaptor + ":2:30:10" 

def check_trim():
        n_trim = float(myngstk.count_reads(trimmed_fastq, args.paired_end))
        rr = float(mypiper.get_stat("Raw_reads"))
        mypiper.report_result("Trimmed_reads", n_trim)

        mypiper.report_result("Trim_loss_rate", round((rr - n_trim) * 100 / rr, 2))

mypiper.run(cmd, trimmed_fastq, follow = check_trim)
# End of Adaptor trimming 




mypiper.stop_pipeline()
