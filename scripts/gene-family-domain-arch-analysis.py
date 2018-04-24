#!/usr/bin/python

#Author: Akshay Yadav
#Version: 1.0

from __future__ import division
import re
import os
import subprocess
import operator
import sys
import argparse

parser = argparse.ArgumentParser(description="Tool for calculating domain composition and domain compostion scores for given set of gene families")
parser._optionals.title="Arguments"

parser.add_argument('--fasta_dir', help="Location of directory containing family fasta files", required=True, dest="family_fasta_file_dir")
parser.add_argument('--output_dir', help="Location of the output directory", required=True, dest="output_directory")
parser.add_argument('--name', help="Name of the dataset used to create output files", required="True", dest="dataset_name")
args = parser.parse_args()


###############################################################################################################

def execute_pfamscan_on_family_fasta(family_fasta_file, family_fasta_dirName, pfamscanout_dirName):
	run_pfamscan=subprocess.Popen(["./Pfam/pfam_scan.pl", "-outfile", pfamscanout_dirName+"/"+family_fasta_file+".pfamdomtblout", "-e_seq","1e-5", "-fasta",family_fasta_dirName+"/"+family_fasta_file , "-dir","Pfam/database/"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	run_pfamscan = run_pfamscan.communicate()


def read_seqid_from_fam_fasta(family_fasta_fileName):
	family_fasta_file=open(family_fasta_fileName,"r")
	family_seqid_dict={}
	for line in family_fasta_file:
		line=line.rstrip()
		if(re.match(r'^\>',line)):
			family_seqid_dict[line[1:]]=1
	
	family_fasta_file.close()
	return(family_seqid_dict)		

def get_domain_order_from_pfamscanout_file(pfamscanout_fileName, pfamscanout_dirName, family_seqid_dict, family_domain_order_fileName):
	pfamscanout_file=open(pfamscanout_dirName+"/"+pfamscanout_fileName,"r")
	seqid_startcoord_domain_dict={}
	for line in pfamscanout_file:
		line=line.strip()
		if(re.match(r'^\#',line) or re.match(r'^$',line)):
			continue;
		linearr=re.split(r'\s+',line)
		seqid=linearr[0]
		start_coord=int(linearr[3])
		domain_name=linearr[6]
		if(seqid_startcoord_domain_dict.has_key(seqid)):
			seqid_startcoord_domain_dict[seqid][start_coord]=domain_name
		else:
			seqid_startcoord_domain_dict[seqid]={}
			seqid_startcoord_domain_dict[seqid][start_coord]=domain_name
	
	pfamscanout_file.close()
	print_domain_order_for_family_sequences(seqid_startcoord_domain_dict, family_seqid_dict, family_domain_order_fileName)

def print_domain_order_for_family_sequences(seqid_startcoord_domain_dict, family_seqid_dict, family_domain_order_fileName):
	family_domain_order_file=open(family_domain_order_fileName,"w")
	for seqid in family_seqid_dict:
		if not(seqid_startcoord_domain_dict.has_key(seqid)):
			family_domain_order_file.write(seqid+" "+"**NULL**"+"\n")
			continue
		startcoord_domain_dict=seqid_startcoord_domain_dict[seqid]
		startcoord_domain_dict_sorted=sorted(startcoord_domain_dict.items(), key=operator.itemgetter(0))
		family_domain_order_file.write(seqid+" ")
		for entry in startcoord_domain_dict_sorted:
			family_domain_order_file.write(entry[1]+" ")
		
		family_domain_order_file.write("\n")
	family_domain_order_file.close()	

def write_family_domain_order_files(family_fasta_dirName, domain_order_dirName, pfamscanout_dirName):
	
	for family_fasta_file in os.listdir(family_fasta_dirName):
		print 'processing family.....{0}'.format(family_fasta_file)
		family_fasta_fileName=family_fasta_dirName+"/"+family_fasta_file

		pfamscanout_fileName=family_fasta_file+".pfamdomtblout"
		family_domain_order_fileName=domain_order_dirName+family_fasta_file+".domain_order"


		family_seqid_dict=read_seqid_from_fam_fasta(family_fasta_fileName)
		if(len(family_seqid_dict)==0):
			continue
		execute_pfamscan_on_family_fasta(family_fasta_file, family_fasta_dirName, pfamscanout_dirName)

		get_domain_order_from_pfamscanout_file(pfamscanout_fileName, pfamscanout_dirName, family_seqid_dict, family_domain_order_fileName)

#####################################################################################################################################################
def process_domain_order_files(domain_order_dirName, dataset_name, output_dirName):
	domain_composition_outfile=open(output_dirName+"/"+dataset_name+".family_domain_compositions","w")
	domain_jaccard_score_outfile=open(output_dirName+"/"+dataset_name+".family_domain_jaccard_scores","w")
	print '\n\nreading domain order files\n\n'
	for order_file in os.listdir(domain_order_dirName):
		domain_order_fileName=domain_order_dirName+"/"+order_file
		family_id=re.split(r'\.',order_file)[0]
		print 'processing family.....{0}'.format(family_id)
		seqid_domainarr_dict=read_domain_order_files(domain_order_fileName)
		
		calculate_family_domain_composition(seqid_domainarr_dict, domain_composition_outfile, family_id)
		calculate_family_domain_jaccard_score(seqid_domainarr_dict, domain_jaccard_score_outfile, family_id)
		
	domain_composition_outfile.close()
def read_domain_order_files(domain_order_fileName):
	domain_order_file =  open(domain_order_fileName, "r")
	seqid_domainarr_dict={}
	for line in domain_order_file:
		line=line.rstrip()
		linearr=re.split(r'\s+',line)
		seqid_domainarr_dict[linearr.pop(0)]=linearr
	
	domain_order_file.close()
	return (seqid_domainarr_dict)

####-- functions for calculating family domain compositions --####
def calculate_family_domain_composition(seqid_domainarr_dict, domain_composition_outfile, family_id):
	domain_seqcount_dict={}
	for seqid in seqid_domainarr_dict:
		domain_dict=get_domain_dict_for_seq(seqid_domainarr_dict[seqid])
		update_domain_seq_count(domain_seqcount_dict, domain_dict)
	print_family_domain_composition(domain_seqcount_dict, len(seqid_domainarr_dict), domain_composition_outfile, family_id)

def get_domain_dict_for_seq(domain_arr):
	domain_dict={}
	for domain in domain_arr:
		domain_dict[domain]=1

	return(domain_dict)

def update_domain_seq_count(domain_seqcount_dict, domain_dict):
	for domain in domain_dict:
		if(domain_seqcount_dict.has_key(domain)):
			domain_seqcount_dict[domain]+=1
		else:
			domain_seqcount_dict[domain]=1
	

def print_family_domain_composition(domain_seqcount_dict, famsize, domain_composition_outfile, family_id):
	domain_seqcount_dict_sorted=sorted(domain_seqcount_dict.items(), key=operator.itemgetter(1))

	domain_composition_outfile.write(family_id+" "+str(famsize)+" ")
	for entry in domain_seqcount_dict_sorted:
		domain_id = entry[0]
		seq_percent = (entry[1]/famsize)*100
		domain_composition_outfile.write(domain_id+"-"+str(seq_percent)+" ")
	
	domain_composition_outfile.write("\n")

####-- functions for calculating family domain jaccard scores --####
def calculate_family_domain_jaccard_score(seqid_domainarr_dict, domain_jaccard_score_outfile, family_id):
	seqid_arr = seqid_domainarr_dict.keys()
	famsize=len(seqid_arr)
	if(famsize<2):
		avg_jaccard_score=1.0
	else:
		avg_jaccard_score=0
		counter=0
		for i in range(0,famsize):
			for j in range(i+1,famsize):
				jaccard_score=get_jaccard_score_for_two_domain_arr(seqid_domainarr_dict[seqid_arr[i]], seqid_domainarr_dict[seqid_arr[j]])
				avg_jaccard_score=avg_jaccard_score+jaccard_score
				counter+=1
		avg_jaccard_score=avg_jaccard_score/counter
	print_avg_jaccard_score_for_family(family_id, avg_jaccard_score, domain_jaccard_score_outfile, famsize)
	

def get_jaccard_score_for_two_domain_arr(domain_arr1, domain_arr2):
	no_of_common_domains=len(set.intersection(*[set(domain_arr1), set(domain_arr2)]))
	total_no_of_domains=len(set.union(*[set(domain_arr1), set(domain_arr2)]))
	jaccard_score=no_of_common_domains/float(total_no_of_domains)
	return (jaccard_score)

def print_avg_jaccard_score_for_family(family_id, avg_jaccard_score, domain_jaccard_score_outfile, famsize):
	domain_jaccard_score_outfile.write(family_id+" "+str(famsize)+" "+str(avg_jaccard_score)+"\n")

def execute_workflow(family_fasta_dirName, dataset_name, output_dirName):
	domain_order_dirName=output_dirName+"/"+"domain_order_results/"
	if not os.path.exists(domain_order_dirName):
    		os.makedirs(domain_order_dirName)

	pfamscanout_dirName=output_dirName+"/"+"pfamscan_results/"
	if not os.path.exists(pfamscanout_dirName):
    		os.makedirs(pfamscanout_dirName)
	
	write_family_domain_order_files(family_fasta_dirName, domain_order_dirName, pfamscanout_dirName)
	process_domain_order_files(domain_order_dirName, dataset_name, output_dirName)

############################################################################################################################################################

family_fasta_dirName=args.family_fasta_file_dir
dataset_name=args.dataset_name
output_dirName=args.output_directory
execute_workflow(family_fasta_dirName, dataset_name, output_dirName)


