import tmhmm
from Bio import SeqIO
from pathvalidate import sanitize_filepath
import pandas as pd
import os.path
import os
import itertools
import shutil

configfile: 
	"config.json"

def summarize(path,reqstate):
	"""
	Summarize a path as a list of (start, end, state) triples. Adapted from https://github.com/dansondergaard/tmhmm.py/blob/master/tmhmm/cli.py
	"""
	for state, group in itertools.groupby(enumerate(path), key=lambda x: x[1]):
		group = list(group)
		start = min(group, key=lambda x: x[0])[0]
		end = max(group, key=lambda x: x[0])[0]
		if not reqstate or state==reqstate:
			yield start, end, state

rule all:
	input:
		config["TMHMM_result"]

checkpoint split_multifasta:
	input:
		config["proteins"]
	output:
		temp(directory("tmp/"))
	run:
		with open(input[0], "rU") as handle:
			for record in SeqIO.parse(handle,"fasta"):
				dirname=os.path.join(output[0],record.id)
				if not os.path.exists(dirname):
					os.makedirs(dirname)
				filepath="{dir}/protein.fasta".format(dir=dirname)
				filepath=sanitize_filepath(filepath)
				with open(filepath, "w") as output_handle:
					SeqIO.write(record, output_handle, "fasta")

rule run_tmhmm:
	input:
		"tmp/{display_id}/protein.fasta"
	output:
		temp("tmp/{display_id}/tmhmm.tsv")
	run:
		with open(input[0], "rU") as handle:
			for record in SeqIO.parse(handle,"fasta"):
				annotation=tmhmm.predict(record.seq,config["TMHMM_model"],compute_posterior=False)
				pos = [f"{start}-{end}" for start, end, state in summarize(annotation, "M")]
				start=config["TMHMM_states"][annotation[0]]
				filepath="tmp/{display_id}/tmhmm.tsv".format(display_id=record.id)
				filepath=sanitize_filepath(filepath)
				with open(filepath, "w") as output_handle:
					print("\t".join(list(map(lambda x: str(x),[record.id, len(pos),start,','.join(pos)]))),file=output_handle)

def get_tmhmm_results(wildcards=None):
	chkpoutdir = checkpoints.split_multifasta.get().output[0]
	return(expand("tmp/{display_id}/tmhmm.tsv",display_id=glob_wildcards(os.path.join(chkpoutdir,"{display_id}/protein.fasta")).display_id))

rule aggregate_tmhmm:
	input:
		get_tmhmm_results
	output:
		config["TMHMM_result"]
	run:
		o=list()
		i=get_tmhmm_results()
#it's a mystery why {input} is tmp/ and does not comprise the return value here...
		for f in i:
			d=pd.read_csv(f,sep="\t")
			o.append(d)
		o=pd.concat(o)
		o.to_csv(output[0],sep="\t",index=False)

onsuccess:
	if  os.path.isdir('tmp'):
		shutil.rmtree('tmp/')	
