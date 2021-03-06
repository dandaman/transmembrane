import tmhmm
from Bio import SeqIO
from pathvalidate import sanitize_filepath
import pandas as pd
import os.path
import os
import itertools
import shutil

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

localrules: all, aggregate_results, split_multifasta

rule all:
	input:
		config["result"]

rule split_multifasta:
	input:
		config["proteins"]
	output:
		dynamic(config["tmp"]+"/{display_id}/protein.fasta")
	run:
		with open(input[0], "rU") as handle:
			for record in SeqIO.parse(handle,"fasta"):
				#dirname=os.path.join(output[0],record.id)
				dirname=os.path.join(config["tmp"],record.id)
				if not os.path.exists(dirname):
					os.makedirs(dirname)
				filepath="{dir}/protein.fasta".format(dir=dirname)
				filepath=sanitize_filepath(filepath)
				with open(filepath, "w") as output_handle:
					SeqIO.write(record, output_handle, "fasta")


rule run_tmhmm:
	input:
		config["tmp"]+"/{display_id}/protein.fasta"
	output:
		config["tmp"]+"/{display_id}/tmhmm.tsv"
	benchmark:
		"benchmark/{display_id}.tmhmm.log"
	threads: 1
	log:
		"log/{display_id}.tmhmm.log"
	run:
		with open(input[0], "rU") as handle:
			for record in SeqIO.parse(handle,"fasta"):
				annotation=tmhmm.predict(str(record.seq).replace("X","A").replace("B","D").replace("Z","E").replace("-",""),config["TMHMM_model"],compute_posterior=False)
				pos = [f"{start}-{end}" for start, end, state in summarize(annotation, "M")]
				start=config["TMHMM_states"][annotation[0]]
				filepath="{tmp}/{display_id}/tmhmm.tsv".format(display_id=record.id,tmp=config["tmp"])
				filepath=sanitize_filepath(filepath)
				with open(filepath, "w") as output_handle:
					print("\t".join(list(map(lambda x: str(x),[record.id, "tmhmm","False", len(pos),start,','.join(pos)]))),file=output_handle)

rule run_memsat_svm:
	input:
		config["tmp"]+"/{display_id}/protein.fasta"
	output:
		config["tmp"]+"/{display_id}/memsat.tsv"
	benchmark:
		"benchmark/{display_id}.memsat_svm.log"
	log:	
		"log/{display_id}.memsat_svm.log"
	threads: config["MEMSAT_cores"]
	shell:
		"mkdir -p {config[tmp]}/{wildcards.display_id}/input 2>{log};"
		"mkdir -p {config[tmp]}/{wildcards.display_id}/output 2>{log};"
		"{config[MEMSAT_executable]} {input[0]} -d {config[MEMSAT_db]} -g 0 -c {config[MEMSAT_cores]} -i {config[tmp]}/{wildcards.display_id}/input -j {config[tmp]}/{wildcards.display_id}/output -e 1 >{log} 2>&1;"
		"perl parse_memsat_svm.pl {config[tmp]}/{wildcards.display_id}/output/protein.memsat_svm > {output} 2>{log}"

rule aggregate_results:
	input:
		dynamic(config["tmp"]+"/{display_id}/tmhmm.tsv"),
		dynamic(config["tmp"]+"/{display_id}/memsat.tsv")
	output:
		config["result"]
	run:
		o=list()
		for f in input:
			d=pd.read_csv(f,sep="\t",names=["protein","method","signal_peptide","membrane_domains","protein_starts","topology"])
			o.append(d)
		o=pd.concat(o)
		o['signal_peptide']=o.signal_peptide.astype('bool')
		o.to_csv(output[0],sep="\t",index=False)

onsuccess:
	dirpath=config["tmp"]
	if os.path.exists(dirpath) and os.path.isdir(dirpath):
		shutil.rmtree(dirpath)
