### Example command for slurm: snakemake --cluster "sbatch -n 4 --qos=regular --constraint=cpu -A m342 -t 1:00:00" --default-resources --jobs 10 --rerun-incomplete --keep-going
### genomad: Identifies virus and plasmid, taxonomic assignment of viral genomes,annotation of proteins 
### CheckV: removes host contamination, identifies viral genome completeness, identifies closed genomes 
configfile: "config.yaml"
#print (config['samples'])

from glob import glob
import pandas as pd
import numpy as np
np.set_printoptions(threshold=sys.maxsize)


rule all:
	input:
		expand("genomad_output/{sample}/{sample}_summary/{sample}_virus.fna", sample=config["samples"]),
		expand("genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa", sample=config["samples"]),
		expand("genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv", sample=config["samples"]),
		expand("checkv_output/{sample}_completeness.tsv", sample=config["samples"]),
		expand("checkv_output/{sample}_completeness.name.tsv", sample=config["samples"],allow_missing=True),
		expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq", sample=config["samples"]),
		expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs", sample=config["samples"]),
		expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs.tsv", sample=config["samples"],allow_missing=True),
###		expand("phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv", sample=config["samples"],allow_missing=True),
		expand("plasmid_output/{sample}_plasmid_summary.name.tsv", sample=config["samples"],allow_missing=True),
		"checkv_output/all_genome_completeness.tsv",
		"phrog_output/all_virus_proteins.csv",
		"plasmid_output/all_plasmids.tsv",
		"checkv_output/all_genome_completeness.high.tsv",
		"checkv_output/all_genome_completeness.low.tsv",
		"phrog_output/all_virus_proteins.high.csv",
		"phrog_output/all_virus_proteins.low.csv",
		"plasmid_output/all_plasmids.high.tsv",
		"plasmid_output/all_plasmids.low.tsv"

rule genomad:
	input:
		"bacterial_genomes/{sample}.fna"
	output:
		virus_fna="genomad_output/{sample}/{sample}_summary/{sample}_virus.fna",
		virus_proteins_faa="genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa",
		plasmid_faa="genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
	params:
		dir=directory("genomad_output/{sample}")
	shell:
		"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {input} {params.dir} genomad_db"

rule checkv:
	input:
		fa="genomad_output/{sample}/{sample}_summary/{sample}_virus.fna"
	output:
		ind="checkv_output/{sample}_completeness.tsv"
	params:
		dir=directory("checkv_output/{sample}"),
		ind="checkv_output/{sample}/completeness.tsv"
	shell:
		"""
		checkv end_to_end {input.fa} {params.dir} -d checkv-db-v1.5 -t 16
		cp {params.ind} {output.ind}
		"""


rule checkv_edit:
	input:
		indf="checkv_output/{sample}_completeness.tsv"
	output:
		name="checkv_output/{sample}_completeness.name.tsv"
	params:
		sample="'{sample}'"
	run:
		checkv_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
		print(params)
		checkv_table['host']=params.sample
		checkv_table.to_csv(output.name, sep='\t')

rule mmseqs:
	input:
		"genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa"
	output:
		fa="mmseqs_target_seq/{sample}/{sample}_virus_proteins.faa",
		file="mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq"
	conda:
		"bin/mmseqs2.yaml"
	shell:
		"""
		cp {input} {output.fa} 
		mmseqs createdb {output.fa} {output.file}
		"""

rule phrog:
	input:
		phrog_db="phrogs_mmseqs_db/phrogs_profile_db",
		target="mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq"
	output:
		mmseqs="mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs",
		tsv="mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs.tsv"
	params:
		tmp="mmseqs_target_seq/{sample}/tmp",
		tsv="mmseqs_target_seq/{sample}_virus_proteins_mmseqs.tsv"
	conda:
		"bin/mmseqs2.yaml"
	shell:
		"""
		mmseqs search {input.phrog_db} {input.target} {output.mmseqs} {params.tmp} -s 7
		mmseqs createtsv {input.phrog_db} {input.target} {output.mmseqs} {output.tsv}
		mmseqs createtsv {input.phrog_db} {input.target} {output.mmseqs} {output.tsv} --full-header 
		cp {output.tsv} mmseqs_target_seq
		echo "file: {params.tsv}"
		"""


rule phrog_edit:
	input:
		indf=lambda wildcards: "mmseqs_target_seq/{sample}_virus_proteins_mmseqs.tsv"
	output:
		indf="phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv"
	params:
		sample="'{sample}'"
	run:
		print("input "+input.indf)
		print("param "+params.sample)
		print("output "+output.indf)
		phrog_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore',header=None)
		phrog_table.columns=['#phrog','host_seq','alnScore','seqIdentity','eVal','qStart','qEnd','qLen','tStart','tEnd','tLen']
		phrog_table[['#phrog','phrog_seq']] =phrog_table['#phrog'].str.split(' ## ',expand=True)
		df_index = pd.read_csv('phrogs_mmseqs_db/PHROG_index.csv')
		df_Bins_Index = pd.merge(phrog_table,df_index, how = 'left', on = '#phrog').fillna('NA')
		df_Bins_Index['host']=params.sample
		np.set_printoptions(threshold=sys.maxsize)
		df_Bins_Index.to_csv(output.indf)

rule plasmid_edit:
	input:
		indf="genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
	output:
		name="plasmid_output/{sample}_plasmid_summary.name.tsv"
	params:
		sample="'{sample}'"
	run:
		plasmid_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
		print(params)
		plasmid_table['host']=params.sample
		plasmid_table.to_csv(output.name, sep='\t')


rule aggregate:
	input:
		genome=glob("checkv_output/{sample}_completeness.name.tsv".format(sample="*")),
		gene=glob("phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv".format(sample="*")),
		plasmid=expand("plasmid_output/{sample}_plasmid_summary.name.tsv", sample=config["samples"],allow_missing=True)
	output:
		genome="checkv_output/all_genome_completeness.tsv",
		gene="phrog_output/all_virus_proteins.csv",
		plasmid="plasmid_output/all_plasmids.tsv"
	shell:
		"""
		cat {input.genome} > {output.genome}
		cat {input.gene} > {output.gene}
		cat {input.plasmid} > {output.plasmid}
		"""

rule filter_viral_genomes:
	input:
		all="checkv_output/all_genome_completeness.tsv"
	output:
		high="checkv_output/all_genome_completeness.high.tsv",
		low="checkv_output/all_genome_completeness.low.tsv"
	run:
		import pandas as pd
		completeness = pd.read_csv(input.all, sep='\t')
		#remove extra headers that got into main dataframe
		completeness=completeness.loc[(completeness['contig_id'] != 'contig_id') & (completeness['contig_length'] != 'contig_length')]
		#fill NAs for aai_completeness and aai_confidence, and make columns numeric
		completeness['aai_completeness']=completeness['aai_completeness'].fillna(0)
		completeness['aai_confidence']=completeness['aai_confidence'].fillna('low')
		completeness['aai_completeness']=pd.to_numeric(completeness['aai_completeness'])
		#higher than 95% aai_completeness and goes to high quality table
		completeness_high=completeness.loc[(completeness['aai_confidence'] == 'high') & (completeness['aai_completeness'] >= 95)] 
		#less than 95% aai_completeness or 95% completeness but hot high confidence goes to low quality table
		completeness_low=completeness.loc[((completeness['aai_confidence'] != 'high')&(completeness['aai_completeness'] >= 95))|(completeness['aai_completeness'] < 95)]
		#write tables
		completeness_high.to_csv(output.high,sep='\t')
		completeness_low.to_csv(output.low,sep='\t')

rule filter_viral_genes:
	input:
		all="phrog_output/all_virus_proteins.csv"
	output:
		high="phrog_output/all_virus_proteins.high.csv",
		low="phrog_output/all_virus_proteins.low.csv"
	run:
		genes = pd.read_csv(input.all)
		#remove extra headers that got into main dataframe
		genes=genes.loc[(genes['#phrog'] != '#phrog') & (genes['host_seq'] != 'host_seq')]
		#e-value cutoff of 10E-6 or less (depends on the length of the query sequence and size of database. probability of occurance by chance alone 
		genes_high=genes.loc[pd.to_numeric(genes['eVal'])<=1e-6]
		genes_low=genes.loc[pd.to_numeric(genes['eVal'])>1e-6]
		#write tables
		genes_high.to_csv(output.high,sep='\t')
		genes_low.to_csv(output.low,sep='\t')


checkpoint filter_plasmid_genes:
	input:
		all="plasmid_output/all_plasmids.tsv"
	output:
		high="plasmid_output/all_plasmids.high.tsv",
		low="plasmid_output/all_plasmids.low.tsv"
	run:
		plasmid = pd.read_csv(input.all,sep='\t')
		plasmid=plasmid.loc[(plasmid['seq_name'] != 'seq_name') & (plasmid['length'] != 'length')]
		plasmid_high=plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])>=0.95]
		plasmid_low=plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])<0.95]
		plasmid_high.to_csv(output.high,sep='\t')
		plasmid_low.to_csv(output.low,sep='\t')
