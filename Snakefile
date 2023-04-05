### Example command for slurm: snakemake --cluster "sbatch -n 4 --qos=regular --constraint=cpu -A m342 -t 1:00:00" --default-resources --jobs 10 --rerun-incomplete --keep-going
### genomad: Identifies virus and plasmid, taxonomic assignment of viral genomes,annotation of proteins 
### CheckV: removes host contamination, identifies viral genome completeness, identifies closed genomes 
from glob import glob
configfile: "config.yaml"
print (config['samples'])

rule all:
	input:
#		expand("genomad_output/{sample}/{sample}_summary/{sample}_virus.fna", sample=config["samples"]),
#		expand("genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa", sample=config["samples"]),
#		expand("checkv_output/{sample}_completeness.tsv", sample=config["samples"]),
#		expand("checkv_output/{sample}_completeness.name.tsv", sample=config["samples"],allow_missing=True),
#               expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq", sample=config["samples"])
#		expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs", sample=config["samples"]),
#		expand("mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs.tsv", sample=config["samples"],allow_missing=True),
#		expand("phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv", sample=config["samples"],allow_missing=True),
		expand("plasmid_output/{sample}_plasmid_summary.name.tsv", sample=config["samples"],allow_missing=True),
		"checkv_output/all_genome_completeness.tsv",
		"phrog_output/all_virus_proteins.csv",
		"plasmid_output/all_plasmids.tsv"

#rule genomad:
#	input:
#		"bacterial_genomes/{sample}.fna"
#	output:
#		virus_fna="genomad_output/{sample}/{sample}_summary/{sample}_virus.fna",
#		virus_proteins_faa="genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa"
#	params:
#		dir=directory("genomad_output/{sample}")
#	shell:
#		"genomad end-to-end --min-score 0.7 --cleanup --splits 8 {input} {params.dir} genomad_db"

#rule checkv:
#	input:
#		fa="genomad_output/{sample}/{sample}_summary/{sample}_virus.fna"
#	output:
#		ind="checkv_output/{sample}_completeness.tsv"
#	params:
#		dir=directory("checkv_output/{sample}"),
#		ind="checkv_output/{sample}/completeness.tsv"
#	shell:
#		"""
#		checkv end_to_end {input.fa} {params.dir} -d checkv-db-v1.5 -t 16
#		cp {params.ind} {output.ind}
#		"""


#rule checkv_edit:
#	input:
#		indf="checkv_output/{sample}_completeness.tsv"
#	output:
#		name="checkv_output/{sample}_completeness.name.tsv"
#	params:
#		sample="'{sample}'"
#	run:
#		import pandas as pd
#		checkv_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
#		print(params)
#		checkv_table['host']=params.sample
#		checkv_table.to_csv(output.name, sep='\t')

#rule mmseqs:
#	input:
#		"genomad_output/{sample}/{sample}_summary/{sample}_virus_proteins.faa"
#	output:
#		fa="mmseqs_target_seq/{sample}/{sample}_virus_proteins.faa",
#		file="mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq"
#	shell:
#		"""
#		cp {input} {output.fa} 
#		mmseqs createdb {output.fa} {output.file}
#		"""

#rule phrog:
#	input:
#		phrog_db="phrogs_mmseqs_db/phrogs_profile_db",
#		target="mmseqs_target_seq/{sample}/{sample}_virus_proteins.target_seq"
#	output:
#		mmseqs="mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs",
#		tsv="mmseqs_target_seq/{sample}/{sample}_virus_proteins_mmseqs.tsv"
#	params:
#		tmp="mmseqs_target_seq/{sample}/tmp"
#	shell:
#		"""
#		mmseqs search {input.phrog_db} {input.target} {output.mmseqs} {params.tmp} -s 7
#		mmseqs createtsv {input.phrog_db} {input.target} {output.mmseqs} {output.tsv}
#		mmseqs createtsv {input.phrog_db} {input.target} {output.mmseqs} {output.tsv} --full-header 
#		cp {output.tsv} mmseqs_target_seq
#		"""


#def find_phrog_files(wildcards):
#        path = "mmseqs_target_seq/{sample}_virus_proteins_mmseqs.tsv"
#        files = glob(path.format(sample="*"))
#        print(files)
#        return files

#rule phrog_edit:
#	input:
#		indf=find_phrog_files
#	output:
#		indf="phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv"
#	params:
#		sample="'{sample}'"
#	run:
#		import pandas as pd
#		for x in input.indf:
#			phrog_table=pd.read_csv(x,sep = '\t', encoding_errors='ignore',header=None)
#			phrog_table.columns=['#phrog','host_seq','alnScore','seqIdentity','eVal','qStart','qEnd','qLen','tStart','tEnd','tLen']
#			phrog_table[['#phrog','phrog_seq']] =phrog_table['#phrog'].str.split(' ## ',expand=True)
#			df_index = pd.read_csv('phrogs_mmseqs_db/PHROG_index.csv')
#			df_Bins_Index = phrog_table.merge(df_index, how = 'inner', on = '#phrog').fillna('NA')
#			df_Bins_Index['host']=params.sample
#			df_Bins_Index.to_csv(output.indf)



rule plasmid_edit:
	input:
		indf="genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
	output:
		name="plasmid_output/{sample}_plasmid_summary.name.tsv"
	params:
		sample="'{sample}'"
	run:
		import pandas as pd
		plasmid_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
		print(params)
		plasmid_table['host']=params.sample
		plasmid_table.to_csv(output.name, sep='\t')


def find_checkv_files(wildcards):
	path = "checkv_output/{sample}_completeness.name.tsv"
	files = glob(path.format(sample="*"))
	print(files)
	return files

rule aggregate:
	input:
		genome=find_checkv_files,
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

#rule get_fastas












