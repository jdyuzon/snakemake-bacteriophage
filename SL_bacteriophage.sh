#!/bin/bash
#SBATCH --account m342
#SBATCH --qos regular
#SBATCH --nodes 2
#SBATCH --constraint cpu
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH -t 02:00:00
#SBATCH -J bacteriophage

#####################
#####rule genomad:###
#####################
source /global/homes/j/jdyuzon/.bashrc
conda activate snakemake-bacteriophage


sample=$1

### Genomad
genomad end-to-end --min-score 0.7 --cleanup --splits 8 bacterial_genomes/${sample}.fna genomad_output/${sample} genomad_db

### CheckV
checkv end_to_end genomad_output/${sample}/${sample}_summary/${sample}_virus.fna checkv_output/${sample} -d checkv-db-v1.5 -t 16
cp checkv_output/${sample}/completeness.tsv checkv_output/${sample}_completeness.tsv

### MMseqs2
conda activate mmseqs2
mkdir mmseqs_target_seq/${sample}
cp genomad_output/${sample}/${sample}_summary/${sample}_virus_proteins.faa mmseqs_target_seq/${sample}/${sample}_virus_proteins.faa 
mmseqs createdb mmseqs_target_seq/${sample}/${sample}_virus_proteins.faa mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq 

### MMseqs2/Phrogs
mmseqs search phrogs_mmseqs_db/phrogs_profile_db \
mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs \
mmseqs_target_seq/${sample}/tmp -s 7

mmseqs createtsv phrogs_mmseqs_db/phrogs_profile_db \
mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs.tsv --full-header 

cp "mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs.tsv mmseqs_target_seq
echo "file: mmseqs_target_seq/${sample}_virus_proteins_mmseqs.tsv"

conda activate snakemake-bacteriophage
python3 bacteriophage_edit.py ${sample}

#rule tiger:
#	input:
#		genome="bacterial_genomes/{sample}.fna",
#		genome_link="tiger_output/{sample}/genome.fa"
#		database=""
#	output:
#		attsites="{sample}.islesFinal.gff"
#	shell:
#	"""
#ln -s bacterial_genomes/${sample}.fna tiger_output/${sample}/genome.fa
#perl /pscratch/sd/j/jdyuzon/BioSoft/TIGER-TIGER/test-tiger.pl -db {database} -fasta tiger_output/${sample}/genome.fa -nickname ${sample}
#perl /pscratch/sd/j/jdyuzon/BioSoft/TIGER-TIGER/bin/typing.pl tiger_output/${sample}/genome.island.nonoverlap.gff 
#perl /pscratch/sd/j/jdyuzon/BioSoft/TIGER-TIGER/resolve.pl tiger > tiger_output/${sample}/resolved.gff
#perl /pscratch/sd/j/jdyuzon/BioSoft/TIGER-TIGER/bin/typing.pl tiger_output/${sample}/resolved.gff
#mv tiger_output/${sample}/islesFinal.gff ${sample}.islesFinal.gff
#	"""


#rule checkv_edit:
#	input:
#		indf="checkv_output/{sample}_completeness.tsv"
#	output:
#		name="checkv_output/{sample}_completeness.name.tsv"
#	params:
#		sample="'{sample}'"
#	run:
#		checkv_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
#		print(params)
#		checkv_table['host']=params.sample
#		checkv_table.to_csv(output.name, sep='\t')
#
#


#rule phrog_edit:
#	input:
#		indf=lambda wildcards: "mmseqs_target_seq/{sample}_virus_proteins_mmseqs.tsv"
#	output:
#		indf="phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv"
#	params:
#		sample="'{sample}'"
#	run:
#		print("input "+input.indf)
#		print("param "+params.sample)
#		print("output "+output.indf)
#		phrog_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore',header=None)
#		phrog_table.columns=['#phrog','host_seq','alnScore','seqIdentity','eVal','qStart','qEnd','qLen','tStart','tEnd','tLen']
#		phrog_table[['#phrog','phrog_seq']] =phrog_table['#phrog'].str.split(' ## ',expand=True)
#		df_index = pd.read_csv('phrogs_mmseqs_db/PHROG_index.csv')
#		df_Bins_Index = pd.merge(phrog_table,df_index, how = 'left', on = '#phrog').fillna('NA')
#		df_Bins_Index['host']=params.sample
#		np.set_printoptions(threshold=sys.maxsize)
#		df_Bins_Index.to_csv(output.indf)


#rule plasmid_edit:
#	input:
#		indf="genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv"
#	output:
#		name="plasmid_output/{sample}_plasmid_summary.name.tsv"
#	params:
#		sample="'{sample}'"
#	run:
#		plasmid_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore')
#		print(params)
#		plasmid_table['host']=params.sample
#		plasmid_table.to_csv(output.name, sep='\t')


#rule aggregate:
#	input:
#		genome=glob("checkv_output/{sample}_completeness.name.tsv".format(sample="*")),
#		gene=glob("phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv".format(sample="*")),
#		plasmid=expand("plasmid_output/{sample}_plasmid_summary.name.tsv", sample=config["samples"],allow_missing=True)
#	output:
#		genome="checkv_output/all_genome_completeness.tsv",
#		gene="phrog_output/all_virus_proteins.csv",
#		plasmid="plasmid_output/all_plasmids.tsv"
#	shell:
#		"""
#		cat {input.genome} > {output.genome}
#		cat {input.gene} > {output.gene}
#		cat {input.plasmid} > {output.plasmid}
#		"""

#rule filter_viral_genomes:
#	input:
#		all="checkv_output/all_genome_completeness.tsv"
#	output:
#		high="checkv_output/all_genome_completeness.high.tsv",
#		low="checkv_output/all_genome_completeness.low.tsv"
#	run:
#		import pandas as pd
#		completeness = pd.read_csv(input.all, sep='\t')
#		#remove extra headers that got into main dataframe
#		completeness=completeness.loc[(completeness['contig_id'] != 'contig_id') & (completeness['contig_length'] != 'contig_length')]
#		#fill NAs for aai_completeness and aai_confidence, and make columns numeric
#		completeness['aai_completeness']=completeness['aai_completeness'].fillna(0)
#		completeness['aai_confidence']=completeness['aai_confidence'].fillna('low')
#		completeness['aai_completeness']=pd.to_numeric(completeness['aai_completeness'])
#		#higher than 95% aai_completeness and goes to high quality table
#		completeness_high=completeness.loc[(completeness['aai_confidence'] == 'high') & (completeness['aai_completeness'] >= 95)] 
#		#less than 95% aai_completeness or 95% completeness but hot high confidence goes to low quality table
#		completeness_low=completeness.loc[((completeness['aai_confidence'] != 'high')&(completeness['aai_completeness'] >= 95))|(completeness['aai_completeness'] < 95)]
#		#write tables
#		completeness_high.to_csv(output.high,sep='\t')
#		completeness_low.to_csv(output.low,sep='\t')

#rule filter_viral_genes:
#	input:
#		all="phrog_output/all_virus_proteins.csv"
#	output:
#		high="phrog_output/all_virus_proteins.high.csv",
#		low="phrog_output/all_virus_proteins.low.csv"
#	run:
#		genes = pd.read_csv(input.all)
#		#remove extra headers that got into main dataframe
#		genes=genes.loc[(genes['#phrog'] != '#phrog') & (genes['host_seq'] != 'host_seq')]
#		#e-value cutoff of 10E-6 or less (depends on the length of the query sequence and size of database. probability of occurance by chance alone 
#		genes_high=genes.loc[pd.to_numeric(genes['eVal'])<=1e-6]
#		genes_low=genes.loc[pd.to_numeric(genes['eVal'])>1e-6]
#		#write tables
#		genes_high.to_csv(output.high,sep='\t')
#		genes_low.to_csv(output.low,sep='\t')


#checkpoint filter_plasmid_genes:
#	input:
#		all="plasmid_output/all_plasmids.tsv"
#	output:
#		high="plasmid_output/all_plasmids.high.tsv",
#		low="plasmid_output/all_plasmids.low.tsv"
#	run:
#		plasmid = pd.read_csv(input.all,sep='\t')
#		plasmid=plasmid.loc[(plasmid['seq_name'] != 'seq_name') & (plasmid['length'] != 'length')]
#		plasmid_high=plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])>=0.95]
#		plasmid_low=plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])<0.95]
#		plasmid_high.to_csv(output.high,sep='\t')
#		plasmid_low.to_csv(output.low,sep='\t')

#rule get_fasta
