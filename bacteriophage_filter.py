from glob import glob
import pandas as pd
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

### filter_viral_genomes:
completeness = pd.read_csv("checkv_output/all_genome_completeness.tsv", sep='\t')
#remove extra headers that got into main dataframe
completeness=completeness.loc[(completeness['contig_id'] != 'contig_id') & (completeness['contig_length'] != 'contig_length')]
#fill NAs for aai_completeness and aai_confidence, and make columns numeric
completeness['aai_completeness']=completeness['aai_completeness'].fillna(0)
completeness['aai_confidence']=completeness['aai_confidence'].fillna('low')
completeness['aai_completeness']=pd.to_numeric(completeness['aai_completeness'])
#higher than 95% aai_completeness and goes to high quality table
completeness.loc[(completeness['aai_confidence'] == 'high') & (completeness['aai_completeness'] >= 95), 'quality'] = "high"
completeness.loc[((completeness['aai_confidence'] != 'high')&(completeness['aai_completeness'] >= 95))|(completeness['aai_completeness'] < 95),'quality']="low"
completeness.to_csv("all_genome_completeness.qual.tsv",sep='\t')

### Y-variable: Annotate completeness file with TIGER (Phage 1/Phage 2), T6SS, eCIS, tailocin using bedtools

### Y-variable: Merge Annotated file with Phrog file

### rule filter_viral_genes:
genes = pd.read_csv("phrog_output/all_virus_proteins.csv")
#remove extra headers that got into main dataframe
genes=genes.loc[(genes['#phrog'] != '#phrog') & (genes['host_seq'] != 'host_seq')]
#e-value cutoff of 10E-6 or less (depends on the length of the query sequence and size of database. probability of occurance by chance alone 
genes.loc[pd.to_numeric(genes['eVal'])<=1e-6,'quality']="high"
genes.loc[pd.to_numeric(genes['eVal'])>1e-6,'quality']="low"
genes.to_csv("phrog_output/all_virus_proteins.qual.csv",sep='\t')

### get contig column
contig_id=genes['host_seq'].replace({' #.*': ''}, regex=True)
contig_id=contig_id.str.rsplit("_", n = 1, expand = True)[0]
genes['contig_id']=contig_id
### merge genomad and phrog data
completeness_genes=pd.merge(completeness,genes, how = 'left', on = ['contig_id','host'], suffixes=('_virus', '_gene')).fillna('NA')
completeness_genes.to_csv("phrog_output/all_genome_completeness.qual.tsv",sep='\t')

plasmid = pd.read_csv("plasmid_output/all_plasmids.tsv",sep='\t')
plasmid=plasmid.loc[(plasmid['seq_name'] != 'seq_name') & (plasmid['length'] != 'length')]
plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])>=0.95, 'quality']="high"
plasmid.loc[pd.to_numeric(plasmid['plasmid_score'])<0.95, 'quality']="low"
plasmid.to_csv("plasmid_output/all_plasmids.qual.tsv",sep='\t')
