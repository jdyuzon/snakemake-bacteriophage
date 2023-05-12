from glob import glob
import pandas as pd
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import sys

print 'sample:', str(sys.argv)

### checkv_edit
checkv_table=pd.read_csv("checkv_output/{sample}_completeness.tsv",sep = '\t', encoding_errors='ignore')
print("'{sample}'")
checkv_table['host']="'{sample}'"
checkv_table.to_csv("checkv_output/{sample}_completeness.name.tsv", sep='\t')


### phrog_edit
print("mmseqs_target_seq/{sample}_virus_proteins_mmseqs.tsv "+input.indf)
print("'{sample}' "+params.sample)
print("phrog_output/{sample}_virus_proteins_mmseqs_Annotated.csv "+output.indf)
phrog_table=pd.read_csv(input.indf,sep = '\t', encoding_errors='ignore',header=None)
phrog_table.columns=['#phrog','host_seq','alnScore','seqIdentity','eVal','qStart','qEnd','qLen','tStart','tEnd','tLen']
phrog_table[['#phrog','phrog_seq']] =phrog_table['#phrog'].str.split(' ## ',expand=True)
df_index = pd.read_csv('phrogs_mmseqs_db/PHROG_index.csv')
df_Bins_Index = pd.merge(phrog_table,df_index, how = 'left', on = '#phrog').fillna('NA')
df_Bins_Index['host']=params.sample
np.set_printoptions(threshold=sys.maxsize)
df_Bins_Index.to_csv(output.indf)


### plasmid_edit
plasmid_table=pd.read_csv("genomad_output/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv",sep = '\t', encoding_errors='ignore')
print("'{sample}'")
plasmid_table['host']="'{sample}'"
plasmid_table.to_csv("plasmid_output/{sample}_plasmid_summary.name.tsv", sep='\t')

