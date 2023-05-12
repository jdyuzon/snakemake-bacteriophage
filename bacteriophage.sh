#!/bin/bash
#SBATCH --account m342
#SBATCH --qos regular
#SBATCH --nodes 2
#SBATCH --constraint cpu
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH -t 04:00:00
#SBATCH -J bacteriophage
###
for genome in bacterial_genomes/*fna
do
        genome=${genome%.fna}
        genome=${genome#bacterial_genomes/}
        echo $genome

        bash SL_bacteriophage.sh $genome 
#       sbatch SL_bacteriophage.sh $genome 
done

###aggregate:
cat checkv_output/*_completeness.name.tsv > checkv_output/all_genome_completeness.tsv
cat phrog_output/*_virus_proteins_mmseqs_Annotated.csv > phrog_output/all_virus_proteins.csv
cat plasmid_output/*_plasmid_summary.name.tsv > plasmid_output/all_plasmids.tsv
