# snakemake-bacteriophage
snakemake pipeline identifies viral genomes, viral genes, and plasmids from bacterial genomes

## 1. Software requirements:
channels:
  - bioconda
  - conda-forge

dependencies:
  - snakemake-minimal >=5.24.1
  - geNomad >=1.5.0
  - mmseqs2 >=13.45111
  - checkv >=1.0.1
  - blast >=2.6.0
  - prokka =1.11
  - pandas
  - bedtools =2.27.1
  - trnascan-se >= 2.0.2
 
Set up the environment:
```
conda env create --name snakemake-bacteriophage --file environment.yaml
```

## 2. List of bacterial genomes to query:
```
less config.yaml
```
For example the config file looks like:
samples:
     Agrobacterium_tumefaciens_12D1_IMG_2846775462: Agrobacterium_tumefaciens_12D1_IMG_2846775462
     Bradyrhizobium_japonicum_USDA-110_IMG_637000038: Bradyrhizobium_japonicum_USDA-110_IMG_637000038
     Bradyrhizobium_sp.-ARR65_IMG_2508501128: Bradyrhizobium_sp.-ARR65_IMG_2508501128
     Candidatus-Liberibacter_asiaticus-str.-gxpsy_IMG_2540341131: Candidatus-Liberibacter_asiaticus-str.-gxpsy_IMG_2540341131
     Clavibacter_michiganensis-michiganensis_CAYO001_IMG_2854244720: Clavibacter_michiganensis-michiganensis_CAYO001_IMG_2854244720
     Dickeya_solani_IFB0223_IMG_2847231513: Dickeya_solani_IFB0223_IMG_2847231513
     Ensifer_meliloti_2011_IMG_2562617130: Ensifer_meliloti_2011_IMG_2562617130
     Escherichia_coli_ST131-ErtS_IMG_2690315815: Escherichia_coli_ST131-ErtS_IMG_2690315815
     Paraburkholderia_kururiensis_M130_IMG_2551306521: Paraburkholderia_kururiensis_M130_IMG_2551306521
     Pectobacterium_atrosepticum_SCRI1043_IMG_637000102: Pectobacterium_atrosepticum_SCRI1043_IMG_637000102
     Pseudomonas_aeruginosa_PAO1161_IMG_2917026592: Pseudomonas_aeruginosa_PAO1161_IMG_2917026592
     Pseudomonas_putida_KT2440_ASM756v2: Pseudomonas_putida_KT2440_ASM756v2
     Pseudomonas_simiae_WCS417_IMG_2585427642: Pseudomonas_simiae_WCS417_IMG_2585427642
     Pseudomonas_syringae-pv.-tomato-str.-DC3000_ASM780v1: Pseudomonas_syringae-pv.-tomato-str.-DC3000_ASM780v1
     Ralstonia_solanacearum_UA-1612_IMG_2841725018: Ralstonia_solanacearum_UA-1612_IMG_2841725018
     Rathayibacter_toxicus_FH128_IMG_2858033217: Rathayibacter_toxicus_FH128_IMG_2858033217
     Rhizobium_etli-bv.-mimosae_IE4771_IMG_2585427632: Rhizobium_etli-bv.-mimosae_IE4771_IMG_2585427632
     Rhizobium_leguminosarum-bv.-viciae_3841_IMG_639633055: Rhizobium_leguminosarum-bv.-viciae_3841_IMG_639633055
     Xanthomonas_axonopodis-pv.-vasculorum_NCPPB-900_IMG_2636416122: Xanthomonas_axonopodis-pv.-vasculorum_NCPPB-900_IMG_2636416122


## 3. Snakefile shows the pipeline:
Proviruses are identified by genomad
CheckV returns complete and clean phage genomes
MMseqs queries the viral proteins against the Phrog database
Results are aggregated into a file for provirus genomes, viral genes, and plasmids

An example command:
```
snakemake --cluster "sbatch -n 1 --qos=regular --constraint=cpu -t 1:00:00" --default-resources --jobs 1 --rerun-incomplete --keep-going
```

## 4. SL_snakemake-bacteriophage.sh
The pipeline can also be run on slurm
```
sbatch SL_snakemake-bacteriophage.sh
```

