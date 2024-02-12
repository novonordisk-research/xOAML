# Code to generate polygenic risk scores (PRS) either genome-wide or at the gene/pathway levels
# Nielsen et al, Nat Comms 2024

# OA_GWAS_excluding_UKBB.txt: OA GWAS (sub-analyses excluding the UKBB) summary statistics shared by the authors of:
# Boer, C. G. et al. Deciphering osteoarthritis genetics across 826,690 individuals from 9 populations. Cell 184, 4784-4818.e17 (2021


################# PRSice (genome-wide) ###############

/nfs_home/users/tsmf/software/PRSice/PRSice.R \
--prsice /nfs_home/users/tsmf/software/PRSice/PRSice_linux \
--base input_data_path/OA_GWAS_excluding_UKBB.txt \
--snp SNP \
--chr CHR  \
--A1 EA  \
--A2 NEA  \
--or  \
--stat OR  \
--pvalue P  \
--target filtered_imputed_ukbb_data_path/ukb_imp_chr# \
--type bed \
--bar-levels 1 \
--fastscore \
--no-regress \
--all-score \
--print-snp \
--thread 10 \
--memory 400Gb \
--seed 42 \
--extract out_path/PRSice_result.valid \
--out out_path/PRSice_result 


################# PRSset (gene/pathway-levels) ###############

# optional parameters used
# --proxy 0.8 \

# tested_set.txt correspond to the outputs of the PRSet_prep.R script:
# Either a list of protein_coding genes (coding_genes.txt) or of pathways (msigdb_subset_for_PRS.txt)
# Subsets of the generated gene-level and pathway-levels risk scores were used as input for the ML models, as described in the paper's methods

# gencode.v43lift37.annotation.gtf from https://www.gencodegenes.org/human/release_39lift37.html

/nfs_home/users/tsmf/software/PRSice/PRSice.R \
--prsice /nfs_home/users/tsmf/software/PRSice/PRSice_linux \
--base input_data_path/OA_GWAS_excluding_UKBB.txt \
--snp SNP \
--chr CHR  \
--bp POS \
--A1 EA  \
--A2 NEA  \
--or \
--stat OR  \
--pvalue P  \
--target filtered_imputed_ukbb_data_path/ukb_imp_chr# \
--type bed \
--gtf input_data_path/gencode.v43lift37.annotation.gtf \
--msigdb input_data_path/tested_set.txt \
--no-regress \
--wind-3 1500 \
--wind-5 5000 \
--print-snp \
--thread 10 \
--seed 42 \
--memory 400Gb \
--extract out_path/PRSet_results.valid \
--out out_path/PRSet_results

