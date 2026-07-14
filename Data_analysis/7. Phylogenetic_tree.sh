### TREEMIX analysis ###
Script by Ariel Rodríguez



# Merge VCF files of the three species: Oophaga granulifera and two outgroups (Oophaga sylvatica and Dendrobates tinctorius)

# This section will extract unlinked SNPs
plink2 --vcf granulifera_tintorius_jointCalling_filteredSNPs_OK.vcf --allow-extra-chr --set-all-var-ids @:#\$r,\$a --remove outgroup_samples.txt --indep-pairwise 50 1 0.7 --bad-ld 

plink2 --vcf granulifera_tintorius_jointCalling_filteredSNPs_OK.vcf --allow-extra-chr --set-all-var-ids @:#\$r,\$a --exclude plink2.prune.out --snps-only --max-alleles 2 --geno 0.1 --remove outgroup_samples.txt --make-bed --recode --out granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK

plink --bfile granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK --allow-extra-chr --recode A --out granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK

plink --bfile granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK --allow-extra-chr --recode --tab --out granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK

# Prepare Treemix input
plink --file granulifera_jointCalling_filteredSNPs_OK.LDprunned.PLINK \
--freq --allow-extra-chr --within populations.granu_ingroup.txt --out granulifera_ingroup_LDprunned

# Compress the output above [unlinked/linked SNPs] and feed it to the python converter
rm granulifera_ingroup.frq.strat.gz
gzip -k granulifera_ingroup_LDprunned.frq.strat
python /data/treemix-1.13/plink2treemix.py granulifera_ingroup_LDprunned.frq.strat.gz granulifera_LDpruned_treemix.gz

# The above is not so important for treemix, here is an alternative with all ingroup SNPs
plink2 --vcf granulifera_tintorius_jointCalling_filteredSNPs_OK.vcf --allow-extra-chr --set-all-var-ids @:#\$r,\$a --snps-only --max-alleles 2 --geno 0.01 --remove outgroup_samples.txt --make-bed --recode --out granulifera_jointCalling_filteredSNPs_OK.PLINK

plink --bfile granulifera_jointCalling_filteredSNPs_OK.PLINK --allow-extra-chr --recode A --out granulifera_jointCalling_filteredSNPs_OK.PLINK

plink --bfile granulifera_jointCalling_filteredSNPs_OK.PLINK --allow-extra-chr --recode --tab --out granulifera_jointCalling_filteredSNPs_OK.PLINK

plink --file granulifera_jointCalling_filteredSNPs_OK.PLINK \
--freq --allow-extra-chr --within populations.granu_ingroup.txt --out granulifera_ingroup_all

rm granulifera_ingroup_all.frq.strat.gz
gzip -k granulifera_ingroup_all.frq.strat
python /data/treemix-1.13/plink2treemix.py granulifera_ingroup_all.frq.strat.gz granulifera_treemix_all.gz

# Also converting to genepop format for divMigrate analyses.

python /data/treemix-1.13/convert.py

# Running TreeMix
cp /data/granulifera_w_outgroups/CASTERsite_granulifera_species_tree.branchlength.tre species_tree.tre
# Since LD prunning was previously done, no need for -k parameter
# since tree is known, no need for -bootstrap param
mkdir test_migrations
for m in {0..5}
   do
   for i in {1..20}
      do
      treemix -i granulifera_LDpruned_treemix.gz -root PAL -o ./test_migrations/treemix.${m}.${i} \
         -bootstrap -k 500 -global -m ${m} -noss -se 
      done 
done

# Open results in R 

setwd("/data/granulifera/all_tissues/treemix/test_migrations")
source("/data/treemix-1.13/src/plotting_funcs.R")
source("/data/treemix-1.13/TreeMix_Explorer/TreeMix-main/TreeMix_functions.R")

# (A) Test migration events #	: https://github.com/carolindahms/TreeMix/blob/main/Step2%264_TreeMix.R
# install.packages("OptM")
# devtools::install_local("/data/treemix-1.13/BITEV2_2.1.2.tar.gz")

library(OptM)
library(dplyr)
library(data.table)
library(BITEV2)
library(plyr)

folder <- "/data/granulifera/all_tissues/treemix/test_migrations"  #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear", plot=TRUE, pdf="Optimal_number_of_Migrations_linearMethods.OK.pdf")    #shows changes in log likelihood for different m, and suggests optimum m as 'change points'
test.optM = optM(folder, method = "Evanno", skip=0, tsv ="Evanno.variance.txt")  #another option is the Evanno method - see optM package description for detailed information on output

#if data is robust and all runs have the same likelihoods, SD will be 0 and this function will give an error as it can't produce the ad hoc statistic. 
#in this case you might want to increase variance by varying -k (SNP block size), change permutation methods etc.

plot_optM(test.optM, method = "Evanno")     #plot the proportion of variation explained by each migration event. Calculates deltaM, which is a second-order rate of change in likelihood weighted by the standard deviation
quit() 

treemix -i granulifera_LDpruned_treemix.gz -root PAL -o granulifera_m-1 -m 1 -noss -tf species_tree_ingroup.tre | tee treemix.log
     
# Final plots and stats
R
setwd("/data/granulifera/all_tissues/treemix")

pdf("granulifera_treemixoutput_M-1.pdf")
plot_tree("granulifera_m-1", ybar=0.5, flip=0, arrow=0.2, lwd=2)
dev.off()

pdf("TreeMix_Drift.pdf") 
treemix.drift(in.file = "granulifera_m-1",     #pairwise matrix for drift estimates with specified order of pops 
              pop.order.color.file = "pop.colors.csv") + 
              title("Drift")     
dev.off()


plot_resid("finalrun_1",       #pairwise matrix for residuals
           pop_order = "poporder.txt") +
           title("Residuals")                         
dev.off()

