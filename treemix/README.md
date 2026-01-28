# TreeMix
Scripts to infer population splits and mixture events from alelle frequency data using TreeMix by Pickrell & Pritchard (2012). This pipeline runs TreeMix with bootstrapping, helps choose number of migration events and creates a consensus tree. It plots a maximum likelihood tree with bootstrap values, drift and residuals and calculates statistics for every migration event.

# Citation for generating TreeMix input file: Thom Nelson's adaptation (2019) of Tomaz Berisa's script to convert plink clustered allele frequencies to treemix format
# Citation of TreeMix pipeline: Dahms, C. (2021) (All scripts and codes)

# Pipeline
# 1. Generate TreeMix Input file
#In commandline cluster, the following code was run to generate input file
module spider vcf tools # load required module
$ module load CC/11.2.0
$ module load GCC/11.2.0
$ module load VCFtools/0.1.16
$ vcftools --vcf file.vcf --plink-tped --out plinkfileoutput # convert vcf to plink format
$ module purge # close all modules
$ module spider PLINK $module spider PLINK/2.00a3.7
$ module load GCC/12.3.0
$ module load  PLINK/2.00a3.7
$ plink --tfile plinkfileoutput --freq --within poporder.list # calculate allele frequency based on poporder.list which contains all the population names as rows
$ gzip plink.frq.strat #compress the output file
$ module purge
$ module load GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
$ module load  Biopython/1.75-Python-3.7.4
$ python2 plink2treemix.py plink.frq.strat.gz treemix.gz #use the python script to generate TreeMix input file

# 2. Build consensus tree with multiple migration events
# Run `Step1_TreeMix.sh` in the command line, providing an input file, maximum number of cores, block size, outgroup (or 'noRoot' for unrooted trees), 
number of bootstrap replicates, path to PHYLIP consense program, output file name, range of migration events (m) and their number of replicates, for example:

$ sh Step1_TreeMix.sh treemix_neutral.gz 10  100 noRoot 500 /sw/eb/sw/PHYLIP/3.697-GCC-9.3.0/bin/consense neutral_18pop 1 10 10 
# in our case the input file is treemix_neutral.gz and we had 10 cores, 100 block size, noRoot as unrooted tree, the specified path for PHYLIP consense program, neutral_18pop as output file, 10 migration events and 10 replicates

#This builds a consensus tree from bootstraps and adds a specified range of m. 
#Tree replicates will be stored in the *test_migrations* folder.

# 3. Test migration edges with OptM in R

#Set working directory to the *test_migrations* folder and run the R package OptM (step A) from the R script `Step2&4_TreeMix.R`.
#This helps identify the optimum number of m.

# 4. Final runs with optimum number of migration edges in terminal/cluster

#Run `Step3_TreeMix.sh` by providing an input file, maximum number of cores, block size, outgroup (alternatively 'noRoot' for unrooted trees), number of bootstrap replicates, number of migrations, output file name, number of independent runs (N), name of consensus tree built in Step 2, and path to consense program, -noss if there are very few samples in some population to reduce bias

$ sh Step3_TreeMix.sh treemix_neutral.gz 10 100 noRoot 500 5 5neutral_pop18 30 neutral_15pop_constree.newick /sw/eb/sw/PHYLIP/3.697-GCC-9.3.0/bin/consense -noss

#Returns trees from chosen number of independent runs with optimum number of m. 
#Final runs of trees will be stored in the *final_runs* folder.

#5. Tree visualization + Migration stats and support 

For this step you will need to have saved the file `TreeMix_functions.R`
Set working directory to the *final_runs* folder, run steps B and C from the `Step2&4_TreeMix.R` script.
From the final runs, compares tree likelihoods, plots ML tree with bootstrap values and migration weights, as well as drift and residuals.

References
Milanesi, M., Capomaccio, S., Vajana, E., Bomba, L., Garcia, J.F., Ajmone-Marsan, P., Colli, L., 2017. BITE: an R package for biodiversity analyses. bioRxiv 181610. doi:10.1101/181610

Pickrell, J., & Pritchard, J. (2012). Inference of population splits and mixtures from genome-wide allele frequency data. Nature Precedings, 1-1.

Zecca, G., Labra, M., & Grassi, F. (2020). Untangling the Evolution of American Wild Grapes: Admixed Species and How to Find Them. Frontiers in Plant Science, 10, 1814.

Dahms, C., Kemppainen, P., Zanella, L. N., Zanella, D., Carosi, A., Merilä, J., & Momigliano, P. (2022). Cast away in the Adriatic: Low degree of parallel genetic differentiation in three‐spined sticklebacks. Molecular Ecology, 31(4), 1234-1253.
