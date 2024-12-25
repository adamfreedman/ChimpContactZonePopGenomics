# EEMS
To make spatially-explicit inferences regarding where on the landscape gene flow is restricted, we use Estimating Effective Migration Surfaces (EEMS). The code repository can be found [here](https://github.com/dipetkov/eems) with methdological details described in [Petkova *et al.* 2016, *Nature Genetics*](https://www.nature.com/articles/ng.3464). 


Before running the main migration surface module, we perform the following initial steps:

## 1. Generate plink files from a vcf file, filtered to exclude sites with a minor allele frequency less than 0.05. For example, with the neutral SNPs, we do:

```
plink --vcf maf05_Phase_3_CHIMP_10k_NEUTRAL.vcf.recode.vcf --out maf05_Phase_3_CHIMP_10k_NEUTRAL --allow-extra-chr
```


## 2. create EEMS input coordinate file
This can be done with scripts/WriteEemsCoordFile.py
```
python WriteEemsCoordFile.py
```

## 3. fix bim file so doesn't have non standard chromosomes
```
awk '{print 1,$2,$3,$4,$5,$6}' datapath.bim  > newdatapath.bim
```

## .4 run bed2iffs
```
/n/home14/afreedman/software/eems/bed2diffs/src/bed2diffs_v1 --bfile maf05_Phase_3_CHIMP_10k_NEUTRAL
```

Parameter files used in running eems with scripts/RunEems.sh are found in the paramter_files directory. We analyze neutral and outlier SNPs separately, with these two categories of SNPs on the capture array being determined based upon whole-genome selection scans performed with captive chimpanzees from Cameroon. To account for the effects of varying the deme size setting, for each set of SNPs we generated five replicates for each of three deme sizes: 100, 200, and 400. Because the run with deme size 100 did not converge with the other two, and displayed substantially lower log-likelihoo traces for the MCMC, this deme size was excluded when integrating results across replicates and deme sizes.

## 5. Run EEMS
We have provided an example SLURM sbatch script to do this in scripts, which can simply be run on a cluster using the SLURM scheduler as follows:
```
sbatch RunEems.sh
```
where the primary command line argument is the parameter file that specifies the path to the data and other settings. 
