from glob import glob
from collections import defaultdict
import sys

"""
the data assigned to real_data below are one of the results files from the gradient forest output, which looks like this:

snpid,rsquare
chr1_12500530,0.297233861141119
chr1_13695307,0.130011209071222
chr1_13793572,0.00567890566778528

...
...

the random

real_data = open(sys.argv[1],'r')    # e.g. FULLDATASET/ALLSNPNLLids.csv'
real_data.readline()



real_snp_dict = {}

for line in real_data:
    snpid,rsquare =  line.replace('"','').strip().split()
    #snpid,rsquare =  line.strip().split(',') # for ellioti
    rsquare = float(rsquare)
    real_snp_dict[snpid] = {'rsquare': rsquare,'rand_greater_rsquare':0}

random_sig_snps_dict = defaultdict(int)


"""
this script assumes data for randomly permuted runs are in the directory where the script is launched.
the random runs do not have a header and enclose the snpid (concatenated chromosome and position) in quotes, e.g.:

"chr1_30781607",0.118502677192657
"chr1_41014695",0.232063166947916
"chr1_60903019",0.0632163580044689
...
...
"""



randoms = glob('rand*csv')
for random in randoms:
    randopen = open(random,'r')
    rand_dict = {}
    for line in randopen:
        snpid,rsquare =  line.replace('"','').strip().split(',')
        rsquare = float(rsquare)
        if snpid in real_snp_dict and rsquare > real_snp_dict[snpid]['rsquare']:
            real_snp_dict[snpid]['rand_greater_rsquare']+=1
        else:
            #print(line)
            random_sig_snps_dict[snpid]+=1


fout = open('gradient_forest_randomenvs_snp_pvalues.tsv','w')
fout.write('snpid\tpvalue\n')
for snp in real_snp_dict:
        fout.write('{}\t{}\n'.format(snp,real_snp_dict[snp]['rand_greater_rsquare']/200))
   
print(min(random_sig_snps_dict.values())/200)
 
fout.close()

