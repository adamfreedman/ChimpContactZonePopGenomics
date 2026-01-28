from glob import glob
import sys
snpid = sys.argv[1]
gene = sys.argv[2]

def BuildPopFreqDict(frqfile):
    fopen = open(frqfile,'r')
    snp_dict = {}
    fields = ['chrom','pos','n_alleles','n_chr','frq_1','frq_2']
    fopen.readline()
    for line in fopen:
        line_dict = dict(zip(fields,line.strip().split()))
        snp_dict['{}:{}'.format(line_dict['chrom'],line_dict['pos'])] = line_dict['frq_1']
    return snp_dict 

sort_order = ['Bankim','DengDeng','Dja','Ebo',
              'LiabelemHighlands', 'Batanga', 'Lobeke',
              'MountCameroon', 'MbamDjerem', 'MountGolep',
              'Minta','Sabongida','Takamanda', 'Wouchaba',
              'Yagba', 'Mone']

sites = open('filtered_site_latlongs.tsv','r')
sites_dict = {}
fields = sites.readline().strip().split()
for site in sites:
    line_dict = dict(zip(fields,site.strip().split()))
    sites_dict[line_dict['Location']] = {'Lat' : line_dict['Lat'],'Long' : line_dict['Long']}


snp_dict = {}
for site in sites_dict:
    if site not in sort_order:
        raise ValueError('{} not in sort order list'.format(site))
    if site == "LiabelemHighlands":
        site = 'Liabelem'
    frqfile = glob("frqfiles/population*{}*.frq".format(site))[0]
    snp_dict["{}".format(site)] = BuildPopFreqDict(frqfile)


fout = open('{}_{}_pop_allele_frequencies.csv'.format(snpid,gene),'w')
fout.write('population,lat,long,snpid,gene,allele_freq\n')


for population in sort_order:
    if population == 'LiabelemHighlands':
        poplabel = 'Liabelem'
    else:
        poplabel = population
    fout.write('{},{},{},{},{},{}\n'.format(population,sites_dict[population]['Lat'],sites_dict[population]['Long'],
               snpid,gene,snp_dict[poplabel][snpid]))

fout.close()     
