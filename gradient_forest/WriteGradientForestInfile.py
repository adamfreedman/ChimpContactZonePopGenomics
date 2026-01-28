from collections import defaultdict
from os import system
import sys
numsnps=int(sys.argv[1]) #7878

path='/n/holylfs05/LABS/informatics/Users/afreedman/africa/chimps/allele_frequency/allsites/'
bankim=open('%spopulation_Bankim_2024.10.29.frq' % path,'r')
dengdeng=open('%spopulation_DengDeng_2024.10.29.frq' % path,'r')
dja=open('%spopulation_Dja_2024.10.29.frq' % path,'r')
ebo=open('%spopulation_Ebo_2024.10.29.frq' % path,'r')
liabelem=open('%spopulation_LiabelemHighlands_2024.10.29.frq' % path,'r')
batanga=open('%spopulation_Batanga_2024.10.29.frq' % path,'r')
lobeke=open('%spopulation_Lobeke_2024.10.29.frq' % path,'r')
mountcameroon=open('%spopulation_MountCameroon_2024.10.29.frq' % path,'r')
mbamdjerem=open('%spopulation_MbamDjerem_2024.10.29.frq' % path,'r')
mountgolep=open('%spopulation_MountGolep_2024.10.29.frq' % path,'r')
minta=open('%spopulation_Minta_2024.10.29.frq' % path,'r')
sabongida=open('%spopulation_Sabongida_2024.10.29.frq' % path,'r')
takamanda=open('%spopulation_Takamanda_2024.10.29.frq' % path,'r')
wouchaba=open('%spopulation_Wouchaba_2024.10.29.frq' % path,'r')
yagba=open('%spopulation_Yagba_2024.10.29.frq' % path,'r')
mone=open('%spopulation_Mone_2024.10.29.frq' % path,'r')


handles=[bankim,dengdeng,dja,ebo,liabelem,
         batanga,lobeke,mountcameroon,mbamdjerem,
         mountgolep,minta,sabongida,takamanda,
         wouchaba,yagba,mone]

for handle in handles:
    handle.readline()

snplist = []

fout=open('chimps_allsites_withMinta_2024.10.28.inp','w')
snpdict=defaultdict(list)
for i in range(numsnps):
    print('processing snp %s' % i)
    allele_dict=defaultdict(list)
    for handle in handles:
        handlelist=handle.readline().strip().split()
        snpid="%s_%s" % (handlelist[0],handlelist[1]) # id is concatenation of chromosome and position
        if snpid not in snplist:
            snplist.append(snpid)
        #snpdict[snpid].append(handlelist[4].split(':')[1])
        snpdict[snpid].append(handlelist[4])

labels=['bankim','dengdeng','dja','ebo','liabelem',
         'batanga','lobeke','mountcameroon','mbamdjerem',
         'mountgolep','minta','sabongida','takamanda',
         'wouchaba','yagba','mone']
for i in range(len(labels)):
    freqs=[]
    for key in snpdict.keys():
        freqs.append(snpdict[key][i])
            
    fout.write('%s\t%s\n' % (labels[i],'\t'.join(freqs)))

header=open('header.txt','w')
header.write('POP\t%s\n' % '\t'.join(snplist))

header.close()
fout.close()
system('cat header.txt chimps_allsites_withMinta_2024.10.28.inp > wheader_chimps_allsites_withMinta_2024.10.28.inp')
