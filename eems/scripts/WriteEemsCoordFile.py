famin = open('maf05_Phase_3_CHIMP_10k_NEUTRAL.fam','r')
coordout = open('maf05_Phase_3_CHIMP_10k_NEUTRAL.coord','w')
envdata = open('Phase_3_10k_coords_env.csv','r')

env_dict = {}
envdata.readline()
for line in envdata:
    linelist = line.strip().split(',')
    env_dict[linelist[0]] = {'X': linelist[3],'Y': linelist[2]}


order_out = open('maf05_Phase_3_CHIMP_10k_NEUTRAL.order','w')

for line in famin:
    id=line.strip().split()[0]
    order_out.write('%s %s\n' % (id,id))
    coordout.write('%s %s\n' % (env_dict[id]['X'],env_dict[id]['Y']))

coordout.close()
order_out.close()
