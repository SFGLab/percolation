import sys
import numpy as np

def values(waointersect, data_name, output):
    length={}
    sl={}
    with open(sys.argv[1]) as fileobject: 
        for linee in fileobject:
            line = linee.strip('\n').split('\t')
            name = line[0]+'_'+line[1]+'_'+line[2]
            length[name] = int(line[2])-int(line[1])
            if (data_name.find('merged') != -1):
                value0 = float(line[6]) #length of overlap
            else:
                value0 = float(line[7]) #length of overlap
            if line[6] == '.': #loops without overlap get score 0
                line[6] = 0
            value1 = float(line[6]) #score
            value2 = value1 * value0 #score * length of overlap = sum of a peak
            if name in sl:
                sl[name][0].append(value0) #length of overlap
                sl[name][1].append(value1) #score
                sl[name][2].append(value2) #sum of a peak
            else:
                sl[name]=[[value0],[value1],[value2]]

    w = open(output,'w')
    w.write('chr\tstart\tend\t'+data_name+'_mean'+'\t'+data_name+'_sum'+'\t'+data_name+'_fraction'+'\t'+data_name+'_min'+'\t'+data_name+'_max'+'\t'+data_name+'_mean.pbp'+'\n')
    for i in sl:
        p=i.split('_')
        #mean of the peak signal, sum of the peak signal, fraction of the loop that is covered by peaks, minimal peak, maximal peak
        w.write(p[0]+'\t'+p[1]+'\t'+p[2]+'\t'+str(np.mean(sl[i][1]))+'\t'+str(np.sum(sl[i][1]))+'\t'+str(np.sum(sl[i][0])/float(length[i]))+'\t'+str(np.min(sl[i][1]))+'\t'+str(np.max(sl[i][1]))+'\t'+str(np.sum(sl[i][2])/float(length[i]))+'\n')
    w.close() 
    
values(sys.argv[1], sys.argv[2], sys.argv[3])

#-wao	Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r.
