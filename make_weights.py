# make a weights file of 1 and 0 
# this is used to select only positive examples for computing the spearmand and the pearson correlations in the predictAndEval.py class 
data=open('diffMatAnna.binarized.combined.stringent.csv','r').read().strip().split('\n') 
outf=open('remapped.stringent.weights.txt','w') 
tasks=['CC','3hr','16hr','48hr','H1','Hk','M5','X3hr.CC','X16hr.3hr','X16hr.CC','X48hr.16hr','X48hr.3hr','X48hr.CC','H1.48hr','H1.16hr','H1.3hr','H1.CC','H1.Hk','X48hr.Hk','X16hr.Hk','X3hr.Hk','H1.M5','X48hr.M5','X16hr.M5','X3hr.M5'] 
header=data[0].split('\t') 
cur_peak=0 
outf.write('peak\t'+'\t'.join(tasks)+'\n')

for line in data[1::]: 
    tokens =line.split('\t') 
    peak_dict=dict() 
    for i in range(3,len(tokens)): 
        cur_header=header[i].replace('Up','').replace('Down','') 
        if cur_header not in peak_dict: 
            peak_dict[cur_header]=[int(tokens[i])] 
        else: 
            peak_dict[cur_header].append(int(tokens[i])) 
    outf.write('peak_'+str(cur_peak)) 
    for task in tasks: 
        taskval=max([abs(i) for i in peak_dict[task]])
        outf.write('\t'+str(taskval))
    outf.write('\n') 
    cur_peak+=1 
    
        
