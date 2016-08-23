low_thresh=0.1
medium_thresh=0.20 

#previous input file 
#peak_entries=open('/srv/scratch/annashch/stemcells/het/anna_all_timepoint_pairs/aggregate_pangwei.adjusted.binarized.filtered.tsv','r').read().split('\n') 
peak_entries=open('diffMatAnna.binarized.combined.stringent.csv','r').read().split('\n') 
while '' in peak_entries: 
    peak_entries.remove('') 

#confidence-file: 
#confidence=open('/srv/scratch/annashch/stemcells/het/confidenceMatAnna_padjCutoff-0.1_peakCutoff-20.csv','r').read().split('\n') 
confidence=open('/srv/scratch/annashch/stemcells/het/confidenceMatAnna.csv','r').read().replace('\"','').split('\n') 
while '' in confidence: 
    confidence.remove('') 


#fold change file so we know if something is upregulated or down-regulated 
#foldchange=open('/srv/scratch/annashch/stemcells/het/foldChangeMatAnna_padjCutoff-0.1_peakCutoff-20.csv','r').read().split('\n') 
foldchange=open('/srv/scratch/annashch/stemcells/het/foldChangeMatAnna.csv','r').read().replace('\"','').split('\n') 
while '' in foldchange: 
    foldchange.remove('') 

#map the header of one file to the header of the other file, the naming conventions are different 
header_map=dict() 
header_map={"X3hr.CCUp":"X3hr.CC",
            "X3hr.CCDown":"X3hr.CC",
            "X16hr.3hrUp":"X16hr.3hr",
            "X16hr.3hrDown":"X16hr.3hr",
            "X16hr.CCUp":"X16hr.CC",
            "X16hr.CCDown":"X16hr.CC",
            "X48hr.16hrUp":"X48hr.16hr",
            "X48hr.16hrDown":"X48hr.16hr", 
            "X48hr.3hrUp":"X48hr.3hr",
            "X48hr.3hrDown":"X48hr.3hr",
            "X48hr.CCUp":"X48hr.CC",
            "X48hr.CCDown":"X48hr.CC",
            "H1.48hrUp":"H1.48hr",
            "H1.48hrDown":"H1.48hr",
            "H1.16hrUp":"H1.16hr",
            "H1.16hrDown":"H1.16hr",
            "H1.3hrUp":"H1.3hr",
            "H1.3hrDown":"H1.3hr",
            "H1.CCUp":"H1.CC",
            "H1.CCDown":"H1.CC",
            "H1.HkUp":"H1.Hk",
            "H1.HkDown":"H1.Hk",
            "X48hr.HkUp":"X48hr.Hk",
            "X48hr.HkDown":"X48hr.Hk",
            "X16hr.HkUp":"X16hr.Hk",
            "X16hr.HkDown":"X16hr.Hk",
            "X3hr.HkUp":"X3hr.Hk",
            "X3hr.HkDown":"X3hr.Hk",
            "H1.M5Up":"H1.M5",
            "H1.M5Down":"H1.M5",
            "X48hr.M5Up":"X48hr.M5",
            "X48hr.M5Down":"X48hr.M5",
            "X16hr.M5Up":"X16hr.M5",
            "X16hr.M5Down":"X16hr.M5",
            "X3hr.M5Up":"X3hr.M5",
            "X3hr.M5Down":"X3hr.M5"}

confidence_dict=dict() 
foldchange_dict=dict() 
confidence_header=('\t'+confidence[0]).split('\t') 
foldchange_header=('\t'+foldchange[0]).split('\t') 
for line in confidence[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[1] 
    startpos=tokens[2] 
    endpos=tokens[3] 
    key=(chrom,startpos,endpos) 
    confidence_dict[key]=dict() 
    for i in range(4,len(confidence_header)): 
        conf_val=float(tokens[i]) 
        confidence_dict[key][confidence_header[i]]=float(conf_val) 
print "built confidence dict" 
for line in foldchange[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[1] 
    startpos=tokens[2] 
    endpos=tokens[3] 
    key=(chrom,startpos,endpos) 
    foldchange_dict[key]=dict() 
    for i in range(4,len(foldchange_header)): 
        foldchange_val=float(tokens[i]) 
        foldchange_dict[key][foldchange_header[i]]=float(foldchange_val) 
print "built fold change dict" 
    

adjusted_dict=dict() 
peak_header=peak_entries[0].split('\t') 
for line in peak_entries[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[0] 
    startpos=tokens[1] 
    endpos=tokens[2] 
    key=(chrom,startpos,endpos) 
    adjusted_dict[key]=dict() 
    for i in range(3,len(peak_header)):
        cur_value=int(tokens[i]) 
        if peak_header[i] not in header_map: 
            #just pass the current value! 
            adjusted_dict[key][peak_header[i]]=cur_value 
        elif cur_value==1: 
            adjusted_dict[key][peak_header[i]]=cur_value 
        else: 
            cur_header=header_map[peak_header[i]]
            corresponding_confidence=confidence_dict[key][cur_header] 
            cur_value=int(tokens[i]) 
            fold_change_value=foldchange_dict[key][cur_header] 
            if corresponding_confidence > medium_thresh: 
                #0! 
                adjusted_dict[key][peak_header[i]]=0 
            elif corresponding_confidence > low_thresh: 
                #-1,if the fold change agrees with the column header 
                if (fold_change_value>0) and (peak_header[i].endswith('Up')): 
                    adjusted_dict[key][peak_header[i]]=-1 
                elif (fold_change_value<0) and (peak_header[i].endswith('Down')): 
                    adjusted_dict[key][peak_header[i]]=-1 
                else: 
                    adjusted_dict[key][peak_header[i]]=0 
            else: 
                #1, if the fold change agrees with the column header 
                if (fold_change_value>0) and (peak_header[i].endswith('Up')): 
                    adjusted_dict[key][peak_header[i]]=1 
                elif (fold_change_value<0) and (peak_header[i].endswith('Down')): 
                    adjusted_dict[key][peak_header[i]]=1 
                else:
                    adjusted_dict[key][peak_header[i]]=0 
#outf=open('thresholded_labels.txt','w')
outf=open('diffMatAnna.binarized.combined.thresholded.stringent.tsv','w')
outf.write('Chrom\tStartPos\tEndPos\tPeak\t'+'\t'.join(peak_header[3::])+'\n')
peak_id=0 
for entry in adjusted_dict: 
    outf.write('\t'.join(entry)+'\t'+'peak_'+str(peak_id))
    peak_id+=1 
    for task in peak_header[3::]: 
        label=adjusted_dict[entry][task] 
        outf.write('\t'+str(label))
    outf.write('\n') 

        


    
