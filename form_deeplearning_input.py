import random
random.seed(1234) 
from fastaFromBed import * 
from generateSplits import *
import argparse
import pdb
from Params import * 

def parse_args():
    parser=argparse.ArgumentParser(description="Provide a file containing the full path name to bed files for all tasks in formatin task_name\tfile_name")
    parser.add_argument("--task_bed",help="full path to file containing the full path name to bed files for all tasks in format task_name\tfile_name")
    parser.add_argument("--out_prefix",help="output path to where fasta,train/test/valid, and labels should be stored")
    parser.add_argument("--browser_bed_dir",help="directory where to store per-task bed files for qc in the browser")
    parser.add_argument("--slide_by",type=int,help="bin sizes to split the genome into",default=250)
    parser.add_argument("--seq_size",type=int,help="length of final sequence to use in training",default=2000)
    parser.add_argument("--sample",type=float,help="fraction of peaks to use in dataset",default=1.0)
    parser.add_argument("--universal_negative_source",type=str,help="bed file of universal negative regions, optional",required=False)
    parser.add_argument("--universal_negative_fraction",type=float,help="fraction of dataset that should consist of universal negative datapoints",default=0.05,required=False)
    parser.add_argument("--low_confidence_task_bed",help="bed file (i.e. from naive overlap) to intersect w/ task_bed argument to determine peaks of low confidence",required=False) 
    return parser.parse_args()


def add_entry(data_windows,entry,val):
    #make sure we're not running off the edge of the chromosome 
    if entry[1]<1:
        return data_windows
    if val ==1:
        #override any existing entry with the current entry if the value is 1 
        data_windows[entry]=val
    #only add a 0 entry if an entry does not currently exist in the dictionary! 
    elif entry not in data_windows:
        data_windows[entry]=val
    return data_windows

def add_universal_negatives(observed_bins,data_windows,universal_negative_source,universal_negative_fraction,slide_by,seq_size):
    flank_size=(seq_size-slide_by)/2
    num_to_add=int(round(len(observed_bins)*universal_negative_fraction))
    neg_set=open(universal_negative_source,'r').read().strip().split('\n')
    print("loaded universal negative set") 
    size_neg_set=len(neg_set) 
    num_added=0
    #we want to generate a dictionary of observed bins to avoid having to do a slow set search
    observed_bin_dict=dict()
    for entry in observed_bins:
        observed_bin_dict[entry]=1
    
    while num_added < num_to_add:
        #select a random peak from the negative set
        line_index=random.randint(0,size_neg_set-1)
        line=neg_set[line_index]
        tokens=line.split('\t')
        chrom=tokens[0]
        start_pos=int(tokens[1])
        end_pos=int(tokens[2]) 
        first_bin_start=slide_by*(start_pos/slide_by)
        last_bin_start=slide_by*(end_pos/slide_by)
        first_bin_label=1
        last_bin_label=1
        if (start_pos >=first_bin_start) & (end_pos <=first_bin_start+slide_by):
            fullyInBin=1
        else:
            fullyInBin=0
        if not fullyInBin:
            if(start_pos - first_bin_start)>(slide_by/2):
                first_bin_label=0
            if(end_pos - last_bin_start) < (slide_by/2):
                last_bin_label=0
        if (first_bin_label==1):
            entry=tuple([chrom,first_bin_start-flank_size,first_bin_start-flank_size+seq_size])
            if ((entry[1]>0) and (entry not in observed_bin_dict)):
                observed_bins.add(entry)
                observed_bin_dict[entry]=1 
                num_added+=1 
                #add as universal negative
                for task in data_windows:
                    data_windows[task][entry]=0
        if (last_bin_label==1):
            entry=tuple([chrom,last_bin_start-flank_size,last_bin_start-flank_size+seq_size])
            if ((entry[1]>0) and (entry not in observed_bin_dict)):
                observed_bin_dict[entry]=1
                observed_bins.add(entry) 
                num_added+=1
                #add as a universal negative
                for task in data_windows:
                    data_windows[task][entry]=0
        for bin_index in range(first_bin_start+slide_by,last_bin_start,slide_by):
            entry=tuple([chrom,bin_index-flank_size,bin_index-flank_size+seq_size])
            if ((entry[1]>0) and ( entry not in observed_bin_dict)):
                observed_bin_dict[entry]=1
                observed_bins.add(entry) 
                num_added+=1
                #add as a universal negative!
                for task in data_windows:
                    data_windows[task][entry]=0             
    return observed_bins,data_windows



def generate_data_windows(task_name,task_file,slide_by,seq_size,sample_rate,observed_bins):
    data=open(task_file,'r').read().strip().split('\n')
    data_windows=dict()
    #get the flank size
    flank_size=(seq_size-slide_by)/2
    for line in data:
        #if sampling is used, generate a random value to determine whether to include this datapoint in the dataset
        if sample_rate < 1.0:
            randval=random.random()
            if randval > sample_rate:
                continue 
        tokens=line.split('\t')[0:3]
        chrom=tokens[0]
        start_pos=int(tokens[1])
        end_pos=int(tokens[2])
        
        #get the first and last bins
        first_bin_start=slide_by*(start_pos/slide_by)
        last_bin_start=slide_by*(end_pos/slide_by)
        #check if >=1/2 the bin is overlapping with the peak 
        first_bin_label=1
        last_bin_label=1
        #check if the peak is fully contained in the bin -- correctly labels 1 for very sharp peaks 
        if (start_pos >= first_bin_start) & (end_pos <=first_bin_start + slide_by):
            fullyInBin=1
        else:
            fullyInBin=0
        if not fullyInBin: 
            if (start_pos - first_bin_start)>(slide_by/2):
                first_bin_label=0
            if (end_pos - last_bin_start) < (slide_by/2):
                last_bin_label=0
                
        entry=tuple([chrom,first_bin_start-flank_size,first_bin_start-flank_size+seq_size])
        data_windows=add_entry(data_windows,entry,first_bin_label)
        entry=tuple([chrom,last_bin_start - flank_size,last_bin_start-flank_size+seq_size])
        data_windows=add_entry(data_windows,entry,last_bin_label)
        
        #get the left & right negative flanks 
        left_negative_flank_start=first_bin_start-slide_by-flank_size 
        left_negative_flank_end=left_negative_flank_start+seq_size 
        right_negative_flank_start=last_bin_start+slide_by-flank_size 
        right_negative_flank_end=right_negative_flank_start+seq_size
        entry=tuple([chrom,left_negative_flank_start,left_negative_flank_end])
        data_windows=add_entry(data_windows,entry,0)
        entry=tuple([chrom,right_negative_flank_start,right_negative_flank_end])
        data_windows=add_entry(data_windows,entry,0)

        #fill in any bins between the first & last bins with a positive label
        for bin_index in range(first_bin_start+slide_by,last_bin_start,slide_by):
            entry=tuple([chrom,bin_index-flank_size,bin_index-flank_size+seq_size])
            data_windows=add_entry(data_windows,entry,1)
    return data_windows

#writes file of labels for each task
def write_label_file(data_dict,data_bins,out_prefix):
    outf=open(out_prefix+'.labels.txt','w')
    tasks=list(data_dict.keys())
    header="Peak"+'\t'+'\t'.join(tasks)
    outf.write(header+'\n')
    for data_bin in data_bins:
        to_write=['_'.join([str(i) for i in data_bin])]
        allzeros=True 
        for task in tasks:
            if data_bin in data_dict[task]:
                to_write.append(str(data_dict[task][data_bin]))
            else:
                to_write.append("0")
        outf.write('\t'.join(to_write)+'\n')
        
        
#test if a peak in the 'test_chroms' or 'validation_chroms', as defined in Params.py 
def split_folds(data_bins,out_prefix):
    outf_train=open(out_prefix+'.train.txt','w')
    outf_test=open(out_prefix+'.test.txt','w')
    outf_valid=open(out_prefix+'.validate.txt','w')
    for data_bin in data_bins:
        peak_name='_'.join([str(i) for i in data_bin])
        chrom=data_bin[0]
        if chrom in test_chroms:
            outf_test.write(peak_name+'\n')
        elif chrom in validation_chroms:
            outf_valid.write(peak_name+'\n')
        else:
            outf_train.write(peak_name+'\n')
        


#generates a browser track of the positive and negative bins 
def dump_bed(task_name,task_window_dict,browser_bed_dir,window_size,full_size):
    outf_pos=open(browser_bed_dir+'/'+task_name+".positives.bed",'w')
    outf_neg=open(browser_bed_dir+'/'+task_name+".negatives.bed",'w')
    counter_pos=0
    counter_neg=0
    delta=(full_size-window_size)/2
    for key in sorted(task_window_dict.iterkeys()):
        startpos=key[1]+delta
        endpos=key[2]-delta
        chrom=key[0]
        entry=chrom+'\t'+str(startpos)+'\t'+str(endpos)
        #print str(entry)
        #negative peak
        if task_window_dict[key]==0:
            outf_neg.write(entry+'\t.\t'+str(counter_neg)+'\t.\n')
            counter_neg+=1
        else:
            #positive peak
            outf_pos.write(entry+'\t.\t'+str(counter_pos)+'\t.\n')
            counter_pos+=1
         
def main():
    args=parse_args() 
    task_bed_file=open(args.task_bed,'r').read().strip().split('\n')
    browser_bed_dir=args.browser_bed_dir 
    data_windows=dict()
    observed_bins=set()  
    
    for line in task_bed_file:
        tokens=line.split('\t')
        task_name=tokens[0]
        task_file=tokens[1]
        print('generating labels for task:'+str(task_name))
        data_windows[task_name]=generate_data_windows(task_name,task_file,args.slide_by,args.seq_size,args.sample,observed_bins)
        #write an output bed file for the task
        dump_bed(task_name,data_windows[task_name],browser_bed_dir,args.slide_by,args.seq_size)
        #update the full set of observed bins        
        observed_bins=observed_bins.union(data_windows[task_name].keys())
        
    #generate a set of universal negatives
    if args.universal_negative_source!=None:
        print("adding universal negative bins!")
        observed_bins,data_windows=add_universal_negatives(observed_bins,data_windows,args.universal_negative_source,args.universal_negative_fraction,args.slide_by,args.seq_size)
       
    #write a labels file
    print('writing labels to file')
    #we apply some sort of filtering, such as removing peaks that are universal negatives 
    write_label_file(data_windows,observed_bins,args.out_prefix)
    
    #generate a train/test/valid split
    print('generating train/test/validation splits') 
    split_folds(observed_bins,args.out_prefix)
    
    #generate fasta file
    print('extracting fasta sequence') 
    make_fasta(refseq,observed_bins,args.out_prefix)
    
if __name__=="__main__": 
    main() 
