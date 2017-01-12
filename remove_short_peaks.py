import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="Provide a bed file and a threshold to remove very narrow Peaks")
    parser.add_argument("--i",help="input bed file")
    parser.add_argument("--o",help="output filtered file")
    parser.add_argument("--thresh",help="minimum allowed peak width")
    return parser.parse_args()

def main():
    args=parse_args() 
    outf=open(args.o,'w')
    thresh=int(args.thresh) 
    data=open(args.i,'r').read().strip().split('\n')
    val_dict=dict() 
    for line in data:
        tokens=line.split('\t')
        chrom=tokens[0]
        default_start=int(tokens[1])
        default_end=int(tokens[2])
        width=default_end-default_start
        if width not in val_dict:
            val_dict[width]=1
        else:
            val_dict[width]+=1 
        if (default_end-default_start) < thresh:
            continue
        else:
            outf.write(line+'\n')
    outf=open(args.o+'.dict','w')
    keys=val_dict.keys()
    keys.sort()
    for key in keys:
        outf.write(str(key)+'\t'+str(val_dict[key])+'\n')
if __name__=="__main__":
    main() 
