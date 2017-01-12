import argparse
def parse_args():
    parser=argparse.ArgumentParser(description="Provide a gappedPeak file")
    parser.add_argument("--i",help="input gapped peak file")
    parser.add_argument("--o",help="output filtered file")
    return parser.parse_args()

def main():
    args=parse_args() 
    outf=open(args.o,'w')
    data=open(args.i,'r').read().strip().split('\n')
    for line in data:
        tokens=line.split('\t')
        chrom=tokens[0]
        default_start=int(tokens[1])
        default_end=int(tokens[2])
        thick_starts=[int(i) for i in tokens[11].split(',')]
        thick_sizes=[int(i) for i in tokens[10].split(',') ]
        #no exons present 
        if max(thick_sizes)==1:
            continue
        for i in range(len(thick_starts)):
            adjusted_start=default_start+thick_starts[i]
            adjusted_end=adjusted_start+thick_sizes[i]
            outf.write(chrom+'\t'+str(adjusted_start)+'\t'+str(adjusted_end)+'\n')
            
if __name__=="__main__":
    main() 
