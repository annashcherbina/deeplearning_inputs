source=open('diffMatAnna.binarized.combined.thresholded.stringent.tsv','r').read().strip().split('\n') 
outf=open('debug.deeplift.input.tsv','w') 
outf.write(source[0]+'\n') 
linelength=2000 
index=0 
for line in source[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[0] 
    startpos=int(tokens[1]) 
    endpos=int(tokens[2]) 
    summit=startpos+0.5*(float(endpos)-startpos)
    new_start=int(max([1,summit-0.5*linelength]))
    new_end=new_start+linelength 
    peak_id=index 
    outf.write(chrom+'\t'+str(new_start)+'\t'+str(new_end)+'\t'+str(peak_id)+'\t'+'\t'.join(tokens[3::])+'\n')
    index+=1
