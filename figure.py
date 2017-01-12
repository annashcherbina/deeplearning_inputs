header=open('header','r').read()
outf=open('figure.SLICED.tsv','w') 
outf.write(header+'\n') 
outf.write('chr6\t31141273\t31140477\t'+'\t'.join(['1']*43)+'\n')
outf.write('chr1\t59239000\t59240000\t'+'\t'.join(['1']*43)+'\n')
outf.write('chr1\t26713000\t26755000\t'+'\t'.join(['1']*43)+'\n')
outf.write('chr19\t45971252\t45978437\t'+'\t'.join(['1']*43)+'\n')
outf.write('chr1\t154377668\t154441926\t'+'\t'.join(['1']*43)+'\n')
outf.write('chr8\t23536206\t23540402\t'+'\t'.join(['1']*43)+'\n')


