#output_data_dir="/srv/scratch/annashch/stemcells/het/anna_data/" 
#output_data_dir="/srv/scratch/annashch/stemcells/het/anna_all_timepoint_pairs/" 
output_data_dir="/srv/scratch/annashch/stemcells/het/anna_code/" 
input_data_dir="/srv/scratch/annashch/stemcells/het/peaks/peaks/" 
tmp_dir='/srv/scratch/annashch/stemcells/het/tmp/' #directory to store temp files generated by samtools etc. 

#samples for the study 
sample_dict=dict() 
sample_dict['CC']=['10-CC_rep2','9-CC_rep1'] 
sample_dict['H1']=['18-H1_rep1','19-H1_rep2','20-H1_rep3']
sample_dict['48hr']=['21-48hr_rep1','7-48hr_rep2','8-48hr_rep3']
sample_dict['3hr']=['4-3hr_rep1','5-3hr_rep2','6-3hr_rep3']
sample_dict['16hr']=['1-16hr_rep1','2-16hr_rep2','3-16hr_rep3'] 
sample_dict['ES']=['11-ES_rep1','12-ES_rep2'] 
sample_dict['Hk']=['13-Hk_rep1','14-Hk_rep2']
sample_dict['M5']=['22-M5_rep1','15-M5_rep2','16-M5_rep3','17-M5_rep4'] 

#remove any peaks with -*log(q-value) < threshold 
idr_threshold=2 


#pairs (stores adjacent timepoints to merge together) 
pairs=[] 
pairs.append(tuple(['CC','3hr']))
pairs.append(tuple(['3hr','16hr']))
pairs.append(tuple(['16hr','48hr']))
pairs.append(tuple(['48hr','H1']))

# number of aligned reads from each replicate BAM file 
# for normalization 
scale_factor=15000000.0
aligned_reads=dict() 
aligned_reads['10-CC_rep2']=10718751
aligned_reads['9-CC_rep1']=12070949
aligned_reads['19-H1_rep2']=5661355
aligned_reads['20-H1_rep3']=8894062
aligned_reads['21-48hr_rep1']=3933662
aligned_reads['7-48hr_rep2']=14410904
aligned_reads['8-48hr_rep3']=11264723
aligned_reads['4-3hr_rep1']=6213189
aligned_reads['5-3hr_rep2']=9207572
aligned_reads['6-3hr_rep3']=6343355
aligned_reads['1-16hr_rep1']=6013770
aligned_reads['2-16hr_rep2']=3647924
aligned_reads['3-16hr_rep3']=12459317
aligned_reads['11-ES_rep1']=279184
aligned_reads['12-ES_rep2']=2694
aligned_reads['13-Hk_rep1']=21290440
aligned_reads['14-Hk_rep2']=11425647
aligned_reads['17-M5_rep4']=14053583
aligned_reads['16-M5_rep3']=43636805
aligned_reads['22-M5_rep1']=27345710
aligned_reads['15-M5_rep2']=42614152
aligned_reads['18-H1_rep1']=7765690

#parameters for sliding window 
windowSize=2000 
windowOverlap=2000
line_length=2000 
jitterMax=15 
jitterProb=0.05 

#parameters for data splitting/ FASTA extraction 
external_data_dir = '/srv/scratch/annashch/deeplearning/heterokaryon/data/'
deeplearning_inputs_dir = '/srv/scratch/annashch/deeplearning/heterokaryon/inputs/'
pluripotency_region_file = external_data_dir + 'pluripotency_gene_regions.bed'
refseq = 'hg19.genome.fa'


#out_fasta = deeplearning_inputs_dir + 'specialized.deeplift.fasta'
#out_splits_train = deeplearning_inputs_dir + 'specialized.deeplift.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'specialized.deeplift.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'specialized.deeplift.validate.txt'
#out_labels = deeplearning_inputs_dir + 'specialized.deeplift.labels.txt'


out_fasta = deeplearning_inputs_dir + 'regression.fasta'
out_splits_train = deeplearning_inputs_dir + 'regression.train.txt'
out_splits_test = deeplearning_inputs_dir + 'regression.test.txt'
out_splits_validate = deeplearning_inputs_dir + 'regression.validate.txt'
out_labels = deeplearning_inputs_dir + 'regression.labels.txt'

#out_fasta = deeplearning_inputs_dir + 'remapped.fasta'
#out_splits_train = deeplearning_inputs_dir + 'remapped.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'remapped.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'remapped.validate.txt'
#out_labels = deeplearning_inputs_dir + 'remapped.labels.txt'



#out_fasta = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.win50.full.fasta'
#out_splits_train = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.win50.full.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.win50.full.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.win50.full.validate.txt'
#out_labels = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.win50.full.labels.txt'


#out_fasta = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.fasta'
#out_splits_train = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.validate.txt'
#out_labels = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.labels.txt'

#out_fasta = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.fasta'
#out_splits_train = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.validate.txt'
#out_labels = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.idr.sliced.labels.txt'

#out_fasta = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.nb.sliced.fasta'
#out_splits_train = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.nb.sliced.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.nb.sliced.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.nb.sliced.validate.txt'
#out_labels = deeplearning_inputs_dir + 'heterokaryon.sliding.pangwei.nb.sliced.labels.txt'

#out_fasta = deeplearning_inputs_dir + 'heterokaryon.sliding.fasta'
#out_splits_train = deeplearning_inputs_dir + 'heterokaryon.sliding.train.txt'
#out_splits_test = deeplearning_inputs_dir + 'heterokaryon.sliding.test.txt'
#out_splits_validate = deeplearning_inputs_dir + 'heterokaryon.sliding.validate.txt'
#out_labels = deeplearning_inputs_dir + 'heterokaryon.sliding.labels.txt'
pluripotent_region_flank = 10000
test_chroms = ['chr2', 'chr3']
validation_chroms = ['chr1', 'chr4']

#augment the negative peaks, otherwise augment the positive peaks 
jitter_minus=[3,4,5,6] 