#QC
#python form_deeplearning_input.py  --task_bed task_inputs_gecco.small.txt --out_prefix /srv/scratch/annashch/deeplearning/gecco/inputs/gecco.small --browser_bed_dir /srv/www/kundaje/annashch/gecco/model_windows --slide_by 500 --seq_size 2000

#FULL DATASET 
#python form_deeplearning_input.py  --task_bed task_inputs_gecco.txt --out_prefix /srv/scratch/annashch/deeplearning/gecco/inputs/gecco --browser_bed_dir /srv/www/kundaje/annashch/gecco/model_windows --slide_by 250 --seq_size 2000 --sample 0.10 --universal_negative_source /srv/scratch/annashch/deeplearning/gecco/encode_dnase/final_merged_set.bed --universal_negative_fraction 0.05


#SUBSAMPLE TO 10%
python form_deeplearning_input.py  --task_bed task_inputs_gecco.txt --out_prefix /srv/scratch/annashch/deeplearning/gecco/inputs/gecco.sampled --browser_bed_dir /srv/www/kundaje/annashch/gecco/model_windows_sampled --slide_by 250 --seq_size 2000 --sample 0.10 --universal_negative_source /srv/scratch/annashch/deeplearning/gecco/encode_dnase/final_merged_set.bed --universal_negative_fraction 0.05

