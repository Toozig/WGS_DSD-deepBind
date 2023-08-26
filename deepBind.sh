
 
#deepBind article: https://www.nature.com/articles/nbt.3300
# requierment: https://github.com/kipoi/kipoi


module load hurcs homer # needed for the refrence file 
REFRENCE='/usr/local/hurcs/homer/4.11/data/genomes/hg38/genome.fa'
conda activate kipoi-DeepBind
MODEL="DeepBind/Homo_sapiens/TF/D00649.002_SELEX_SOX9"
INTERVALS="example/intervals_file"


OUTPUT=""


kipoi get-example $MODEL -o example



kipoi predict $MODEL --dataloader_args="{\"intervals_file\": $INTERVALS, \"fasta_file\": $REFRENCE}" -o $OUTPUT
# check the results
head '/tmp/DeepBind|Homo_sapiens|TF|D00649.002_SELEX_SOX9.example_pred.tsv'

