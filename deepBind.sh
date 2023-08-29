
 
#deepBind article: https://www.nature.com/articles/nbt.3300
# requierment: https://github.com/kipoi/kipoi



intervals=$1
fasta_file=$2
# conda activate kipoi-DeepBind
MODEL="DeepBind/Homo_sapiens/TF/D00649.002_SELEX_SOX9"
MODEL_NAME="${MODEL##*_}"
echo $faste_file
FASTA_NAME=`basename $fasta_file .fa`
output="${MODEL_NAME}_${FASTA_NAME}.tsv" 

kipoi get-example $MODEL -o example



kipoi predict $MODEL --dataloader_args="{\"intervals_file\": $intervals, \"fasta_file\": $fasta_file}" -o "tmp_out.tsv"
# check the results

sed -i '1d' tmp_out.tsv 
 paste $intervals  tmp_out.tsv  > $output
rm -f tmp_out
head  $output 
echo "output can be found at- $output"

