#!/bin/bash

# Script to run DeepBind model on given intervals and a FASTA file

# Usage: bash deepBind_SOX9.sh intervals.bed fasta_file.fa

# Required input arguments:
#   - intervals.bed: BED file containing genomic intervals for prediction
#   - fasta_file.fa: FASTA file containing genomic sequences corresponding to the intervals

# DeepBind article: https://www.nature.com/articles/nbt.3300
# Requirement: https://github.com/kipoi/kipoi

# Extract command-line arguments
intervals=$1  # BED file containing genomic intervals
fasta_file=$2  # FASTA file containing genomic sequences

# Name of the pre-trained DeepBind model to use
# MODEL="DeepBind/Homo_sapiens/TF/D00649.002_SELEX_SOX9" #SOX9 mode
MODEL="DeepBind/Homo_sapiens/TF/D00411.003_SELEX_GATA4" #GATA4 model
# MODEL="DeepBind/Homo_sapiens/TF/D00296.006_SELEX_AR" #santy check model
# Extract the model name from the MODEL string
MODEL_NAME="${MODEL##*_}"

# Print the FASTA file name
echo "FASTA File: $fasta_file"

# Extract the base name of the intervals BED file (remove the '.bed' extension)
FASTA_NAME=`basename $intervals .bed`

# Define the output file name
output="${MODEL_NAME}_${FASTA_NAME}.tsv"

# Obtain an example dataset for the model
kipoi get-example $MODEL -o example

# Predict using the DeepBind model
# Arguments:
#   - intervals_file: BED file containing genomic intervals
#   - fasta_file: FASTA file containing genomic sequences
# Output:
#   - tmp_out.tsv: Temporary output file
kipoi predict $MODEL --dataloader_args="{\"intervals_file\": \"$intervals\", \"fasta_file\": \"$fasta_file\"}" -o "tmp_out.tsv"

# Remove the header (first line) from the temporary output file
sed -i '1d' tmp_out.tsv

# Combine the original intervals BED file and the predicted scores
# The resulting file will have columns: [original_interval_info, predicted_scores]
paste $intervals  tmp_out.tsv  > $output

# Remove the temporary output file
rm -f tmp_out

# Display the first few lines of the output file
head  $output

# Print the location of the final output file
echo "Output can be found at: $output"
