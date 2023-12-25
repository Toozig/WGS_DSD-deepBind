import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
    
ID = 'seq_id'
SEQ_LEN = 'seq_length'

def print_usage():
    print("Usage:")
    print("python script_name.py FASTA_FILE WINDOW_SIZE SHIFT_SIZE")
    print("Example:")
    print("python my_script.py input.fasta 100 50")

def validate_input(args):
    if len(args) != 4:
        print("Error: Incorrect number of arguments.")
        print_usage()
        sys.exit(1)
    
    fasta_file = args[1]
    window_size = args[2]
    shift_size = args[3]
    
    try:
        window_size = int(window_size)
        shift_size = int(shift_size)
    except ValueError:
        print("Error: WINDOW_SIZE and SHIFT_SIZE must be integers.")
        print_usage()
        sys.exit(1)
    
    return fasta_file, window_size, shift_size


def get_sequence_lengths(fasta_filename):
    """
    Get the lengths of sequences from a FASTA file.

    Parameters:
        fasta_filename (str): Path to the input FASTA file.

    Returns:
        dict: A dictionary with sequence IDs as keys and sequence lengths as values.
    """
    sequence_lengths = {record.id: len(record.seq) for record in SeqIO.parse(fasta_filename, "fasta")}
    return sequence_lengths

def generate_bed_file(seq_id, seq_length, window_size, shift_size,boundry = 50):
    """
    Generate a DataFrame with BED file records based on sequence information.

    Parameters:
        seq_id (str): Sequence ID.
        seq_length (int): Length of the sequence.
        window_size (int): Window size for the BED segments.
        shift_size (int): Shift size for sliding the window.

    Returns:
        pd.DataFrame: A DataFrame containing BED records.
    """
    start = np.arange(start=boundry, stop=seq_length, step=shift_size)
    end = start + window_size
    
    seq_id_vec = np.full(end.shape[0], seq_id)
    bed = pd.DataFrame(data=zip(seq_id_vec, start, end), columns=["seq_id", "start", "end"])
    bed = bed[bed['end'] <= seq_length - boundry ]
    return bed

def main(fasta_filename, window_size, shift_size, output):
    """
    Generate BED file segments based on sequence lengths and specified parameters.

    Parameters:
        fasta_filename (str): Path to the input FASTA file.
        window_size (int): Window size for the BED segments.
        shift_size (int): Shift size for sliding the window.
        output (str): Path to the output BED file.
    """
    seq_dict = get_sequence_lengths(fasta_filename)
    sequence_lengths_df = pd.DataFrame(seq_dict.items(), columns=[ID, SEQ_LEN])
    bed_df = sequence_lengths_df.apply(lambda x: generate_bed_file(x[ID], x[SEQ_LEN], window_size, shift_size), axis=1)
    bed_df = pd.concat(bed_df.tolist(), ignore_index=True)
    bed_df.to_csv(output, sep='\t', index=False, header=None)
    return bed_df



if __name__ == "__main__":
    fasta_file, window_size, shift_size = validate_input(sys.argv)
    
    output_filename = f"{fasta_file.split('.')[0]}_win_{window_size}_shift_{shift_size}.bed"
    print("Input validated successfully.")
    print(f"FASTA File: {fasta_file}")
    print(f"Window Size: {window_size}")
    print(f"Shift Size: {shift_size}")
    print(f"Output Filename: {output_filename}")
    print("Prepering BED file")
    main(fasta_file, window_size, shift_size, output_filename)
    print(f"segment file generated successfully.")