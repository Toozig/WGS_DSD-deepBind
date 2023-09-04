import numpy as np
import pandas as pd
import seaborn as sns
import bisect
from bin.deepBind_functions import *
from seaborn import color_palette
from dna_features_viewer import GraphicFeature, GraphicRecord


# COLINDEX
START = 1
END = 2
PRED_SCORE = 3
# COL NAMES
FLOAT_COLS = ['START','END','PREDICTION_SCORE']
PRED_COLS = ['SEQ_ID'] + FLOAT_COLS




def point_to_segment(segment_df: pd.DataFrame, point_arr: np.ndarray) -> np.ndarray:
    """
    Maps points to segments in a continuous range.

    Parameters:
    segment_df (pd.DataFrame): DataFrame with 'START' and 'END' columns representing continuous, non-overlapping ranges.
    point_arr (np.ndarray): Array of points to map to segments.

    Returns:
    np.ndarray: Array with segment indices for each point.
    """
    flat_ranges = segment_df.flatten()
    intersection = (point_arr.apply(lambda x: bisect.bisect_left(flat_ranges, x)) // 2)
    return intersection


def process_prediction_file(file_path: str, seq_id: str) -> pd.DataFrame:
    """
    use only for predictions where each segment is overlapping (only) with it's previous segment
    Process a prediction file -
        1) open to pandas DF
        2) adds cols names
        3) create a DF with each overlapping segments on the same rows


    Parameters:
    file_path (str): Path to the prediction file in TSV format.

    Returns:
    pd.DataFrame: Processed DataFrame containing relevant data.

    Example:
    >>> processed_df = process_prediction_file('SOX9_en13_seq_win_16_shift_8.tsv')
    >>> print(processed_df.head())
    """
    pred_df = pd.read_csv(file_path, sep='\t', header=None)
    pred_df.columns = PRED_COLS
    pred_df.loc[:, FLOAT_COLS] = pred_df.loc[:, FLOAT_COLS].astype(float)
    pred_df = pred_df[pred_df.SEQ_ID == seq_id]
    
    shifted = pred_df.shift(-1)
    shifted.columns = shifted.columns + '_shift'
    
    shift_df = pd.concat([pred_df, shifted], axis=1).dropna()
    
    sanity = (shift_df.START_shift > shift_df.START) & (shift_df.START_shift < shift_df.END)
    # print(sanity.all())
    return shift_df.reset_index(drop=True)

def prepare_range_df(file_path: str,  seq_id: str, seq_length: int) -> pd.DataFrame:
    """
    Prepare a range DataFrame from the output of the previous function.
    returns a nX2 df with START & END for each segment and it's mean score (each segments have 2 prediction score)
    

    Parameters:
    shift_df (pd.DataFrame): Processed DataFrame from the process_prediction_file.
    seq_length (int): Length of the sequence.

    Returns:
    pd.DataFrame: Range DataFrame with 'START', 'END', and 'SCORE' columns.

    Example:
    >>> range_df = prepare_range_df(processed_df, sequence_length)
    >>> print(range_df.head())
    """
    shift_df = process_prediction_file(file_path,seq_id)

    # Prepare the first and last lines
    first_line = shift_df.loc[0, ['START', 'START_shift', 'PREDICTION_SCORE']].rename({'START_shift': 'END', 'PREDICTION_SCORE': 'SCORE'})

    first_line['END'] = first_line['END'] - 1
    last_line = shift_df.loc[shift_df.shape[0] - 1, ['END', 'END_shift', 'PREDICTION_SCORE_shift']].rename({'END': 'START',
     'END_shift': 'END', 'PREDICTION_SCORE_shift': 'SCORE'})
    first_line['START'] = first_line['START'] + 1
 
    # Create new columns for the range DataFrame
    new_start = shift_df.START_shift
    new_end = shift_df.END
    new_score = (shift_df.PREDICTION_SCORE + shift_df.PREDICTION_SCORE_shift) / 2
    
    # Create the main portion of the range DataFrame
    new_df = pd.DataFrame({'START': new_start, 'END': new_end, 'SCORE': new_score})
    new_df['END'] = new_df['END'] - 1
    
    # Prepare lines for the start and end of the sequence
    zero_line_s = pd.Series([0, first_line.START - 1, 0], index=first_line.index)
    zero_line_e = pd.Series([last_line.END + 1, seq_length - 1, 0], index=first_line.index)
    
    # Concatenate all lines to create the final range DataFrame
    range_df = pd.concat([pd.DataFrame(zero_line_s).T, pd.DataFrame(first_line).T, new_df, pd.DataFrame(last_line).T, pd.DataFrame(zero_line_e).T]).reset_index(drop=True)
    
    return range_df


def plot_get_bs_feature(seq_df, bs_col, bs_dict, cmap='pastel'):
    """
    Generate DNA sequence binding site features for visualization.

    This function processes a sequence DataFrame, a binding site column,
    to create DNA sequence binding site features suitable for visualization using the
    `dna_features_viewer` library.

    Parameters:
    seq_df (pd.DataFrame): Sequence DataFrame.
    bs_col (str): Column containing binding site label.

    Returns:
    list: List of GraphicFeature objects representing binding site features.

    Example:
    >>> colors = ["#FFB3BA", "#FFFFBA", "#Baffc9"]
    >>> bs_features = plot_get_bs_feature(sequence_df, "Binding_Sites", colors)
    """
    colors = sns.color_palette(cmap)

    features = []
    bs_list = seq_df[bs_col].unique()
    bs_list = bs_list[bs_list != '']

    for i in range(len(bs_list)):
        bs = bs_list[i]
        bs_df = seq_df[seq_df[bs_col] == bs]
        seq = bs_dict[bs]

        while bs_df.shape[0]:
            start = bs_df.index.min()
            end = start + len(seq) - 1
            features.append(GraphicFeature(start=start, end=end,
                                       strand=+1, color=colors[i], label=bs))
            bs_df = bs_df.drop(range(start,end + 1))
    
    return features



def add_binding_site_label(seq_df, binding_site, binding_site_name, binding_site_col):
    """
    Label a binding site in a sequence DataFrame.

    Parameters:
    seq_df (pd.DataFrame): Sequence DataFrame.
    binding_site (str): Binding site to label.
    binding_site_name (str): Label for the binding site.
    binding_site_col (str): Column for the label.

    Returns:
    pd.DataFrame: Sequence DataFrame with added label.

    Example:
    sequence_df = add_binding_site_label(sequence_df, "ACTG", "Promoter", "Region_Label")
    """
    cur_seq = "".join(seq_df.base)
    s_idx = cur_seq.find(binding_site)
    save = 0
    while cur_seq:
        s_idx = cur_seq.find(binding_site)
        if s_idx == -1:
            break
        edge = s_idx + save + len(binding_site) - 1
        seq_df.loc[s_idx + save : edge , binding_site_col] = binding_site_name
        save = edge + 1
        cur_seq = cur_seq[s_idx + len(binding_site) - 1:]

    return seq_df


def deepBind_score_df(seq_id, seq_dict, score_file):
    """
    Retrieve deepBind scores and annotate a DNA sequence.

    This function takes a sequence ID, a dictionary of sequences, and a score file path.
    It calculates deepBind scores for the sequence and annotates the sequence DataFrame
    with these scores.

    Parameters:
    seq_id (str): Identifier for the sequence in seq_dict.
    seq_dict (dict): Dictionary containing sequences.
    score_file (str): Path to the score file.

    Returns:
    pd.DataFrame: Annotated sequence DataFrame with deepBind scores.

    Example:
    >>> seq_dict = {"seq1": "ACTGAGCTAG", "seq2": "CGTAGCTA"}
    >>> score_file = "scores.csv"
    >>> annotated_df = deepBind_score_df("seq1", seq_dict, score_file)
    """

    # Retrieve the sequence based on seq_id
    seq = seq_dict[seq_id]
    
    # Prepare a DataFrame for the sequence
    segment_df = prepare_range_df(score_file, seq_id, len(seq))
    seq_df = pd.DataFrame(list(seq)).reset_index().rename(columns={"index": 'POS', 0: 'base'})
    
    # Calculate deepBind scores and annotate the sequence DataFrame
    score_idx = point_to_segment(segment_df[['START', 'END']].to_numpy(), seq_df.POS)
    seq_df['SCORE'] = segment_df.SCORE[score_idx].tolist()
    
    return seq_df

def get_binding_site_label_df(seq_df : pd.DataFrame, bs_dict : dict, label_col_name: str) -> pd.DataFrame:
    result = seq_df[['POS', 'base']].copy()
    result[label_col_name] = ''
    for b in bs_dict.keys():
        result = add_binding_site_label(result,bs_dict[b],b,label_col_name)
    return result



import numpy as np
import matplotlib.pyplot as plt

def color_rows_in_matrix(mat, labels, cmap_name='tab10', fig_width=10, fig_height=6):
    num_rows, num_cols = mat.shape
    cmap = plt.get_cmap(cmap_name, num_rows)

    colored_matrix = np.zeros((num_rows, num_cols, 4))  # Include alpha (transparency) channel
    
    uniq_lables = list(set([i.split('_')[0] for i in labels]))
    color_dict = {uniq_lables[i]: cmap(i) for i in range(len(uniq_lables)) }

    for i in range(num_rows):
        row_color = color_dict[labels[i].split('_')[0]]
        alpha_values = mat[i]  # Use matrix values for alpha (opacity)
        
        for j in range(num_cols):
            colored_matrix[i, j, :3] = row_color[:3]  # Assign RGB color
            colored_matrix[i, j, 3] = alpha_values[j]  # Assign alpha value from the matrix

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.imshow(colored_matrix, aspect=10)
    ax.axis('off')
    
    if labels is not None:
        handles = [plt.Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor=cmap(i)) for i, label in enumerate(uniq_lables)]
        legend = ax.legend(handles=handles, title="Labels", bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.add_artist(legend)  # Add the legend back to the plot

    plt.show()


