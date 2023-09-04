# Introduction

This Jupyter notebook will explore the use of the "deepBind" model on enhancer13. <br> I will present the data being used, the input and output of the model, as well as the graphical results it produces.

## Data - Enhancer13

Within a 5.2 kb genomic segment lies Enhancer 13, or eSR-A, a vital player in the orchestration of sexual development. It exclusively exerts its influence in human embryonic testes but remains dormant in ovaries. eSR-A's core, a 1514 bp fragment packed with SOX9 and SF1 binding sites, activates vigorously when exposed to SF1 and SOX9. <br>

This enhancer's significance transcends humans, with a parallel counterpart, Enh13, discovered in mice. In mice, Enh13 deletion leads to complete sex reversal in XY individuals.  <br>

The stakes are high for eSR-A, as its deletion in humans results in 46,XY sex reversal, and duplication leads to 46,XX (ovo) testicular disorders of sexual development. This discovery underscores the intricate regulation of genes like SOX9 in sexual development and highlights its conservation and critical role in human and mouse biology.

## segment of Enhancer 13:



    
![png](deepBind_analysis_files/deepBind_analysis_4_0.png)
    


## Model Overview

DeepBind is a neural network model used for predicting DNA- and RNA-binding proteins. <br>It was trained on a large dataset of protein sequences and their corresponding binding information. <br>The model aims to accurately identify regions of proteins that interact with DNA or RNA molecules.

Read more about the model here - [Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning](https://www.nature.com/articles/nbt.3300)
## Input and Output

The "deepBind" model requires DNA sequences in a fasta format, together with a BED file containing regions matching the sequence in the fasta file.


NOTE: in the fasta file, the sequence must start with 'chr' in order the model to work

Sample of FASTA file:

```
>chren13_WT
CAAAACATCCAGGTGGGCTTCAAAACAGGAAGAGAAAAAAAAGAGAGAGAAACGAAAGGAAAGAAAGAAAAGCCCAGAGTGAAGTTT
```

Sample of BED file:

```
chren13_WT    48    64
chren13_WT    56    72
chren13_WT    64    80
chren13_WT    72    88
chren13_WT    80    96
chren13_WT    88    104
```



The output of the model will be the predicted binding affinity between each enhancer sequence and the transcription factors. <br>
 This will be a numerical value that represents the likelihood of binding.

 sample of the output:
 ```
 chren13_WT	48	64	0.8379223
chren13_WT	56	72	1.0102623
chren13_WT	64	80	1.0469577
chren13_WT	72	88	1.0469577
chren13_WT	80	96	1.0469577
chren13_WT	88	104	1.0469577
 ```


## Experiment Overview

**Objective:** To assess the deepBind score for the binding sites of enhancer 13.

### Experimental Details

In our experiments, we employed the following components:

- **Enhancer 13 Segment:** We utilized a specific segment of enhancer 13, as indicated above.


- **Modified Enhancer 13 :** An altered version of enhancer 13 segment was used, featuring a 3-base pair deletion on the SOX9 binding site (bs).

### Analysis Method

We performed our analysis on a 16-base pair window with an 8-base pair shift over the sequence. deepBind analysis was conducted on three transcription factors (TFs) to evaluate their binding to the enhancer:

1. **SOX9**
2. **GATA4**
3. **AR (Androgen Receptor)**

The AR transcription factor was included as a control in our experiments.




## DeepBind Result Processing

### Result Analysis

In the analysis of the deepBind results, the following procedure was applied:

1. **Binding Score Assignment:** Each window in the sequence was assigned a binding score for a specific transcription factor (TF).

2. **Average Score Calculation:** For every position within the enhancer, we calculated the average score of its binding with the respective TF. 

This process ultimately resulted in the assignment of a binding score to each nucleotide within the enhancer sequence.

### Visualization

The graphical representation of the result is provided below:




    
![png](deepBind_analysis_files/deepBind_analysis_8_0.png)
    

