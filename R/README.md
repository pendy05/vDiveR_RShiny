[To the main README >](https://github.com/pendy05/DiveR/blob/main/README.md)

# R code for viral diversity dynamics plots
These five R scripts are to plot the viral diversity and showcase the patterns of the dynamics of the protein sequence change and the frequency distribution of each of the diversity motifs of viruses.<br><br>

## Table of contents
- [Plot 1: Entropy and Incidence of Total Variants](#plot-1-entropy-and-incidence-of-total-variants)
- [Plot 2: Correlation of Entropy](#plot-2-correlation-of-entropy)
- [Plot 3: Dynamics of Diversity Motifs (Proteome)](#plot-3-dynamics-of-diversity-motifs-proteome)
- [Plot 4: Dynamics of Diversity Motifs (Protein)](#plot-4-dynamics-of-diversity-motifs-protein)
- [Plot 5: Frequency Distribution Violin Plots (Protein)](#plot-5-frequency-distribution-violin-plots-protein)

<br>

## Plot 1: Entropy and Incidence of Total Variants
**R script:** [plot_entropy_incidence.R](https://github.com/pendy05/DiveR/blob/main/R/plot_entropy_incidence.R) 
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132232581-ed1a8d49-aefb-4a7c-89bd-09053065e198.jpg" width="1000" height="500">
</p>

**Description:** <br>
Entropy (black) and incidence of total variants (pink) were measured for each aligned nonamer (nine amino acids) position (1-9, 2-10, etc.) of the proteins. The entropy values indicate the level of variability at the corresponding nonamer positions, with zero representing completely conserved positions (total variants incidence of 0%). A reference for each entropy (black dotted line) and incidence of total variants (pink dotted line) were indicated herein for comparison. 
<br><br>

## Plot 2: Correlation of Entropy
**R script:** [plot_correlation.R](https://github.com/pendy05/DiveR/blob/main/R/plot_correlation.R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132232809-2bc5a305-5570-4528-ba41-9819a6f67a09.jpg" width="500" height="500">
</p>

**Description:** <br>
Relationship between incidence of total variants and entropy for viral proteome nonamer positions.
<br><br>

## Plot 3: Dynamics of Diversity Motifs (Proteome)
**R script:** [plot_dynamics_proteome.R](https://github.com/pendy05/DiveR/blob/main/R/plot_dynamics_proteome.R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132233240-e63f6f3f-c0cf-4e73-9879-adcf4dac4d44.jpg" width="500" height="600">
</p>

**Description:** <br>
Nonamers (peptide sequences of nine) are classified into four different motifs, namely index, major, minor and unique, based on their incidences. Nonatypes defines as the distinct nonamers for a given position. The diversity of the position was depicted by the decline of the index incidences (black) and the increase of total variant incidences (pink).
<br><br>

## Plot 4: Dynamics of Diversity Motifs (Protein)
**R script:** [plot_dynamics_protein.R](https://github.com/pendy05/DiveR/blob/main/R/plot_dynamics_protein.R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132233808-fa6dd8b4-12c8-4823-ad3f-01271c30550f.jpg" width="500" height="600">
</p>

**Description:** <br>
Nonamers (peptide sequences of nine) are classified into four different motifs, namely index, major, minor and unique, based on their incidences. Nonatypes defines as the distinct nonamers for a given position. The diversity of the position was depicted by the decline of the index incidences (black) and the increase of total variant incidences (pink).
<br><br>

## Plot 5: Frequency Distribution Violin Plots (Protein)
**R script:** [plot_conservationLevel.R](https://github.com/pendy05/DiveR/blob/main/R/plot_conservationLevel.R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132236375-0639ab97-1452-46ed-b18e-26540a96ead8.jpg" width="1000" height="600">
</p>

**Description:** <br>
The nonamer positions of the proteome and the individual proteins were defined as completely conserved (black) ( 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 = 100% ), highly conserved (blue) (90% ≤ 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 < 100%), mixed variable (green) ( 20% < 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 90%), highly diverse (purple) (10% < 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 20%) and extremely diverse (pink) ( 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 10% ). 

