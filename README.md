# DiveR (Diversity dynamics visualization in R)
DiveR is a graphical user interface (GUI)-based web application hosted on R Shiny for the visualization of various diversity dynamics. The application allow users to visualize the diversity motifs (index and its total variants – major, minor and unique) for elucidation of the underlying inherent dynamics by generating multiple plots: <br>
1. entropy and incidence of total variants, 
2. relationship between entropy and total variants for k-mer positions of a viral protein, 
3. dynamics of diversity motifs for viral proteome, 
4. dynamics of diversity motifs (protein) and 
5. distribution of conservation levels of protein k-mer positions for viral proteome and protein(s).

Please visit our [DiveR web server](https://protocol-viral-diversity.shinyapps.io/DiveR/)!

## Background
Sequence diversity, as a result of various evolutionary forces, challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. A publicly available tool, [Diversity Motif Analyser (DiMA)](https://github.com/PU-SDS/DiMA), has been developed to facilitate the dissection of protein sequence diversity dynamics for viruses, the understanding of which is critical to develop effective intervention strategies. In order to ease those without programming knowledge to perform analytics on the DiMA JSON output, we, herein, provide DiveR, an extension to DiMA in the form of a graphical user interface (GUI)-based  web application hosted on R Shiny. 

<!--
To-do list
1. Installation (Please refer to this: https://github.com/rhondabacher/methylscaper/)
2. Dependencies 
3. Q: is the server support big data -> May test with SARS-CoV-2 spike protein that you will work on for Melike's project. 
Just a note: If not support, then you may suggest a way to user - install and run locally (Please refer to this: https://github.com/alpreyes/GENAVi)
-->

## Contents
### Scripts
- app.R : call ui.R and server.R
- ui.R: user interface (frontend)
- server.R: server side; deal with reactivity, calculations, plotting etc (backend)

### Folders
- www/: consists of images and sample dataset
- R/: consists of R functions embedded as part of DiveR 
- R-individual-scripts/: consists of individual R scripts used for producing plots and table in DiveR (Note: similar content as R/ but if you would like to access the R functions independently without using the DiveR web server, please refer to scripts in this folder!)

## Bug or Questions
Please leave a message on our [GitHub Issues](https://github.com/pendy05/DiveR/issues) or contact Asif M. Khan (asif@perdanauniversity.edu.my).

## Authors
The following individuals have contributed code to DiveR:
- Pendy Tok
- Li Chuin Chong
- Evgenia Chikina
