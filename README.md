# DiveR (Diversity dynamics Visualization in R)
Sequence diversity, as a result of various evolutionary forces, challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. A publicly available tool, Diversity Motif Analyser (DiMA; https://github.com/PU-SDS/DiMA), has been developed to facilitate the dissection of protein sequence diversity dynamics for viruses, the understanding of which is critical to develop effective intervention strategies. Herein, we provide DiveR, an extension to DiMA in the form of a graphical user interface (GUI)-based  web application hosted on R Shiny to ease the visualisation of the various diversity dynamics. This is particularly useful for those without programming knowledge to perform analytics on the DiMA JSON output. DiveR allows users to visualize the diversity motifs (index and its total variants – major, minor and unique) for elucidation of the underlying inherent dynamics. The application generates multiple plots: <br>
1. entropy and incidence of total variants, 
2. relationship between entropy and total variants for k-mer positions of a viral protein, 
3. dynamics of diversity motifs for viral proteome, 
4. dynamics of diversity motifs (protein) and 
5. distribution of conservation levels of protein k-mer positions for viral proteome and protein(s).

## Contents
### Scripts
- ui.R: user interface (frontend)
- server.R: server side; deal with reactivity, calculations, plotting etc (backend)

### Folders
- www/: consists of images and sample dataset
- R/: consists of individual plotting scripts
