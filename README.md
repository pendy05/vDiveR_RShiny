# DiveR (Diversity dynamics Visualization in R)
Sequence diversity, as a result of various evolutionary forces, challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. A publicly available tool, Diversity Motif Analyser (DiMA; https://github.com/PU-SDS/DiMA), has been developed to facilitate the dissection of protein sequence diversity dynamics for viruses, which is critical to develop effective intervention strategies. The development of user-friendly tools may aid biologists to analyse the sequence diversity; however, the need for bioinformatics or programming skills to interpret the analysis challenges them. Hence, we provide an extension to DiMA in the form of a graphical user interface (GUI) application hosted on R Shiny to ease the complex workflows, particularly to those without programming knowledge.

## Contents
### Scripts
- ui.R: user interface (frontend)
- server.R: server side; deal with reactivity, calculations, plotting etc (backend)

### Folders
- www/: consists of images and sample dataset
- R/: consists of individual plotting scripts
