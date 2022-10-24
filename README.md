# DiveR (Diversity dynamics visualization in R)
DiveR is a graphical user interface (GUI)-based web application hosted on R Shiny for the visualization of various diversity dynamics. The application allow users to visualize the diversity motifs (index and its variants – major, minor and unique) for elucidation of the underlying inherent dynamics. DiveR presents a total of eight tabs: 1) homepage, 2) data description, with tabs 3) to 7) presenting five plots depicting sequence variability dynamics and lastly 8) help page tab. DiveR generates five plots for k-mer positions of a viral protein/proteome:
1. entropy and incidence of total variants, 
2. relationship between entropy and total variants for <i>k</i>-mer positions of a viral protein, 
3. dynamics of diversity motifs for viral proteome, 
4. dynamics of diversity motifs for individual protein and 
5. distribution of conservation levels (completely conserved, highly con-served, mixed variable, highly diverse, and extremely diverse) of <i>k</i>-mer positions for viral proteome and protein(s).

Please visit our [DiveR web server](https://protocol-viral-diversity.shinyapps.io/DiveR/)!
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/168590154-9870f6d9-12d2-47d2-b7ed-8150104dce91.jpg" width=80% height=80%>
</p>
Figure: Overview of the DiveR R shiny app. DiveR consists of eight tabs where each tab includes a brief description of the implemented functionalities. Users must provide either an aligned sequence file or a DiMA JSON output file in tab 2 to initiate the visualization of dynamics in sequence change in the form of plots presented in tabs 3 to 7.

## Table of Contents
- [Background](#background)
- [Local Setup](#local-setup)
- [DiveR Contents](#diver-contents)
  - [Scripts](#scripts)
  - [Folders](#folders)

## Background
Viruses are one of the main contributors to the global burden of infectious-related mortality and disability. Sequence diversity, as a result of various evolutionary forces, challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. A publicly available tool, [Diversity Motif Analyser (DiMA)](https://github.com/PU-SDS/DiMA), has been developed to facilitate the dissection of protein sequence diversity dynamics for viruses, the understanding of which is critical to develop effective intervention strategies. DiMA quantifies the sequence diversity using Shannon’s entropy for each aligned overlapping <i>k</i>-mer positions, distributes the <i>k</i>-mers into four diversity motifs (index, major, minor and unique) and stores this information in JSON format. However, interpretation and analysis of data stored in JSON data might be a challenging task to biologists who have limited or no knowledge of bioinformatics or programming background. Herein, we present [DiveR](https://protocol-viral-diversity.shinyapps.io/DiveR), a DiMA wrapper implemented as a web-based application, hosted on R Shiny server, to ease the visualization of outputs from DiMA. 

<!--
To-do list
1. Installation (Please refer to this: https://github.com/rhondabacher/methylscaper/)
2. Dependencies 
3. Q: is the server support big data -> May test with SARS-CoV-2 spike protein that you will work on for Melike's project. 
Just a note: If not support, then you may suggest a way to user - install and run locally (Please refer to this: https://github.com/alpreyes/GENAVi)
-->

## Local Setup
### Requirement: Python Virtual Environment
1. Install “virtualenv” python application cmd package
```
$ sudo pip install virtualenv
```
2. Create a new virtual environment in a folder called “python_env” within your project directory
```
$ virtualenv python_env
```
3. Activate the virtualenv
```
# Windows
$ source python_env/Scripts/activate
```
4. Verify that you have activated the correct version of Python (show show absolute path to your current working directory)
 ```
 $ which python
 ```
5. Install the numpy & dima-cli python package (NOTE: numpy is a must download package as it is required for reticulate to move data between R and Python) 
```
$ pip install numpy dima-cli==4.1.1 pandas typing-extensions
```


## DiveR Contents
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
