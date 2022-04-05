#Entropy and incidence of total variants
library(ggplot2)
library(gridExtra) #tutorial: https://ggplot2.tidyverse.org/reference/facet_grid.html 
library(facetscales) #https://stackoverflow.com/a/54074323/13970350
library(dplyr)

#load data from csv file
data<-read.csv("HCV_proteins.csv")
df <- data.frame(data) 
#NOTE: user is required to set the protein order here
#leave the protein order as default; follow the order in csv file
proteinOrder=""

#NOTE: user is required to input the number of host here
host <-1
#determine number of host
#scale the amino acid position for each protein
if (host ==1){ #single host
  if (proteinOrder ==""){
    a<-table(df$proteinName)
    proteinName<-as.vector(names(a))
    position<-as.vector(a)
    df$size_f = factor(df$proteinName,levels = proteinName)
    scales_x<-mapply(function(x,y){
      x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
    }, proteinName,position)
  }else{
    #order the proteins based on user input
    level<-strsplit(proteinOrder, ',')[[1]]
    position<-c()
    for (i in level){
      position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
    }
    df$size_f = factor(df$proteinName,levels = level)
    scales_x<-mapply(function(x,y){
      x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
    }, level,position)
  } 
}else{ # multihost
  #categorise data based on host
  df$Host<- factor(df$Host)
  #count the aa length for each proteins (each host is expected to have same number of proteins with same length)
  df_sub<-df[df$Host==unique(df$Host[1]),]
  if (proteinOrder ==""){
    a<-table(df_sub$proteinName)
    proteinName<-as.vector(names(a))
    position<-as.vector(a)
    df$size_f = factor(df$proteinName,levels = proteinName)
    scales_x<-mapply(function(x,y){
      x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
    }, proteinName,position)
  }else{
    #order the proteins based on user input
    level<-strsplit(proteinOrder, ',')[[1]]
    position<-c()
    for (i in level){
      position<-append(position,table(df_sub$proteinName)[names(table(df_sub$proteinName)) == i])
    }
    df$size_f = factor(df$proteinName,levels = level)
    scales_x<-mapply(function(x,y){
      x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
    }, level,position)
  }
}

#---------------zero entropy-----------------#
#check if the kmer positions are zero entropy
#Reference for getting all the positions of kmer based on kmer size and the starting position of kmer: 
#https://stackoverflow.com/questions/31120552/generate-a-sequence-of-numbers-with-repeated-intervals
#determine the kmer size
kmer_size=9

#identify the kmer position with zero entropy
pos_withZeroEntropy<- unique(sequence(kmer_size) + rep(df[df$entropy==0,]$position-1, rep(kmer_size,length(df[df$entropy==0,]$position))))
#fill in "TRUE" for position with zero entropy
df$zeroEntropy <- ifelse(df$position %in% pos_withZeroEntropy, TRUE, FALSE)

#---------------set maximum y limit-----------------#
#if the max y value in data < 10 => maxy = 10
#if max y value in data > 10 => maxy = ceiling(max y value)
if (max(df$entropy) <= 10){
  maxy <-10.0
}else if (max(df$entropy) > 10){
  maxy <- ceiling(max(df$entropy))
}

#----------------plotting--------------------#
#detect if low support present
if (TRUE %in% df$lowSupport){
  df$lowSupportPos <- -0.3
  df$lowSupportPos[df$lowSupport ==TRUE]<- -0.5
  plot1<-ggplot(df) + 
    geom_vline(mapping = aes(xintercept = position), color='#FFECAF',alpha = ifelse(df$zeroEntropy == TRUE, 1, 0.0))+
    geom_area(df,mapping = aes(x = position, y = entropy,color= "Nonamer Entropy", linetype="Nonamer Entropy"), show.legend=F)+
    geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), size= 0.2)+
    geom_point(df,mapping = aes(x = position,y=lowSupportPos),col=ifelse(df$lowSupportPos==-0.5, 'black', ifelse(df$lowSupportPos==-0.3, 'white', 'white')), alpha=ifelse(df$lowSupportPos==-0.5, 1, ifelse(df$lowSupportPos==-0.3, 0,0)),pch=17)+
    geom_line(df,mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), size= 0.3 )+ 
    geom_hline(df,mapping = aes( yintercept=98* maxy / 100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), size= 0.2)+ 
    labs(y = "Nonamer entropy (bits)\n",x= "\nNonamer position (aa)",color = "#f7238a")+
    #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
    scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)\n",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")), 
                       breaks = seq(0.0, maxy, length.out = 5),labels= sprintf(seq(0.0, maxy, length.out = 5), fmt = "%.1f")) + 
    theme_classic() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line.y.right = element_line(color = "#f7238a"), 
      axis.ticks.y.right = element_line(color = "#f7238a"),
      axis.title.y.right =  element_text(color = "#f7238a"),
      legend.position="bottom"
    )+
    scale_colour_manual("",
                        values = c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"="#f7238a",
                                   "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"="black",
                                   "Nonamer Entropy" = "black",
                                   "Total Variants"="#f7238a"), 
                        guide = guide_legend(override.aes=aes(fill=NA)))+
    scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                      "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                      "Nonamer Entropy" = 1,
                                      "Total Variants"=1))
  
  #number of host
  if (host == 1){ #one host
    plot1 +facet_grid_sc(col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "x")
  }else{ # multi host
    plot1 +facet_grid_sc(rows = vars(df$Host),col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "both")
  }
  
}else{
  plot1<-ggplot(df) +
    geom_vline(mapping = aes(xintercept = position), color='#FFECAF',alpha = ifelse(df$zeroEntropy == TRUE, 1, 0.0))+
    geom_area(df,mapping = aes(x = position, y = entropy,color= "Nonamer Entropy", linetype="Nonamer Entropy"),show.legend = F)+
    geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), size= 0.2)+
    geom_line(df,mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), size= 0.3 )+ 
    geom_hline(df,mapping = aes( yintercept=98 * maxy /100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), size= 0.2)+ 
    labs(y = "Nonamer entropy (bits)\n",x= "\nNonamer position (aa)",color = "#f7238a")+
    #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
    scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)\n",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")),
                       breaks = seq(0.0, maxy, length.out = 5),labels= sprintf(seq(0.0, maxy, length.out = 5), fmt = "%.1f")) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y.right = element_line(color = "#f7238a"),
      axis.ticks.y.right = element_line(color = "#f7238a"),
      axis.title.y.right =  element_text(color = "#f7238a"),
      legend.position="bottom",
    )+
    scale_colour_manual("",
                        values = c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"="#f7238a",
                                   "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"="black",
                                   "Nonamer Entropy" = "black",
                                   "Total Variants"="#f7238a"), 
                        guide = guide_legend(override.aes=aes(fill=NA)))+
    scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                      "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                      "Nonamer Entropy" = 1,
                                      "Total Variants"=1))
  #number of host
  if (host == 1){ #one host
    plot1 +facet_grid_sc(col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "x")
  }else{ # multi host
    plot1 +facet_grid_sc(rows = vars(df$Host),col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "both")
  }
}

#save image, modify based on need
ggsave("plot-Entropy-and-incidence-of-totalvariants.jpg",dpi=600,height= 7.5,width = 15,unit="in",device = 'jpg')

