#Entropy and incidence of total variants
library(ggplot2)
library(gridExtra) #tutorial: https://ggplot2.tidyverse.org/reference/facet_grid.html 
library(dplyr)

#load data from csv file
data<-read.csv("../www/DiMA_HCV.csv")
df <- data.frame(data) 
#NOTE: user is required to set the protein order here
#leave the protein order as default; follow the order in csv file
proteinOrder=""

#NOTE: user is required to input the number of host here
host <- 1
#determine number of host
#scale the amino acid position for each protein
if (host ==1){ #single host
  if (proteinOrder ==""){
    a<-table(df$proteinName)
    proteinName<-as.vector(names(a))
    position<-as.vector(a)
    df$size_f = factor(df$proteinName,levels = proteinName)
  }else{
    #order the proteins based on user input
    level<-strsplit(proteinOrder, ',')[[1]]
    position<-c()
    for (i in level){
      position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
    }
    df$size_f = factor(df$proteinName,levels = level)
  } 
}else{ # multihost
  #categorise data based on host
  df$host<- factor(df$host)
  #count the aa length for each proteins (each host is expected to have same number of proteins with same length)
  df_sub<-df[df$host==unique(df$host[1]),]
  if (proteinOrder ==""){
    a<-table(df_sub$proteinName)
    proteinName<-as.vector(names(a))
    position<-as.vector(a)
    df$size_f = factor(df$proteinName,levels = proteinName)
  }else{
    #order the proteins based on user input
    level<-strsplit(proteinOrder, ',')[[1]]
    position<-c()
    for (i in level){
      position<-append(position,table(df_sub$proteinName)[names(table(df_sub$proteinName)) == i])
    }
    df$size_f = factor(df$proteinName,levels = level)
  }
}

#---------------set zero entropy region-----------------#
#check if the kmer positions are zero entropy
#Reference for getting all the positions of kmer based on kmer size and the starting position of kmer: 
#https://stackoverflow.com/questions/31120552/generate-a-sequence-of-numbers-with-repeated-intervals
#https://stackoverflow.com/questions/41637518/adding-a-shaded-rectangle-to-an-existing-plot-using-geom-rect
#determine the kmer size
kmer_size=9

#identify the kmer position with zero entropy
#position: startPosition with zero entropy
#end: startPosition + kmersize - 1
df_zeroEntropy<- df %>%
  dplyr::group_by(proteinName)%>%
  dplyr::summarize(
    end = position[which(entropy == min(entropy))+kmer_size-1], #end: startPosition + kmersize - 1
    position = position[which(entropy == min(entropy))] #extract those starting positions with zero entropy
  )%>%  
  as.data.frame()

#concatenate df with df_zeroEntropy
df <-merge(df,df_zeroEntropy,id="position",all=T)
#replace NAN in column "zero entropy" with FALSE
df$end[is.na(df$end)]<- -1

#---------------set maximum y limit-----------------#
#if the max y value in data < 10 => maxy = 10
#if max y value in data > 10 => maxy = ceiling(max y value)
if (max(df$entropy) <= 10){
  maxy <-10.0
}else if (max(df$entropy) > 10){
  maxy <- ceiling(max(df$entropy))
}

breaks_fun <- function(x) {
    seq(0,max(x),50)
}

limits_fun <- function(x) {
    c(0,max(x))
}

#----------------plotting--------------------#
#detect if low support present
if (TRUE %in% df$lowSupport){
  df$lowSupportPos <- -0.3
  df$lowSupportPos[df$lowSupport ==TRUE]<- -0.5
  plot1<-ggplot(df) + 
    geom_rect(inherit.aes = FALSE,
              aes(xmin=position, xmax=end, ymin=-Inf, ymax=+Inf), 
              fill='#FFECAF', alpha=ifelse(df$end == -1, 0, 0.5))+
    geom_area(mapping = aes(x = position, y = entropy,color= "k-mer Entropy", linetype="k-mer Entropy"), show.legend=F)+
    geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth = 0.2)+
    geom_point(mapping = aes(x = position,y=lowSupportPos),col=ifelse(df$lowSupportPos==-0.5, 'black', ifelse(df$lowSupportPos==-0.3, 'white', 'white')), alpha=ifelse(df$lowSupportPos==-0.5, 1, ifelse(df$lowSupportPos==-0.3, 0,0)),pch=17, linewidth = 0.2)+
    geom_line(mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), linewidth = 0.3 )+ 
    geom_hline(mapping = aes( yintercept=98* maxy / 100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), linewidth = 0.2)+ 
    labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
    scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
    #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
    scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")), 
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
                                   "k-mer Entropy" = "black",
                                   "Total Variants"="#f7238a"), 
                        guide = guide_legend(override.aes=aes(fill=NA)))+
    scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                      "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                      "k-mer Entropy" = 1,
                                      "Total Variants"=1))
  
  #number of host
  if (host == 1){ #one host
    plot1 +facet_grid(col=vars(df$size_f),scales = 'free',space = "free",switch = "x")
  }else{ # multi host
    plot1 +facet_grid(rows = vars(df$host),col=vars(df$size_f),scales = 'free',space = "free",switch = "both")
  }
  
}else{
  plot1<-ggplot(df) +
    geom_rect(inherit.aes = FALSE,
              aes(xmin=position, xmax=end, ymin=-Inf, ymax=+Inf), 
              fill='#FFECAF', alpha=ifelse(df$end == -1, 0, 0.5))+
  geom_area(mapping = aes(x = position, y = entropy,color= "k-mer Entropy", linetype="k-mer Entropy"),show.legend = F)+
    geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth = 0.2)+
    geom_line(mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), linewidth = 0.3 )+ 
    geom_hline(mapping = aes( yintercept=98 * maxy /100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), linewidth = 0.2)+ 
    labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
    scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
    #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
    scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")),
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
                                   "k-mer Entropy" = "black",
                                   "Total Variants"="#f7238a"), 
                        guide = guide_legend(override.aes=aes(fill=NA)))+
    scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                      "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                      "k-mer Entropy" = 1,
                                      "Total Variants"=1))
  #number of host
  if (host == 1){ #one host
    plot1 +facet_grid(col=vars(df$size_f),scales = 'free',space = "free",switch = "x")
  }else{ # multi host
    plot1 +facet_grid(rows = vars(df$host),col=vars(df$size_f),scales = 'free',space = "free",switch = "both")
  }
}

#save image, modify based on need
ggsave("plot-Entropy-and-incidence-of-totalvariants.jpg",dpi=600,height= 7.5,width = 15,unit="in",device = 'jpg')

