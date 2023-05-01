plot_entropy_incidence<-function(df,line_dot_size,wordsize,host,scales_x,proteinOrder,kmer_size){
  #detect if low support present
  if (proteinOrder !=""){
    #order the proteins based on user input
    level<-unlist(lapply(strsplit(proteinOrder,','),trimws))
    #set protein order as factor
    df$proteinName<-factor(df$proteinName, levels=level)
    df$size_f = factor(df$proteinName,levels = level)
  }
  
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
  
  #------------------plotting---------------------#
  #detect if low support present
  if (TRUE %in% df$lowSupport){
    df$lowSupportPos <- -0.3
    df$lowSupportPos[df$lowSupport ==TRUE]<- -0.5
    plot1<-ggplot(df) + 
    geom_rect(inherit.aes = FALSE,
             aes(xmin=position, xmax=end, ymin=-Inf, ymax=+Inf), 
             fill='#FFECAF', alpha=ifelse(df$end == -1, 0, 0.5))+
      geom_area(mapping = aes(x = position, y = entropy,color= "k-mer Entropy", linetype="k-mer Entropy"), show.legend=F)+
      geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth = (line_dot_size/10))+
      geom_point(mapping = aes(x = position,y=lowSupportPos),col=ifelse(df$lowSupportPos==-0.5, 'black', ifelse(df$lowSupportPos==-0.3, 'white', 'white')), alpha=ifelse(df$lowSupportPos==-0.5, 1, ifelse(df$lowSupportPos==-0.3, 0,0)),pch=17, linewidth =line_dot_size)+
      geom_line(mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), linewidth = (line_dot_size/10) )+ 
      geom_hline(mapping = aes( yintercept=98* maxy /100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), linewidth = (line_dot_size/10))+ 
      labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
      scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
      #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
      scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")), 
                         breaks = seq(0.0, maxy, length.out = 5),labels= sprintf(seq(0.0, maxy, length.out = 5), fmt = "%.1f")) +
      theme_classic(base_size = wordsize) + 
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
    if (host == 1){#1 host
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
      geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), linewidth = (line_dot_size/10))+
      geom_line(mapping = aes(x = position, y = totalVariants.incidence * maxy / 100, color = "Total Variants",linetype="Total Variants"), linewidth = (line_dot_size/10))+
      geom_hline(mapping = aes( yintercept=98* maxy /100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), linewidth = (line_dot_size/10))+
      labs(y = "k-mer entropy (bits)\n",x= "\nk-mer position (aa)",color = "#f7238a")+
      scale_x_continuous(limits = limits_fun,breaks = breaks_fun)+
      #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
      scale_y_continuous(sec.axis = sec_axis(~ . * 100 / maxy , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")),
                         breaks = seq(0.0, maxy, length.out = 5),labels= sprintf(seq(0.0, maxy, length.out = 5), fmt = "%.1f")) +
      theme_classic(base_size = wordsize) +
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
    if (host == 1){#one host
      plot1 +facet_grid(col=vars(df$size_f),scales = 'free',space = "free",switch = "x")
    }else{# multi host
      plot1 +facet_grid(rows = vars(df$host),col=vars(df$size_f),scales = 'free',space = "free",switch = "both")
    }
  }
}
