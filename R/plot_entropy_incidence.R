plot_entropy_incidence<-function(df,line_dot_size,wordsize,host,scales_x,proteinOrder){
  #detect if low support present
  if (proteinOrder !=""){
    #order the proteins based on user input
    level<-unlist(lapply(strsplit(proteinOrder,','),trimws))
    #set protein order as factor
    df$proteinName<-factor(df$proteinName, levels=level)
    df$size_f = factor(df$proteinName,levels = level)
  }
  if (TRUE %in% df$lowSupport){
    print("lowSupport TRUE")
    df$lowSupportPos <- -0.3
    df$lowSupportPos[df$lowSupport ==TRUE]<- -0.5
    plot1<-ggplot(df) + 
      geom_area(df,mapping = aes(x = position, y = entropy,color= "Nonamer Entropy", linetype="Nonamer Entropy"), show.legend=F)+
      geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"), size= (line_dot_size/10))+
      geom_point(df,mapping = aes(x = position,y=lowSupportPos),col=ifelse(df$lowSupportPos==-0.5, 'black', ifelse(df$lowSupportPos==-0.3, 'white', 'white')), alpha=ifelse(df$lowSupportPos==-0.5, 1, ifelse(df$lowSupportPos==-0.3, 0,0)),pch=17, size=line_dot_size)+
      geom_line(df,mapping = aes(x = position, y = totalVariants.incidence *10 / 100, color = "Total Variants",linetype="Total Variants"), size= (line_dot_size/10) )+ 
      geom_hline(df,mapping = aes( yintercept=98*10/100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), size= (line_dot_size/10))+ 
      labs(y = "Nonamer entropy (bits)\n",x= "\nNonamer position (aa)",color = "#f7238a")+
      #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
      scale_y_continuous(sec.axis = sec_axis(~ . * 100 / 10 , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")), 
                         breaks = c(0.0,2.5,5.0,7.5,10.0),labels=c("0.0","2.5","5.0","7.5","10.0")) + 
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
                                     "Nonamer Entropy" = "black",
                                     "Total Variants"="#f7238a"), 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                        "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                        "Nonamer Entropy" = 1,
                                        "Total Variants"=1))
    
    
    #number of host
    if (host == 1){#1 host
      plot1 +facet_grid_sc(col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "x")
    }else{# multi host
      plot1 +facet_grid_sc(rows = vars(df$host),col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "both")
    }
    
  }else{
    
    plot1<-ggplot(df) +
      geom_area(df,mapping = aes(x = position, y = entropy,color= "Nonamer Entropy", linetype="Nonamer Entropy"),show.legend = F)+
      geom_hline(mapping = aes(yintercept=9.2, color = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)", linetype = "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"),size= (line_dot_size/10))+
      geom_line(df,mapping = aes(x = position, y = totalVariants.incidence *10 / 100, color = "Total Variants",linetype="Total Variants"),  size= (line_dot_size/10))+
      geom_hline(df,mapping = aes( yintercept=98*10/100, color = "Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)",linetype ="Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"), size= (line_dot_size/10))+
      labs(y = "Nonamer entropy (bits)\n",x= "\nNonamer position (aa)",color = "#f7238a")+
      #how to second y-axis: https://whatalnk.github.io/r-tips/ggplot2-rbind.nb.html
      scale_y_continuous(sec.axis = sec_axis(~ . * 100 / 10 , name = "Total variants (%)",breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100")),
                         breaks = c(0.0,2.5,5.0,7.5,10.0),labels=c("0.0","2.5","5.0","7.5","10.0")) +
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
                                     "Nonamer Entropy" = "black",
                                     "Total Variants"="#f7238a"), 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_linetype_manual("",values=c("Reference: Maximum Total Variants (98%) for HIV-1 Clade B (Env Protein)"=5,
                                        "Reference: Maximum Entropy (9.2) for HIV-1 Clade B (Env Protein)"=5,
                                        "Nonamer Entropy" = 1,
                                        "Total Variants"=1))
    #number of host
    if (host == 1){#one host
      plot1 +facet_grid_sc(col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "x")
    }else{# multi host
      print("multihost")
      plot1 +facet_grid_sc(rows = vars(df$host),col=vars(df$size_f),scales = list(x = scales_x),space = "free",switch = "both")
    }
  }
}
