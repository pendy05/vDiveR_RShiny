plot_correlation<-function(df,line_dot_size,wordsize,host){
  plot2<-ggplot(df)+geom_point(mapping = aes(x=totalVariants.incidence,y=entropy),alpha=1/3,size=line_dot_size)+
    labs(y = "k-mer entropy (bits)\n",x= "\nTotal variants (%)")+
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
    scale_y_continuous(limits = c(0, ceiling(max(df$entropy))), breaks = seq(0, ceiling(max(df$entropy)), 0.5))+
    theme_classic(base_size = wordsize)+
    theme(
      panel.border = element_rect(colour = "#000000", fill=NA, size=1)
    )
  
  #plot the scatter plot with density
  if (host == 1){ #one host
    plot2
  }else{ #multiple host
    plot2+facet_grid(rows = vars(df$host),space = "free",switch = "x")
  }
}