library(ggpubr)
plot_plot3<- function(data,line_dot_size,wordsize,host){
  
  #modify data structure
  plot3_data<-data.frame()
  group_names<-c("Index","Major","Minor","Unique","Total variants","Distinct variants")
  
  for (i in 7:12){
    tmp<-data.frame(proteinName=data[1],position=data[2],incidence=data[i],total_variants=data[11],Group=group_names[i-6],Multiindex=data[13])
    names(tmp)[3]<-"Incidence"
    names(tmp)[4]<-"Total_Variants"
    plot3_data<-rbind(plot3_data,tmp)
  }
  plot3b_data<-plot3_data
  minor<-rbind(plot3_data[plot3_data$Group == "Index",],plot3_data[plot3_data$Group == "Total variants",])
  uniq<-rbind(plot3_data[plot3_data$Group == "Index",],plot3_data[plot3_data$Group == "Total variants",])
  minor$motif<- "Minor"
  uniq$motif<-"Unique"
  
  #plot3a
  plot3_data<-plot3_data%>%mutate(motif = case_when(
    plot3_data$Group == "Index" ~ "Major",
    plot3_data$Group == "Total variants" ~ "Major",
    plot3_data$Group == "Major" ~ "Major",
    plot3_data$Group == "Minor"  ~ "Minor",
    plot3_data$Group == "Unique" ~ "Unique",
    plot3_data$Group == "Distinct variants" ~ "Distinct variants"
    
  ))
  
  plot3_data<- rbind(plot3_data,minor,uniq)
  plot3_data$motif<-factor(plot3_data$motif,levels = c("Major","Minor","Unique","Distinct variants"))
  
  if (host == 1){ #one host
    ROW=1
  }else{ 
    ROW=2
  }
  
  plot3a<-ggplot()+geom_point(plot3_data,mapping=aes(x=Total_Variants,y=Incidence,color=Group),alpha=1/3,size=line_dot_size)+
    geom_point(plot3_data,mapping = aes(x =Total_Variants,y=Incidence),col=ifelse(plot3_data$multiIndex== TRUE & plot3_data$Group== "Index", 'red', ifelse(plot3_data$multiIndex== FALSE, 'white', 'white')), alpha=ifelse(plot3_data$multiIndex ==TRUE & plot3_data$Group== "Index", 1, ifelse(plot3_data$multiIndex== TRUE, 0,0)),pch=1,size=line_dot_size,stroke=1.05)+
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
    theme_classic(base_size = wordsize)+
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      strip.text.x = element_blank(),
      legend.position="bottom")+
    labs(y= "Incidence (%)", x="\nTotal variants (%)")+
    facet_wrap(~ motif,ncol = 1)+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=line_dot_size/10+1.5),nrow = ROW))+
    scale_colour_manual('',breaks=c("Index","Total variants","Major","Minor","Unique","Distinct variants"),
                        values = c("Index"="black", "Total variants"="#f7238a","MultiIndex"="red",
                                   "Major"="#37AFAF" , "Minor"="#42aaff","Unique"="#af10f1", "Distinct variants"="#c2c7cb"))
  #host label
  if("host" %in% colnames(data)){
    plot3a<-plot3a+ggtitle(unique(data$host))+theme(plot.title = element_text(hjust = 0.5))
  }
  
  #plot3b
  index<-plot3b_data[plot3b_data$Group %in% "Index",]
  nonatypes<-plot3b_data[plot3b_data$Group %in% "Distinct variants",]
  variants<-plot3b_data[plot3b_data$Group %in% c("Major","Minor","Unique"),]
  variants$x<-"x"
  variants_max_yaxis<-ceiling((max(as.numeric(variants$Incidence))/10))*10
  
  plot3b_index<-ggplot(index, aes(x=Group,y=Incidence))+geom_violin(color="black",fill="black",scale="width")+geom_boxplot(width=0.08,alpha=0.20,fill="white",outlier.shape=NA,color="white")+
    ylim(c(0,100))+
    labs(y=NULL,x="Index")+
    theme_classic(base_size = wordsize)+
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank())+ 
    scale_color_grey()
  
  plot3b_tv<-ggplot(index, aes(x=Group,y=Total_Variants))+geom_violin(color="#f7238a",fill="#f7238a")+geom_boxplot(width=0.08,alpha=0.20,color="black",fill="white",outlier.shape=NA)+
    ylim(c(0,100))+
    labs(y=NULL,x="Total Variants")+
    theme_classic(base_size = wordsize)+
    theme(
      plot.margin = margin( t=5,
                            b=5,
                            r = -0.25),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text=element_text(colour="white"),
      axis.text.x  = element_blank(),
      axis.ticks = element_blank())
  
  plot3b_nonatype<-ggplot(nonatypes, aes(x=Group,y=Incidence))+geom_violin(color="#c2c7cb",fill="#c2c7cb")+geom_boxplot(width=0.08,alpha=0.20,fill="white",outlier.shape=NA)+
    ylim(c(0,100))+
    labs(y=NULL,x="Distinct variants")+
    theme_classic(base_size = wordsize)+
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text=element_text(colour="white"),
      axis.text.x  = element_blank(),
      axis.ticks = element_blank())
  
  plot3b_variants<-ggplot(variants)+geom_violin(mapping=aes(x=x,y=Incidence,color=Group,fill=Group))+
    geom_boxplot(mapping = aes(x=x,y=Incidence),width=0.08,alpha=0.20,fill="white",outlier.shape=NA)+
    scale_y_continuous(position = "right",limits = c(0,variants_max_yaxis))+
    theme_classic(base_size = wordsize-2)+
    labs(y=NULL,x=NULL)+
    facet_wrap(Group ~ .,ncol=1,strip.position ="right")+ 
    theme(strip.placement = "outside",
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x=element_blank (), 
          axis.ticks.x=element_blank (),
          panel.spacing = unit(0, "lines"),
          strip.background=element_blank (),
          axis.text.y = element_text(size=wordsize/2+1),
          plot.margin = margin(t=5,l = -0.1)
    )+
    scale_colour_manual('',values = c( "Major"="#37AFAF","Minor"="#42aaff","Unique"="#af10f1" ))+
    scale_fill_manual('',values = c( "Major"="#37AFAF","Minor"="#42aaff","Unique"="#af10f1" ))+
    guides(fill="none",color='none') 
  
  plot3b_variants<-annotate_figure(plot3b_variants,bottom = text_grob("a\n",size=ceiling(wordsize/2-1) ,color = "white"))
  
  plot3b<-ggarrange(plot3b_index, plot3b_tv,plot3b_variants, plot3b_nonatype, ncol=4,widths = c(1,0.9,0.95,1))
  plot3b<-annotate_figure(plot3b,left = text_grob("Incidence (%)",size=wordsize,rot=90,hjust=0.3))
  
  #plot3
  ggarrange(plot3a,plot3b,ncol=1,heights = c(1,0.3))
}

plot_dynamics_proteome<-function(data,line_dot_size,wordsize,host){
  if (length(unique(data$proteinName)) >1){
    #single host
    if (host == 1){
      plot_plot3(data,line_dot_size,wordsize,host)
    }else{#multihost
      
      #split the data into multiple subsets (if multiple hosts detected)
      plot3_list<-split(data,data$host)
      plot3_multihost<-lapply(plot3_list,plot_plot3,line_dot_size,wordsize,host)
      
      #create spacing between multihost plots
      theme = theme(plot.margin = unit(c(0.5,1.0,0.1,0.5), "cm"))
      do.call("grid.arrange", c(grobs=lapply(plot3_multihost,"+",theme), ncol = length(unique(data$host))))
      
    }
  }
}
