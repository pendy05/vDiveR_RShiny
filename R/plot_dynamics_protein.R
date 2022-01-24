plot4_5<-function(data,line_dot_size,wordsize,host,proteinOrder){
  group_names<-c("Index","Major","Minor","Unique","Total variants","Nonatypes")
  #modify the data structure
  plot4_data<-data.frame()
  #print(data)
  #single host
  if (host == 1){
    for (i in 7:12){
      tmp<-data.frame(proteinName=data[1],position=data[2],incidence=data[i],total_variants=data[11],Group=group_names[i-6],Multiindex=data[13])
      names(tmp)[3]<-"Incidence"
      names(tmp)[4]<-"Total_Variants"
      plot4_data<-rbind(plot4_data,tmp)
    }
  }else{ #multiple host
    for (i in 7:12){
      tmp<-data.frame(proteinName=data[1],position=data[2],incidence=data[i],total_variants=data[11],Group=group_names[i-6],Multiindex=data[13],host=data[14])
      names(tmp)[3]<-"Incidence"
      names(tmp)[4]<-"Total_Variants"
      plot4_data<-rbind(plot4_data,tmp)
    }
  }
  
  if (proteinOrder !=""){
    #order the proteins based on user input
    #level<-strsplit(proteinOrder, ',')[[1]]
    level<-unlist(lapply(strsplit(proteinOrder,','),trimws))
    #set protein order as factor
    #print(level)
    plot4_data$proteinName<-factor(plot4_data$proteinName, levels=level)
    plot4_data$size_f = factor(plot4_data$proteinName,levels = level)
  }
  plot5_data<-plot4_data
  #print(plot4_data$proteinName)
  #plotting
  plot4<-ggplot()+geom_point(plot4_data,mapping=aes(x=Total_Variants,y=Incidence,color=Group),alpha=1/3,size=line_dot_size)+
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
    labs(y = "Incidence (%)",x= NULL)+
    theme_classic(base_size = wordsize)+
    theme(
      legend.background = element_rect(fill = "transparent"),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      legend.position = "bottom"
    )+ guides(colour = guide_legend(override.aes = list(alpha = 1,size=line_dot_size/10+1.5),keywidth = 1,keyheight = 1,nrow=1,byrow=TRUE))+
    scale_colour_manual('',values = c("Index"="black","Total variants"="#f7238a", "Major"="#37AFAF","Minor"="#42aaff","Unique"="#af10f1","Nonatypes"="#c2c7cb" ))
  plot4<-plot4+facet_grid(col=vars(plot4_data$proteinName))
  
  #host label
  if("host" %in% colnames(data)){
    plot4<-plot4+ggtitle(unique(data$host))+
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (length(unique(data$proteinName)) <=10){
    #prepare the data for each subplot of plot5
    index<-plot5_data[plot5_data$Group %in% c("Index"),]
    major<-plot5_data[plot5_data$Group %in% c("Major"),]
    minor<-plot5_data[plot5_data$Group %in% c("Minor"),]
    unique<-plot5_data[plot5_data$Group %in% c("Unique"),]
    nonatypes<-plot5_data[plot5_data$Group %in% c("Nonatypes"),]
    variants_max_yaxis<-ceiling((max(major$Incidence,minor$Incidence,unique$Incidence)/10))*10
    
    #plot 5
    plot5_index<-ggplot(index, aes(x=proteinName, y=Incidence))+
      geom_violin(fill="black",trim = FALSE, color="black",alpha=0.9)+ylim(0,100)+ylab("Index nonamer (%)")+xlab("") +theme_bw() + 
      geom_boxplot(outlier.shape = NA,width=0.05, color="white",alpha=0.15,fill="white")+ 
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.ticks.x = element_blank())
    
    plot5_tv<-ggplot(index, aes(x=proteinName, y=Total_Variants))+
      geom_violin(fill="#f7238a",trim = FALSE, color="#f7238a",alpha=0.9)+ylim(0,100)+ylab("Total variant (%)")+xlab("") +theme_bw() + 
      geom_boxplot(outlier.shape = NA,width=0.05, color="black",alpha=0.15,fill="white")+ 
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0.1,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.ticks.x = element_blank())
    
    plot5_major<-ggplot(major, aes(x=proteinName, y=Incidence)) + 
      geom_violin(fill="#37AFAF",trim = FALSE, color="#37AFAF")+ylim(0,variants_max_yaxis)+ylab("Major variant (%)")+xlab("")+theme_bw() + 
      geom_boxplot(outlier.shape = NA,width=0.04, color="black", alpha=0.15,fill="white")+ 
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.ticks.x = element_blank(),
            axis.text.y  = element_text(face="bold"))
    
    plot5_minor<-ggplot(minor, aes(x=proteinName, y=Incidence))+
      geom_violin(fill="#42aaff",trim = FALSE,color="#42aaff")+ylim(0,variants_max_yaxis)+ylab("Minor variants (%)")+xlab("") +theme_bw() +  
      geom_boxplot(outlier.shape = NA,width=0.04, color="black", alpha=0.15,fill="white")+
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.ticks.x = element_blank(),
            axis.text.y  = element_text(face="bold"))
    
    plot5_unique<-ggplot(unique, aes(x=proteinName, y=Incidence)) + 
      geom_violin(fill="#af10f1",trim = FALSE, color="#af10f1")+ylim(0,variants_max_yaxis)+ylab("Unique variants (%)")+xlab("")+theme_bw() + 
      geom_boxplot(outlier.shape = NA,width=0.05, color="black", alpha=0.15,fill="white")+  
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.ticks.x = element_blank(),
            axis.text.y  = element_text(face="bold"))
    
    plot5_nonatypes<-ggplot(nonatypes, aes(x=proteinName, y=Incidence)) + 
      geom_violin(fill="#c2c7cb",trim = FALSE, color="#c2c7cb")+ylim(0,100)+ylab("Nonatype (%)")+xlab("")+theme_bw()+ 
      geom_boxplot(outlier.shape = NA,width=0.05, color="black", alpha=0.15,fill="white") + 
      theme_classic(base_size = wordsize)+
      theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
            axis.ticks.x = element_blank())
    plot5<-ggarrange(plot5_index,plot5_tv,plot5_nonatypes,plot5_major,plot5_minor,plot5_unique,ncol=3,nrow=2)
    
  }else{
    plot5_data$Group[plot5_data$Group == "Index"] <- "Index nonamer" 
    plot5_data$Group[plot5_data$Group == "Major"] <- "Major variant" 
    plot5_data$Group[plot5_data$Group == "Minor"] <- "Minor variants" 
    plot5_data$Group[plot5_data$Group == "Unique"] <- "Unique variants" 
    
    plot5_data$Group<-factor(plot5_data$Group, levels=c("Index nonamer","Total variants", "Nonatypes", "Major variant", "Minor variants", "Unique variants"))
    variants<-subset(plot5_data, Group=="Major variant" | Group=="Minor variants" | Group=="Unique variants")
    max_ylim<-ceiling((max(variants$Incidence)/10))*10
    
    scales_y <- list(
      "Index nonamer" = scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)),
      "Total variants" = scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)),
      "Nonatypes" = scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)),
      "Major variant" = scale_y_continuous(limits = c(0, max_ylim), breaks = seq(0, max_ylim, 10)),
      "Minor variants" = scale_y_continuous(limits = c(0, max_ylim), breaks = seq(0, max_ylim, 10)),
      "Unique variants" = scale_y_continuous(limits = c(0, max_ylim), breaks = seq(0, max_ylim, 10))
    )
    
    plot5<-ggplot()+
      geom_violin(data=plot5_data,aes(x=proteinName,y=Incidence, fill=Group, color=Group), trim=FALSE)+
      theme_classic(base_size = 8)+xlab("Protein")+ylab("Incidence (%)\n")+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position="none")+
      facet_grid_sc(rows = vars(Group),switch="y",scales = list(y = scales_y))+
      scale_colour_manual('',values = c("Index nonamer"="black","Total variants"="#f7238a", "Major variant"="#37AFAF","Minor variants"="#42aaff","Unique variants"="#af10f1","Nonatypes"="#c2c7cb" ))+
      scale_fill_manual('',values = c("Index nonamer"="black","Total variants"="#f7238a", "Major variant"="#37AFAF","Minor variants"="#42aaff","Unique variants"="#af10f1","Nonatypes"="#c2c7cb" ))
    
  }
  #plot4_5
  ggarrange(plot4,plot5,ncol=1,heights = c(1,0.5))
  
}

plot_dynamics_protein<-function(data,line_dot_size,wordsize,host,proteinOrder){
  if (host == 1){
    plot4_5(data,line_dot_size,wordsize,host,proteinOrder)
  }else{ #multihost
    #split the data into multiple subsets (if multiple hosts detected)
    plot4_list<-split(data,data$host)
    plot4_multihost<-lapply(plot4_list,plot4_5,line_dot_size,wordsize,host,proteinOrder)
    
    #create spacing between multihost plots
    theme = theme(plot.margin = unit(c(0.5,1.0,0.1,0.5), "cm"))
    do.call("grid.arrange", c(grobs=lapply(plot4_multihost,"+",theme), ncol = length(unique(data$host))))
  }
}