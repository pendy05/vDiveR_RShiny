#################
#  User Input   #  
#################
plotting7<-function(proteindata,line_dot_size, wordsize,host,proteinOrder){
  print('inside plotting 7')
  print(line_dot_size)
  p<-ggplot(proteindata,
            environment = environment()) +
    geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
    geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
                size=line_dot_size)+
    labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
    scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
    theme_classic(base_size = wordsize)+
    theme(
      #plot.title = element_text(hjust=-0.05),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.position="bottom"
    )+
    scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
                        values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
    guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
    theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
    coord_cartesian(clip = "off") #allow ggtext outside of the plot
  ggplotly(p)
}

#latest version for plotly
plot_plot7conservation_plotly<- function(data,line_dot_size,wordsize,host,proteinOrder){
  print("7.2")
  #add word 'protein' in front of each protein name
  data$proteinName<-paste("Protein",data$proteinName)
  #create data for proteome bar "All" from existing data
  data1<-data
  data1$proteinName <- "All"
  data1$level <- "All"
  #set up the order of proteins in plot from left to right
  if (proteinOrder ==""){ #follow the default order in csv file
    level<-c("All",unique(data$proteinName))
  }else{ #order the proteins based on user input
    level<-c("All",unlist(lapply(strsplit(s,','),trimws)))
  }
  print("data")
  print(data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  #combine proteome bar with protein bars
  data<-rbind(data1,data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  print("n7.3")
  #calculation for total and percentage of conservation levels for each protein
  #sum up the total positions for each conservation level of proteins
  #print(data)
  plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
  print('n7.35')
  names(plot7_data)[3]<-"Total"
  print('n7.38')

  C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")

  print("n7.48")
  #sort the dataframe
  plot7_data[order(plot7_data$proteinName),]
  plot7_data$Total<- as.integer(plot7_data$Total)
  #get the percentage of each conservation level for each protein
  plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)
  

  print("n7.49")
  #set conservation level in specific order (CC,HC,MV,HD,ED)
  plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
  proteinTotal<-length(unique(plot7_data[['proteinName']]))

  #rearrange protein order
  plot7_data<-plot7_data[order(factor(plot7_data$proteinName,levels=level)),]

  print("7.5")

    #plotting
    plot7<-ggplot(data) +
      geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
      geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
                  size=line_dot_size)+
      labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
      scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
      theme_classic(base_size = wordsize)+
      theme(
        #plot.title = element_text(hjust=-0.05),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position="bottom"
      )+
      scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
                          values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
      guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
      theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
      coord_cartesian(clip = "off") #allow ggtext outside of the plot
    list(plot7)
}

#latest version for ggplot2
#older version where geom_richtext() function is included
#why geom_richtext() is not includede in latest version? coz plotly does not support this geom_richtext()
# plotting<-function(proteindata, proteinlabel,line_dot_size, wordsize,host,proteinOrder){
#   p<-ggplot(proteindata,
#             environment = environment()) +
#     geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
#     geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
#                 size=line_dot_size)+
#     labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
#     scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
#     theme_classic(base_size = wordsize)+
#     theme(
#       #plot.title = element_text(hjust=-0.05),
#       legend.key = element_rect(fill = "transparent", colour = "transparent"),
#       legend.position="bottom"
#     )+
#     scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
#                         values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
#     geom_richtext(data = proteinlabel, aes(x=proteinName,label = Label, y=c(rep(105,length(proteinName))), label.size=0, label.color="transparent"),
#                   position = position_dodge(width=0.1),size=((wordsize/2)-2),color="black", fill="white",hjust=0,angle=90) + 
#     guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
#     theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
#     coord_cartesian(clip = "off") #allow ggtext outside of the plot
#   list(p)
# }
plot_plot7conservation<- function(data,line_dot_size,wordsize,host,proteinOrder,conservationLabel){
  print("7.2")
  #add word 'protein' in front of each protein name
  data$proteinName<-paste("Protein",data$proteinName)
  #create data for proteome bar "All" from existing data
  data1<-data
  data1$proteinName <- "All"
  data1$level <- "All"
  #set up the order of proteins in plot from left to right
  if (proteinOrder ==""){ #follow the default order in csv file
    level<-c("All",unique(data$proteinName))
  }else{ #order the proteins based on user input
    proteins<-paste("Protein",unlist(lapply(strsplit(proteinOrder,','),trimws)))
    level<-c("All",proteins)
  }
  print(level)
  print("data")
  print(data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  #combine proteome bar with protein bars
  data<-rbind(data1,data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  print("n7.3")
  #calculation for total and percentage of conservation levels for each protein
  #sum up the total positions for each conservation level of proteins
  print(data)
  plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
  print('n7.35')
  names(plot7_data)[3]<-"Total"
  print('n7.38')
  
  C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")
  if (conservationLabel == 1){ #full label
    #check the presence of conservation level: insert value 0 if it is absent
    for ( conservation in C_level){ #conservation level
      for (name in level){ #proteinName
        #print(conservation,name)
        if (!(conservation %in% plot7_data[plot7_data$proteinName==name,]$ConservationLevel)){
          plot7_data<-rbind(plot7_data,c(name,conservation,0))
        }}}
  }
  print("n7.48")
  #sort the dataframe
  plot7_data[order(plot7_data$proteinName),]
  plot7_data$Total<- as.integer(plot7_data$Total)
  #get the percentage of each conservation level for each protein
  plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)
  
  #gather the protein label in multicolor
  plot7_data<-plot7_data%>%mutate(Label = case_when(
    plot7_data$ConservationLevel == "Completely conserved (CC)" ~ paste0(sprintf("<span style =
    'color:#000000;'>CC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    
  ))
  print("n7.49")
  #set conservation level in specific order (CC,HC,MV,HD,ED)
  plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
  proteinTotal<-length(unique(plot7_data[['proteinName']]))
  #combine all conservation level labels into one for each protein
  Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
  #rearrange protein order
  plot7_data<-plot7_data[order(factor(plot7_data$proteinName,levels=level)),]
  #get number of protein for labelling
  nProtein<-nrow(Proteinlabel)
  print("nProtein")
  print(nProtein)
  print("7.5")
  #print(Proteinlabel)

    #plotting
    plot7<-ggplot(data) +
      geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
      geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
                  size=line_dot_size)+
      labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
      scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
      theme_classic(base_size = wordsize)+
      theme(
        #plot.title = element_text(hjust=-0.05),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position="bottom"
      )+
      scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
                          values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
      geom_richtext(data = Proteinlabel, aes(x=proteinName,label = Label, y=c(rep(105,nProtein)), label.size=0, label.color="transparent"),
                    position = position_dodge(width=0.1),size=((wordsize/2)-2),color="black", fill="white",hjust=0,angle=90) + 
      guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
      theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
      coord_cartesian(clip = "off") #allow ggtext outside of the plot
    plot7
}


#################
#    Sample     #
#################

plot_plot7sample<- function(data,line_dot_size,wordsize,host,proteinOrder,conservationLabel){
  print("7.2")
  
  #add word 'protein' in front of each protein name
  data$proteinName<-paste("Protein",data$proteinName)
  #create data for proteome bar "All" from existing data
  data1<-data
  data1$proteinName <- "All"
  data1$level <- "All"
  #set up the order of proteins in plot from left to right
  print("7.3")
  if (proteinOrder ==""){ #follow the default order in csv file
    print("7.35")
    level<-c("All",unique(data$proteinName))
  }else{ #order the proteins based on user input
    
    level<-c("All",unlist(lapply(strsplit(s,','),trimws)))
  }
  print("7.4")
  #print(level)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  #combine proteome bar with protein bars
  data<-rbind(data1,data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  print("7.45")
  #calculation for total and percentage of conservation levels for each protein
  #sum up the total positions for each conservation level of proteins
  plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
  names(plot7_data)[3]<-"Total"
  print("7.5")
  C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")
  #check the presence of conservation level: insert value 0 if it is absent
  #print(conservationLabel)
  if (conservationLabel == 1){ #full label
    #check the presence of conservation level: insert value 0 if it is absent
    for ( conservation in C_level){ #conservation level
      for (name in level){ #proteinName
        #print(conservation,name)
        if (!(conservation %in% plot7_data[plot7_data$proteinName==name,]$ConservationLevel)){
          plot7_data<-rbind(plot7_data,c(name,conservation,0))
        }}}
  }
  print("7.6")
  #sort the dataframe
  plot7_data[order(plot7_data$proteinName),]
  plot7_data$Total<- as.integer(plot7_data$Total)
  #get the percentage of each conservation level for each protein
  plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)
  
  #gather the protein label in multicolor
  plot7_data<-plot7_data%>%mutate(Label = case_when(
    plot7_data$ConservationLevel == "Completely conserved (CC)" ~ paste0(sprintf("<span style =
    'color:#000000;'>CC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
    
  ))
  
  #set conservation level in specific order (CC,HC,MV,HD,ED)
  plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
  
  #combine all conservation level labels into one for each protein
  Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
  #get number of protein for labelling
  print("label")
  print(Proteinlabel)
  nProtein<-nrow(Proteinlabel)
  print("nProtein")
  print(nProtein)
  print("7.5")
  #plotting
  plot7<-ggplot(data) +
    geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
    geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
                size=line_dot_size)+
    labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
    scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
    theme_classic(base_size = wordsize)+
    theme(
      #plot.title = element_text(hjust=-0.05),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.position="bottom"
    )+
    scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
                        values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
    geom_richtext(data = Proteinlabel, aes(x=proteinName,label = Label, y=c(rep(105,nProtein)), label.size=0, label.color="transparent"),
                  position = position_dodge(width=0.1),size=((wordsize/2)-2),color="black", fill="white",hjust=0,angle=90) + 
    guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
    theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
    coord_cartesian(clip = "off") #allow ggtext outside of the plot
  plot7
  #grid.arrange(plot7,top=unique(data$host))
}
plot_conservationLevel_sample<-function(data,line_dot_size,wordsize,host,proteinOrder,conservationLabel){
  print("plot7reactive")
  data<- data.frame(data)
  data<-data%>%mutate(ConservationLevel = case_when(
    data$index.incidence == 100 ~ "Completely conserved (CC)",
    data$index.incidence >= 90 ~ "Highly conserved (HC)",
    data$index.incidence >= 20 ~ "Mixed variable (MV)",
    data$index.incidence >= 10  ~ "Highly diverse (HD)",
    data$index.incidence < 10 ~ "Extremely diverse (ED)"
    
  ))
  print("7.0")
  #single host
  if (host == 1){
    plot_plot7sample(data,line_dot_size,wordsize,host,proteinOrder,conservationLabel)
  }else{#multihost
    data$host = factor(data$host)
    #split the data into multiple subsets (if multiple hosts detected)
    plot7_list<-split(data,data$host)
    plot7_multihost<-lapply(plot7_list,plot_plot7sample,line_dot_size,wordsize,host,proteinOrder,conservationLabel)
    
    #create spacing between multihost plots
    theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
    do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$host))))
    
  }
}