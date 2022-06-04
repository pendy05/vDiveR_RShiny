#load packages
if (!require("tools")) install.packages("tools")       # Install tools package
if (!require("gghalves")) install.packages("gghalves")
library("tools")  
library("gghalves")
library(Rcpp)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggtext)
if (!require("pacman")) install.packages("pacman")
if (!require("gghalves")) install.packages("gghalves")
pacman::p_load(ggplot2, Hmisc)




#read data
data<-read.csv("../www/DiMA_HCV.csv")

data<-data%>%mutate(ConservationLevel = case_when(
  data$index.incidence == 100 ~ "Completely conserved (CC)",
  data$index.incidence >= 90 ~ "Highly conserved (HC)",
  data$index.incidence >= 20 ~ "Mixed variable (MV)",
  data$index.incidence >= 10  ~ "Highly diverse (HD)",
  data$index.incidence < 10 ~ "Extremely diverse (ED)"
))

#determine the protein order; leave it as it is if you would like to follow default order in csv file
#NOTE: change this line of code below for protein order
proteinOrder=""
#NOTE: modify the following lines based on need
#conservationLabel: "1" for full label of protein; "0" for simplified label of protein
conservationLabel<-1


#plotting function
plot_plot7<- function(data){
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
    level<-c("All",strsplit(proteinOrder, ',')[[1]])
  }
  
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  #combine proteome bar with protein bars
  data<-rbind(data1,data)
  #determine the protein order
  data$level = factor(data$proteinName, levels=level)
  #calculation for total and percentage of conservation levels for each protein
  #sum up the total positions for each conservation level of proteins
  plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
  names(plot7_data)[3]<-"Total"
  
  C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")
  #check the presence of conservation level: insert value 0 if it is absent
  if (conservationLabel == 1){ #full label
    #check the presence of conservation level: insert value 0 if it is absent
    for ( conservation in C_level){ #conservation level
      for (name in level){ #proteinName
        if (!(conservation %in% plot7_data[plot7_data$proteinName==name,]$ConservationLevel)){
          plot7_data<-rbind(plot7_data,c(name,conservation,0))
        }}}
  }
  
  #sort the dataframe
  plot7_data[order(plot7_data$proteinName),]
  plot7_data$Total<- as.integer(plot7_data$Total)
  #get the percentage of each conservation level for each protein
  plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)
  
  #gather the protein label in multicolor
  plot7_data<-plot7_data%>%mutate(Label = case_when(
    plot7_data$ConservationLevel == "Completely conserved (CC)" ~ paste0(sprintf("<span style =
    'color:#000000;'>CC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
    plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
  ))
  
  #set conservation level in specific order (CC,HC,MV,HD,ED)
  plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
  
  #combine all conservation level labels into one for each protein
  Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
  #get number of protein for labelling
  nProtein<-nrow(Proteinlabel)
  
  #plotting
  ggplot(data, aes(x=level,y=index.incidence)) +
    # gghalfves
    geom_half_boxplot(outlier.shape = NA) +
    geom_half_point(aes(col = ConservationLevel), side = "r", 
                    position = position_jitter(width = 0, height=-0.7),alpha=0.6) +
    ylim(0,105) +
    labs(x=NULL, y="Index incidence (%)\n", fill="Conservation level")+
    theme_classic()+
    theme(
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.position = 'bottom',
      plot.margin = unit(c(5, 1, 1, 1), "lines"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)
    ) +
    scale_colour_manual('Conservation Level',
                        breaks = c("Completely conserved (CC)",
                                   "Highly conserved (HC)",
                                   "Mixed variable (MV)",
                                   "Highly diverse (HD)",
                                   "Extremely diverse (ED)"),
                        values = c("Completely conserved (CC)"="black",
                                   "Highly conserved (HC)"="#0057d1",
                                   "Mixed variable (MV)"="#02d57f",
                                   "Highly diverse (HD)"="#8722ff", 
                                   "Extremely diverse (ED)"="#ff617d")) +  
    geom_richtext(data = Proteinlabel, 
                  aes(x=proteinName,label = Label, y=c(rep(105,nProtein)),
                      label.size=0, label.color="transparent"),
                  position = position_dodge(width=0.1), 
                  size=2.6, color="black", hjust=0, angle=90) + 
    guides(color = guide_legend(override.aes = list(size = 2), nrow=2))+
    coord_cartesian(clip = "off")+ #allow ggtext outside of the plot
    ggtitle(unique(data$host))

}

#NOTE: modify this line of code below to determine number of host
host = 1
#single host
if (host == 1){
  plot_plot7(data)
}else{ #multihost
  
  #split the data into multiple subsets (if multiple hosts detected)
  plot7_list<-split(data,data$Host)
  plot7_multihost<-lapply(plot7_list,plot_plot7)
  
  #create spacing between multihost plots
  theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
  do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$Host))))  
}

# #save plot as 600dpi image
# #set the directory to save to
# setwd("C:\\Users\\Desktop")
 ggsave(filename="plot-Frequency-Distribution-Violin-Plots(Protein).jpg",width = 10, height = 7.5, unit="in",device='jpg', dpi=600)




  
