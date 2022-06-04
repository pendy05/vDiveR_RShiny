#Correlation plot
#load package
library(ggplot2)

#read in data
data<-read.csv("../www/DiMA_HCV.csv")

#load thed data into a dataframe
df <- data.frame(data) 

#plot the scatter plot with density
plot2<-ggplot(df)+geom_point(mapping = aes(x=totalVariants.incidence,y=entropy),alpha=1/3,size=3)+
  labs(y = "k-mer entropy (bits)\n",x= "\nTotal variants (%)")+
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
  scale_y_continuous(limits = c(0, ceiling(max(df$entropy))), breaks = seq(0, ceiling(max(df$entropy)), 0.5))+
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "#000000", fill=NA, size=1)
  )
#NOTE: modify this line of code below to determine number of host
host<-1

if (host == 1){ #single host
  #plot the scatter plot with density
  plot2
}else{ #multiple host
  plot2+facet_grid_sc(rows = vars(df$host),space = "free",switch = "x")
}
#save the image with 600dpi, modify based on need
ggsave(filename="plot-Correlation-Plot.jpg",width = 5.5, height = 5.5, unit="in",device='jpg', dpi=600)
