plot_wp<-function(data,wordsize){
  colnames(data) <- c('region','count')
  world_map <- map_data("world")
  gg <- ggplot(world_map, aes(x = long, y = lat, group = group)) + geom_polygon(fill="lightgray", colour = "#888888")
  pathogens.map <- left_join(data, world_map, by = "region")
  gg <- gg + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
             scale_fill_gradient(low = "#FFFFFF", high = "#E63F00") +
             theme(plot.background = element_rect(fill = "transparent", colour = NA),
                   panel.border = element_blank(), panel.grid = element_blank(),
                   axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                   legend.position = "right") +
             theme_classic(base_size = wordsize)
  gg
}

plot_tm<-function(data, wordsize, scale){
  gg <- ggplot(data, aes(x = time, y = sum_count)) + geom_bar(stat = "identity", position = "dodge") + 
               ylab('Number of protein sequence records') +
               scale_x_date(date_breaks = "1 months",labels = date_format("%b-%Y")) +
               theme_classic(base_size = wordsize) + 
               theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
                     axis.title.x = element_blank())
  if(scale == 'log'){
    gg <- gg + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)))
  }
  gg
}