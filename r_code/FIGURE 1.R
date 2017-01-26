### LOAD LIBRARIES
library('ggplot2')
library('dplyr')
library('reshape2')



### FIGURE 1: STEM DENSITY 

### LOAD DATA
liana<-read.csv("data/stem_density_woodystems.csv",header = TRUE)
str(liana)

### STEMS PER METRE SQUARE
liana_summary <- as.data.frame(liana %>%
      group_by(logging.site) %>%
      summarise(liana.sprout=sum(liana.sprout.origin)/440,liana.seed=sum(liana.seed.origin)/440,free.standing.stems=sum(free.standing.woody.plants)/440))

### CONVERT TO LONG FORMAT
liana_summ_long <- melt(liana_summary,id.vars = "logging.site", variable.name = "regen.type",value.name="abund.m2")

### SET ORDER OF BARS ON GRAPH
liana_summ_long$logging.site<-factor(liana_summ_long$logging.site, levels=c("felling gap","secondary skid trail","primary skid trail","log landing")) # disturbance category
liana_summ_long$regen.type<-factor(liana_summ_long$regen.type,levels=c("free.standing.stems","liana.seed","liana.sprout")) #regen type

ggplot(liana_summ_long,aes(x=logging.site,y=abund.m2,fill=regen.type)) +
      geom_bar (stat="identity") +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      ylab(expression(paste("stems"," - ",m^2))) +
      scale_fill_grey() +
      scale_x_discrete(" ") +
      scale_y_continuous(limits=c(0,4.5)) +
      theme(axis.title.x = element_text(size=16)) +
      theme (axis.title.y=element_text(size=16)) +
      theme (axis.text.y  = element_text(size=14)) +
      theme (axis.text.x  = element_text(size=14)) +
      theme(legend.title=element_blank()) +
      theme(legend.text=element_text(size=14))



