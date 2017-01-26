### LOAD LIBRARIES
library('ggplot2')
library('dplyr')
library('reshape2')

### FIGURE 2: FAMILY ABUNDANCE

### LOAD DATA
liana_fam<-read.csv("data/family_data_final.csv",header = TRUE)
str(liana_fam)


### CONVERY NA'S TO ZEROS
liana_fam[is.na(liana_fam)] <- 0

### DROP REPLICATE ID AND MELT TO LONG FORMAT
liana_fam2 <- liana_fam %>%
      select(-c(replicate.id))
liana_fam_long <- melt(liana_fam2,id.vars = c("logging.site","regen.type"),variable.name="bot.family",value.name = "abund")

liana_fam_sum <- liana_fam_long %>%
      group_by(logging.site,regen.type,bot.family) %>%
      summarise(stem.count.m2=sum(abund/440))

### set levels
liana_fam_sum$logging.site <- factor(liana_fam_sum$logging.site, levels=c("felling gap","secondary skid trail","primary skid trail","log landing")) # logging sites

liana_fam_sum$regen.type <- factor(liana_fam_sum$regen.type, levels=c("Sprout","Seed")) # logging sites

bot.fam <- sort(as.character(unique(liana_fam_sum$bot.family)))
liana_fam_sum$bot.family <- factor(liana_fam_sum$bot.family, levels=bot.fam) # logging sites



ggplot(liana_fam_sum, aes(bot.family,(stem.count.m2),fill=regen.type)) + 
      geom_bar(stat="identity") +
      facet_wrap(~logging.site,ncol=2) +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      ylab(expression(paste("stems"," - ",m^2))) +
      scale_x_discrete(name=" ") +
      theme(axis.text.x=element_text(size=12,angle=90,vjust=0)) +
      theme(axis.text.y=element_text(size=14)) +
      theme(axis.title.y = element_text (size=16)) +
      theme(axis.title.x = element_text (size=16)) +
      theme(strip.text.x = element_text(size = 14)) +
      scale_fill_grey() +
      theme(legend.title=element_blank()) +
      theme(legend.text=element_text(size=14))

