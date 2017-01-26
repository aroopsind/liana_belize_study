#########################################################################################
#########################################################################################

### LOAD LIBRARIES
###PACKAGES
library('ggplot2')
library('dplyr')
library('MASS')
library('reshape2')
library('BiodiversityR')



#########################################################################################
#########################################################################################

### LOAD DATA
liana<-read.csv("data/stem_density_woodystems.csv",header = TRUE)
str(liana)

### TOTAL LIANA ABUNDANCE [SPROUTS + SEEDS] - METRE SQUARE
liana$liana.total <- (liana$sprout.origin.lianas + liana$seed.origin.lianas)


#########################################################################################
#########################################################################################

### NEGATIVE BINOMIAL MODEL: INCLUDES OVERDISPERSION PARAMETER

liana.nb <- glm.nb(liana.total ~ logging.site -1,data=liana) #negative binomial model
summary(liana.nb) #model summary
est.coeff.nb <- cbind(Estimate=coef(liana.nb),confint(liana.nb)) # mean coef estimates + 95% CI
AIC(liana.nb)




#### PREDICTED VALUES FROM NEGATIVE BINOMIAL MODEL
pred.data <- data.frame(logging.site=factor(1:4,levels=1:4,labels=levels(liana$logging.site)))
alphahat_SE <- predict(liana.nb,pred.data,type="link",se.fit = TRUE)

pred.data2 <- cbind(pred.data,alphahat_SE$fit,alphahat_SE$se.fit)
colnames(pred.data2)[2:3] <- c("fit","se.fit")
pred.data.sum <- within(pred.data2,{
      pred.abund <- exp(fit)
      LL <- exp(fit-1.96 * se.fit)
      UL <- exp(fit+1.96 * se.fit)
})

pred.data.sum$logging.site <- factor(pred.data.sum$logging.site, levels=c("felling gap","secondary skid trail","primary skid trail","log landing")) # set levels
### PLOT OUT PREDICTIONS @ 440 M2
ggplot(pred.data.sum,aes(x=logging.site,y=pred.abund)) +
      geom_point(size=6) +
      theme_bw()+
      geom_segment(aes(x = logging.site, y = UL, xend = logging.site, yend = LL), size=3,data = pred.data.sum) +
      theme(axis.text.x=element_text(size=14)) +
      theme(axis.text.y=element_text(size=12)) +
      ylab (expression(paste("Predicted abundance ("," 440 ", m ^{2},")"))) +
      xlab(expression("")) +
      scale_y_continuous(limits=c(0,90)) 
      
      
      

### TEST ASSUMPTION OF OVERDISPERSION BY COMPARING WITH POSSION MODEL WITH NEGATIVE BINOMIAL FIT USING AIC
liana.pois <- glm(liana.total ~ logging.site -1,data = liana, family = "poisson")
summary(liana.pois)
AIC(liana.pois) # AIC value is twice as large as the negative binomial model



#########################################################################################
#########################################################################################

### DIVERSITY ANALYSIS
liana_data <- read.csv("data/liana_data_abund.csv") # liana abundance by botanical family and origin (seed/sprout)
str(liana_data)

###  NA'S= ZEROS
liana_data[is.na(liana_data)] <- 0

### CREATE UNIQUE LOG SITE + REP ID.
liana_data$site.rep <- paste(liana_data$logging.site,liana_data$replicate.id,sep=".")

### DROP UNNECESSARY COLUMNS
liana_data2 <- liana_data[,!(colnames(liana_data) %in% c("logging.site","regen.type","replicate.id"))]


### SUMMARISE ABUNDANCE BY FAMILY + LOG SITE [IGNORE REGEN TYPE]
liana_data_long <- melt(liana_data2,id.vars=c("site.rep"),variable.name = "bot.fam",value.name = "abund")

liana_data_sum <- liana_data_long %>%
      group_by(site.rep,bot.fam) %>%
      summarise(abund=sum(abund))


### CONVERT BACK TO WIDE FORMAT
liana_wide <- dcast(liana_data_sum, site.rep ~ bot.fam,value.var = "abund" )



### ADD LOG SITE ID AND REP ID
liana_wide$log.site <- sapply(strsplit(liana_wide$site.rep, split='.', fixed=TRUE), function(x) (x[1]))
liana_wide$rep.id <- sapply(strsplit(liana_wide$site.rep, split='.', fixed=TRUE), function(x) (x[2]))

### SET LEVELS OF LOGGING SITE

liana_wide$log.site <-  factor(liana_wide$log.site, levels =  c("felling gap","secondary skid trail","primary skid trail","log landing"))
levels(liana_wide$log.site)

### ABUNDANCE DATA
abund_data <- liana_wide[!colnames(liana_wide) %in% c("log.site","rep.id","site.rep","Unknown")]
abund_data_final <- as.matrix(abund_data);colnames(abund_data_final)<-NULL;rownames(abund_data_final)<-NULL

site.env <- liana_wide[colnames(liana_wide) %in% "log.site"]

##############################################################################
##############################################################################
### RENYI'S ENTROPHY

### compare diversity between subsets of the data
renyi <- renyicomp(abund_data,y=site.env,factor="log.site",permutations = 100)      


##exract out mean values from renyi entrophy
renyi.val.mean <- renyi[,,"mean"]


### plot out renyi values with ggplot
reny.long <- melt(renyi.val.mean)
### SET LEVELS OF LOGGING SITE
reny.long$log.site <- factor(reny.long$log.site, levels = c("felling gap","secondary skid trail","primary skid trail","log landing"))

ggplot(reny.long, aes(x=as.factor(scale), y=value, colour=log.site, group=log.site)) + 
      geom_line() + geom_point(size=3,shape=18) +
      xlab(expression(bold("Alpha"))) +
      ylab (expression(bold("Renyi diversity profile"))) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(axis.text.x=element_text(size=12)) +
      theme(axis.text.y=element_text(size=12)) +
      theme(axis.title =element_text(size=14)) +
      theme_bw() +
      theme(legend.title=element_blank())

##############################################################################
##############################################################################
### FIGURE 1: STEM DENSITY 

### LOAD DATA
liana<-read.csv("data/stem_density_woodystems.csv",header = TRUE)
str(liana)

### STEMS PER METRE SQUARE
liana_summary <- as.data.frame(liana %>%
                                     group_by(logging.site) %>%
                                     summarise(sprout.origin.lianas=sum(sprout.origin.lianas)/440,seed.origin.lianas=sum(seed.origin.lianas)/440,other.woody.plants=sum(other.woody.plants)/440))

### CONVERT TO LONG FORMAT
liana_summ_long <- melt(liana_summary,id.vars = "logging.site", variable.name = "regen.type",value.name="abund.m2")

### SET ORDER OF BARS ON GRAPH
liana_summ_long$logging.site<-factor(liana_summ_long$logging.site, levels=c("felling gap","secondary skid trail","primary skid trail","log landing")) # disturbance category
liana_summ_long$regen.type<-factor(liana_summ_long$regen.type,levels=c("other.woody.plants","seed.origin.lianas","sprout.origin.lianas")) #regen type

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





##############################################################################
##############################################################################
### FIGURE 2: FAMILY ABUNDANCE

### LOAD DATA
liana_fam<-read.csv("data/liana_data_abund.csv",header = TRUE)
str(liana_fam)


### CONVERY NA'S TO ZEROS
liana_fam[is.na(liana_fam)] <- 0

### DROP REPLICATE ID AND MELT TO LONG FORMAT
liana_fam2 <- liana_fam[-3]

liana_fam_long <- melt(liana_fam2,id.vars = c("logging.site","regen.type"),variable.name="bot.family",value.name = "abund")

liana_fam_sum <- liana_fam_long %>%
      group_by(logging.site,regen.type,bot.family) %>%
      summarise(stem.count.m2=sum(abund/440))

### set levels
liana_fam_sum$logging.site <- factor(liana_fam_sum$logging.site, levels=c("felling gap","secondary skid trail","primary skid trail","log landing")) # logging sites

liana_fam_sum$regen.type <- factor(liana_fam_sum$regen.type, levels=c("seed.origin.lianas","sprout.origin.lianas")) # logging sites

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
