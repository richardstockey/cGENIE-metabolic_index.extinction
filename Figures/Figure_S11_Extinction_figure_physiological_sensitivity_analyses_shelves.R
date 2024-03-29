library(readr) 
library(ggplot2)
library(RNetCDF)
library(RColorBrewer)
library(reshape2)
library(viridis)
library(AquaEnv)
library(dplyr)
library(questionr)
library(ncdf4)


###################
## low Ao, high Eo
###################

###################
## Ordovician
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_highEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                                       ecotypes.summary.hotCO2, 
                                                       by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
                group_by(O2.pal, iteration) %>%
                summarise(total.extinction = mean(extinction)) %>%
                group_by(O2.pal) %>%
                summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
                          total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.ordovician <- plot.summary
plot.summary.ordovician$Tectonics <- "Ordovician"


###################
## Permian
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_highEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.permian <- plot.summary
plot.summary.permian$Tectonics <- "Permian"

###################
## Paleocene
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_highEo.Rdata")

CO2.base.pal <- 1

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.paleocene <- plot.summary
plot.summary.paleocene$Tectonics <- "Paleocene"

###################
## Combine and plot
###################

total.extinction.sum <- rbind(plot.summary.ordovician,
                              plot.summary.permian,
                              plot.summary.paleocene)
total.extinction.sum$O2.pal <- as.numeric(paste(total.extinction.sum$O2.pal))

low_high <- ggplot(total.extinction.sum, aes(ymin=total.extinction.05*100, ymax=total.extinction.95*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_ribbon(aes(group=Tectonics), alpha=1, fill="grey100")+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics),  size=.1, color="grey70",  alpha=.4)+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics, ymin=total.extinction.25*100, ymax=total.extinction.75*100),   alpha=1)+ 
  theme_bw()+  
  scale_fill_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  scale_color_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  ylab("Global ecophysiotype extinction (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.5,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlim(c(0,1))+ylim(c(0,100))+ 
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1))+
  annotate(geom="text", x=.23, y=99, label="A) low Ao, high Eo", fontface="italic", size=10.3)

###################
## high Ao, high Eo
###################

###################
## Ordovician
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_highEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.ordovician <- plot.summary
plot.summary.ordovician$Tectonics <- "Ordovician"


###################
## Permian
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_highEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.permian <- plot.summary
plot.summary.permian$Tectonics <- "Permian"

###################
## Paleocene
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_highEo.Rdata")

CO2.base.pal <- 1

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.paleocene <- plot.summary
plot.summary.paleocene$Tectonics <- "Paleocene"

###################
## Combine and plot
###################

total.extinction.sum <- rbind(plot.summary.ordovician,
                              plot.summary.permian,
                              plot.summary.paleocene)
total.extinction.sum$O2.pal <- as.numeric(paste(total.extinction.sum$O2.pal))

high_high <- ggplot(total.extinction.sum, aes(ymin=total.extinction.05*100, ymax=total.extinction.95*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_ribbon(aes(group=Tectonics), alpha=1, fill="grey100")+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics),  size=.1, color="grey70",  alpha=.4)+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics, ymin=total.extinction.25*100, ymax=total.extinction.75*100),   alpha=1)+ 
  theme_bw()+  
  scale_fill_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  scale_color_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  ylab("Global ecophysiotype extinction (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.5,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlim(c(0,1))+ylim(c(0,100))+ 
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1))+
  annotate(geom="text", x=.23, y=99, label="B) high Ao, high Eo", fontface="italic", size=10.3)


###################
## low Ao, low Eo
###################

###################
## Ordovician
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_lowEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.ordovician <- plot.summary
plot.summary.ordovician$Tectonics <- "Ordovician"


###################
## Permian
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_lowEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.permian <- plot.summary
plot.summary.permian$Tectonics <- "Permian"

###################
## Paleocene
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_lowAo_lowEo.Rdata")

CO2.base.pal <- 1

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.paleocene <- plot.summary
plot.summary.paleocene$Tectonics <- "Paleocene"

###################
## Combine and plot
###################

total.extinction.sum <- rbind(plot.summary.ordovician,
                              plot.summary.permian,
                              plot.summary.paleocene)
total.extinction.sum$O2.pal <- as.numeric(paste(total.extinction.sum$O2.pal))
low_low <- ggplot(total.extinction.sum, aes(ymin=total.extinction.05*100, ymax=total.extinction.95*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_ribbon(aes(group=Tectonics), alpha=1, fill="grey100")+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics),  size=.1, color="grey70",  alpha=.4)+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics, ymin=total.extinction.25*100, ymax=total.extinction.75*100),   alpha=1)+ 
  theme_bw()+  
  scale_fill_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  scale_color_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  ylab("Global ecophysiotype extinction (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.5,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlim(c(0,1))+ylim(c(0,100))+ 
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1))+
  annotate(geom="text", x=.23, y=99, label="C) low Ao, low Eo", fontface="italic", size=10.3)



###################
## high Ao, low Eo
###################

###################
## Ordovician
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_lowEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.ordovician <- plot.summary
plot.summary.ordovician$Tectonics <- "Ordovician"


###################
## Permian
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_lowEo.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.permian <- plot.summary
plot.summary.permian$Tectonics <- "Permian"

###################
## Paleocene
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles_highAo_lowEo.Rdata")

CO2.base.pal <- 1

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 1
ecotypes.summary.baseCO2vshotCO2$extinction[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == F] <- 0


plot.summary <- ecotypes.summary.baseCO2vshotCO2 %>%
  group_by(O2.pal, iteration) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal) %>%
  summarise(total.extinction.05 = quantile(total.extinction, 0.05), total.extinction.95 = quantile(total.extinction, 0.95), 
            total.extinction.25 = quantile(total.extinction, 0.25), total.extinction.75 = quantile(total.extinction, 0.75))

plot.summary.paleocene <- plot.summary
plot.summary.paleocene$Tectonics <- "Paleocene"

###################
## Combine and plot
###################

total.extinction.sum <- rbind(plot.summary.ordovician,
                              plot.summary.permian,
                              plot.summary.paleocene)
total.extinction.sum$O2.pal <- as.numeric(paste(total.extinction.sum$O2.pal))

high_low <- ggplot(total.extinction.sum, aes(ymin=total.extinction.05*100, ymax=total.extinction.95*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_ribbon(aes(group=Tectonics), alpha=1, fill="grey100")+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics),  size=.1, color="grey70",  alpha=.4)+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics, ymin=total.extinction.25*100, ymax=total.extinction.75*100),   alpha=1)+ 
  theme_bw()+  
  scale_fill_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  scale_color_manual(name = "Tectonic\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  ylab("Global ecophysiotype extinction (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.5,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+xlim(c(0,1))+ylim(c(0,100))+ 
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1))+
  annotate(geom="text", x=.23, y=99, label="D) high Ao, low Eo", fontface="italic", size=10.3)


plots <- ggarrange2(low_high, high_high, low_low,high_low, ncol=2)

ggsave("Figures/Extinction figure physiological sensitivity analyses shelves (Figure S11).pdf", plots, width=9.8*2, height=9*2, units="in")

