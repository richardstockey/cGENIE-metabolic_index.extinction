library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)
library(deeptime)
library(pals)
library(RColorBrewer)

###################
## Ordovician
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$cells.y[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 0

ecotypes.summary.baseCO2vshotCO2.long <- as.numeric()
for (i in 1:20){
  
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y < i] <- 1
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y >= i] <- 0
  
  ecotypes.summary.baseCO2vshotCO2$cell.threshold <- i
  
  ecotypes.summary.baseCO2vshotCO2.long <- rbind(ecotypes.summary.baseCO2vshotCO2.long, ecotypes.summary.baseCO2vshotCO2)
  print(i)
}

plot.summary <- ecotypes.summary.baseCO2vshotCO2.long %>%
  group_by(O2.pal, iteration, cell.threshold) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal, cell.threshold) %>%
  summarise(total.extinction.mean = mean(total.extinction))

plot.summary$O2.pal <- as.numeric(paste(plot.summary$O2.pal))
plot.summary$cell.threshold <- as.numeric(paste(plot.summary$cell.threshold))

Ord <- ggplot(plot.summary, aes(y=total.extinction.mean*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_line(aes(color=cell.threshold, group=cell.threshold), size=1, linetype=1)+
  theme_bw()+  
  scale_color_gradientn(colours=parula(100)[0:90], guide = "colourbar", name = "Cell\nThreshold")+
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
  annotate(geom="text", x=.155, y=99, label="A) Ordovician", fontface="italic", size=10.3)


###################
## Permian
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles.Rdata")

CO2.base.pal <- 8

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$cells.y[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 0

ecotypes.summary.baseCO2vshotCO2.long <- as.numeric()
for (i in 1:20){
  
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y < i] <- 1
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y >= i] <- 0
  
  ecotypes.summary.baseCO2vshotCO2$cell.threshold <- i
  
  ecotypes.summary.baseCO2vshotCO2.long <- rbind(ecotypes.summary.baseCO2vshotCO2.long, ecotypes.summary.baseCO2vshotCO2)
  print(i)
}

plot.summary <- ecotypes.summary.baseCO2vshotCO2.long %>%
  group_by(O2.pal, iteration, cell.threshold) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal, cell.threshold) %>%
  summarise(total.extinction.mean = mean(total.extinction))

plot.summary$O2.pal <- as.numeric(paste(plot.summary$O2.pal))
plot.summary$cell.threshold <- as.numeric(paste(plot.summary$cell.threshold))

Permian <- ggplot(plot.summary, aes(y=total.extinction.mean*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_line(aes(color=cell.threshold, group=cell.threshold), size=1, linetype=1)+
  theme_bw()+  
  scale_color_gradientn(colours=parula(100)[0:90], guide = "colourbar", name = "Cell\nThreshold")+
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
  annotate(geom="text", x=.12, y=99, label="B) Permian", fontface="italic", size=10.3)


###################
## Paleocene
###################

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles.Rdata")

CO2.base.pal <- 1

ecotypes.summary.baseCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal)

ecotypes.summary.hotCO2 <- filter(ecotypes.summary, CO2.pal == CO2.base.pal*4)

ecotypes.summary.baseCO2vshotCO2 <- merge(ecotypes.summary.baseCO2, 
                                          ecotypes.summary.hotCO2, 
                                          by=c("ecotype", "O2.pal", "iteration"), all.x = T)

ecotypes.summary.baseCO2vshotCO2$cells.y[is.na(ecotypes.summary.baseCO2vshotCO2$cells.y) == T] <- 0

ecotypes.summary.baseCO2vshotCO2.long <- as.numeric()
for (i in 1:20){
  
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y < i] <- 1
  ecotypes.summary.baseCO2vshotCO2$extinction[ecotypes.summary.baseCO2vshotCO2$cells.y >= i] <- 0
  
  ecotypes.summary.baseCO2vshotCO2$cell.threshold <- i
  
  ecotypes.summary.baseCO2vshotCO2.long <- rbind(ecotypes.summary.baseCO2vshotCO2.long, ecotypes.summary.baseCO2vshotCO2)
  print(i)
}

plot.summary <- ecotypes.summary.baseCO2vshotCO2.long %>%
  group_by(O2.pal, iteration, cell.threshold) %>%
  summarise(total.extinction = mean(extinction)) %>%
  group_by(O2.pal, cell.threshold) %>%
  summarise(total.extinction.mean = mean(total.extinction))

plot.summary$O2.pal <- as.numeric(paste(plot.summary$O2.pal))
plot.summary$cell.threshold <- as.numeric(paste(plot.summary$cell.threshold))

Pal <- ggplot(plot.summary, aes(y=total.extinction.mean*100, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_line(aes(color=cell.threshold, group=cell.threshold), size=1, linetype=1)+
  theme_bw()+  
  scale_color_gradientn(colours=parula(100)[0:90], guide = "colourbar", name = "Cell\nThreshold")+
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
  annotate(geom="text", x=.15, y=99, label="C) Paleocene", fontface="italic", size=10.3)


###################
## Combine and plot
###################

sum <- ggarrange2(Ord, Permian, Pal, ncol=2)

ggsave("Figures/Extinction figure cell threshold 1EFD shelves (Figure S7).pdf", sum, width=9.8*2, height=9*2, units="in")

