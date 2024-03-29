library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)

# Note that the plot made using this script was slightly edited for publication using Adobe Illustrator in 
# order to relevel overlapping ribbons. 

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

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles.Rdata")

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

load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves.ecotype.sensitivity.no.poles.Rdata")

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

ggplot(total.extinction.sum, aes(ymin=total.extinction.05*100, ymax=total.extinction.95*100, x=O2.pal*100))+
  annotate(geom="rect", xmin=.1*100, xmax=.5*100, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3*100, y=-1, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_ribbon(aes(group=Tectonics), alpha=1, fill="grey100")+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics),  size=.1, color="grey70",  alpha=.4)+ 
  geom_ribbon(aes(group=Tectonics, fill=Tectonics, ymin=total.extinction.25*100, ymax=total.extinction.75*100),   alpha=1)+ 
  theme_bw()+  
  scale_fill_manual(name = "Continental\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  scale_color_manual(name = "Continental\nconfiguration", values=c( rgb(153,192,141, maxColorValue = 255), rgb(254,191,101, maxColorValue = 255), rgb(239,88,69, maxColorValue = 255)), labels=c( "Ordovician", "Paleocene", "Permian"))+ 
  ylab("Global ecophysiotype extinction (%)")+xlab(expression("Atmospheric O"[2]*" (% PAL)"))+
  theme(plot.margin = margin(.5,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.98),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+coord_cartesian(xlim=c(0,100), ylim=c(0,100))+
  guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1))

ggsave("Figures/Extinction figure continental configs 1EFD shelves (Figure 5).pdf", width=9.8, height=9, units="in")

