library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)
library(deeptime)

## Permian 
load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0251.CO2.O2.1EFD ecotypes summary shelves short.no.poles.Rdata")

ecotypes.summary.short$eq.temp <-as.character(round(ecotypes.summary.short$eq.temp, digits=1))
ecotypes.summary.short$O2.pal <-as.numeric(paste(ecotypes.summary.short$O2.pal))

permian <- ggplot(ecotypes.summary.short, aes(y=ecotypes.number/10, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_smooth(aes( y=ecotypes.number/10, x=O2.pal, group=CO2.pal), color="grey40", size=1, linetype=2, se=FALSE, span=.7)+ 
  geom_point(shape=21, aes(fill=eq.temp), size=6, alpha=0.99)+theme_bw()+scale_fill_viridis(name = "Mean equatorial\ntemperature (°C)", discrete=T)+ 
  ylab("Metabolic habitat viability (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.1,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.35),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+coord_cartesian(xlim=c(0,1), ylim=c(0,100))+
  annotate(geom="text", x=.12, y=100.5, label="A) Permian", fontface="italic", size=10.3)
permian

## Paleocene
load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/p0055.CO2.O2.1EFD ecotypes summary shelves short.no.poles.Rdata")

ecotypes.summary.short$eq.temp <-as.character(round(ecotypes.summary.short$eq.temp, digits=1))
ecotypes.summary.short$O2.pal <-as.numeric(paste(ecotypes.summary.short$O2.pal))

paleocene <- ggplot(ecotypes.summary.short, aes(y=ecotypes.number/10, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_smooth(aes( y=ecotypes.number/10, x=O2.pal, group=CO2.pal), color="grey40", size=1, linetype=2, se=FALSE, span=.7)+ 
  geom_point(shape=21, aes(fill=eq.temp), size=6, alpha=0.99)+theme_bw()+scale_fill_viridis(name = "Mean equatorial\ntemperature (°C)", discrete=T)+ 
  ylab("Metabolic habitat viability (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.1,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.35),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+coord_cartesian(xlim=c(0,1), ylim=c(0,100))+
  annotate(geom="text", x=.15, y=100.5, label="B) Paleocene", fontface="italic", size=10.3)
paleocene

## Ordovician shallow e-folding 
load("/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/fm0450.CO2.O2.0.34EFD ecotypes summary shelves short.no.poles.Rdata")

ecotypes.summary.short$eq.temp <-as.character(round(ecotypes.summary.short$eq.temp, digits=1))
ecotypes.summary.short$O2.pal <-as.numeric(paste(ecotypes.summary.short$O2.pal))

shallow <- ggplot(ecotypes.summary.short, aes(y=ecotypes.number/10, x=O2.pal))+
  annotate(geom="rect", xmin=.1, xmax=.5, ymin=-Inf, ymax=Inf, fill="grey85")+
  annotate(geom="text", x=.3, y=0, label = expression("Early Paleozoic O"[2]), color="grey0", size=7)+
  geom_smooth(aes( y=ecotypes.number/10, x=O2.pal, group=CO2.pal), color="grey40", size=1, linetype=2, se=FALSE, span=.7)+ 
  geom_point(shape=21, aes(fill=eq.temp), size=6, alpha=0.99)+theme_bw()+scale_fill_viridis(name = "Mean equatorial\ntemperature (°C)", discrete=T)+ 
  ylab("Metabolic habitat viability (%)")+xlab(expression("Atmospheric O"[2]*" (PAL)"))+
  theme(plot.margin = margin(.1,1,1,1,"cm"),panel.border = element_rect(fill=NA,color="black", size=2,linetype="solid"),
        axis.line = element_line(lineend = 'square'), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=34),
        axis.text = element_text( size=26, color="black"),
        legend.title = element_text(size=24),
        legend.text = element_text( size=18),
        axis.ticks.length = unit(5, "points"),
        legend.justification=c(1,1), legend.position=c(.98,.35),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+coord_cartesian(xlim=c(0,1), ylim=c(0,100))+
  annotate(geom="text", x=.48, y=100.5, label="C) Ordovician: Shallow remineralization", fontface="italic", size=10.3)
shallow

plots <- ggarrange2(permian, paleocene, shallow, ncol=2)

ggsave("Figures/Viability figure supplementary ensembles shelves (Figure S2).pdf", plots, width=19.3, height=18, units="in")

