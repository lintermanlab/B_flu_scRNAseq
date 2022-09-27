
study.layout <-
  cbind("Assay"=c(
    "Spectral cytometry",
    "Spectral cytometry",
    "Spectral cytometry"),
    
    "day"=as.numeric(c(0,7,14)))

study.layout <- as.data.frame(study.layout)
study.layout$day <- as.numeric(as.character(study.layout$day))
study.layout$Assay <- factor(study.layout$Assay)

study <- ggplot(study.layout, aes(x=day, y=Assay, col = Assay,fill=Assay)) + geom_point(shape=25, size=3) + 
  scale_color_manual(values="black") +
  scale_fill_manual(values="black") +
  scale_x_continuous(breaks=c(0,7,14), labels=c(0,7,42))+
  theme_bw() + theme(legend.position = "none", 
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     text = element_text(size=12),
                     axis.title.y = element_blank()) +
  annotate(geom= "segment",x = 0,yend = 1.2,xend = 0, y= 1.5, size = 1, arrow = arrow(length = unit(6, "points"))) +
  annotate(geom="text", x=1, y=1.35, label="Trivalent flu vaccine", hjust="left") +
  labs(subtitle = "Healthy volunteers\n18-36 years old (n=11)\n66-89 years old (n=8)")

study
ggsave(filename = "Figure_3A_design.svg", width=3,height=2,dpi=300)