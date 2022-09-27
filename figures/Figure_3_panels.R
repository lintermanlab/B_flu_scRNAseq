library(tidyverse)

flowCytoDat.HA_FCRL5 <- openxlsx::read.xlsx(xlsxFile = "../cohort_2016_17/data/FlowValidation/Fig2_FlowCytoValidation.xlsx")

FCRL5panel <- flowCytoDat.HA_FCRL5 %>%
  mutate(day = dplyr::case_when(Sample=="B1"~ "d0",
                                Sample=="B2"~ "d7",
                                Sample=="B3"~ "d42")) %>%
  mutate(day = factor(day, levels = c("d0","d7" , "d42"))) %>%
  mutate(Group  = factor(Age, levels = rev(c("Young", "Old")))) %>%
  mutate(Group = dplyr::case_when(Group=="Young"~ "22-36yo",
                                  Group=="Old"~"67-86yo", 
                                  is.na(Group)==TRUE ~"NA")) %>%
  mutate(Group  = factor(Group, levels = c("22-36yo", "67-86yo"))) %>%
  mutate(`%`= as.numeric(`HA+.FcRL5+.of.live`)) %>%
  ggplot(aes(x=day, y=`%`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=ID),
            col="lightgrey",position=position_dodge(width = 0.3)) +
  geom_point(aes(fill=Group,group=ID), shape=21,
             size=2,position=position_dodge(width = 0.3)) +
  
  scale_fill_manual(values = c("white", "darkgrey")) +
  scale_color_manual(values = c("black", "black")) + 
  ylim(0,0.15) +
  ggpubr::stat_compare_means(aes(group=ID), label.y = c(0.11, 0.13), comparisons = list(c(1,2),c(1,3)), method="wilcox", paired=T) + 
  facet_grid(.~Group) + 
  labs(y="HA+FCRL5+ B cells\n[% of live lymhocytes]") + 
  ggpubr::theme_pubr(legend = "none") +
  theme(axis.title.x = element_blank(),strip.background = element_rect(fill = "NA"))


FCRL5panelByDay <- flowCytoDat.HA_FCRL5 %>%
  mutate(day = dplyr::case_when(Sample=="B1"~ "d0",
                                Sample=="B2"~ "d7",
                                Sample=="B3"~ "d42")) %>%
  mutate(day = factor(day, levels = c("d0","d7" , "d42"))) %>%
  mutate(Group  = factor(Age, levels = rev(c("Young", "Old")))) %>%
  mutate(Group = dplyr::case_when(Group=="Young"~ "22-36yo",
                                  Group=="Old"~"67-86yo", 
                                  is.na(Group)==TRUE ~"NA")) %>%
  mutate(Group  = factor(Group, levels = c("22-36yo", "67-86yo"))) %>%
  mutate(`%`= as.numeric(`HA+.FcRL5+.of.live`)) %>%
  ggplot(aes(x=Group, y=`%`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=ID),
            col="lightgrey",position=position_dodge(width = 0.3)) +
  geom_point(aes(fill=Group,group=ID), shape=21,
             size=2,position=position_dodge(width = 0.3)) +
  
  scale_fill_manual(values = c("white", "darkgrey")) +
  scale_color_manual(values = c("black", "black")) + 
  ylim(0,0.15) +
  ggpubr::stat_compare_means(aes(group=ID),label.y=0.13, comparisons = list(c(1,2)), method="wilcox", paired=F) + 
  facet_grid(.~day) + 
  labs(y="HA+FCRL5+ B cells\n[% of live lymhocytes]") + 
  ggpubr::theme_pubr(legend = "none") +
  theme(axis.title.x = element_blank(),strip.background = element_rect(fill = "NA")) +
  ggpubr::rotate_x_text()



cowplot::plot_grid(FCRL5panel, FCRL5panelByDay)

ggsave(file = "Figure_3_HA_FCRL5_Bcells.svg", width=6, height=6, dpi = 300)
