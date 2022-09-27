library(MultiAssayExperiment)
library(SingleCellExperiment)
library(tidyverse)
library(openxlsx)
library(VennDiagram)

load("../cohort_2016_17/ALFNA16_MultiAssayExperiment.RData")
load(file = "../cohort_2016_17/data/SCE_QC_pass_finalised.RData")
set2 <- read.xlsx(xlsxFile = "EJCannotated_FcRL5_SampleIDs.xlsx", sheet = 2)


## Aurora cohort ##
aurora.PIDs <- c("518M","523S", "546S", "642X", "651G","526W",
                 "544Q","559G", "630J", "637R", "557E", "602D",
                 "622A","627F","643Y","652H", "501T","507A","508B")


summary(aurora.PIDs %in% sce$PID)
intersect(aurora.PIDs, sce$PID)

length(intersect(set2$`Alyssa.ext.set.(day.0.+.7)`, sce$PID))

length(intersect(set2$`Alyssa.ext.set.(day.0.+.7)`, aurora.PIDs))

venn.diagram(
  x = list(levels(factor(sce$PID)), aurora.PIDs, set2$`Alyssa.ext.set.(day.0.+.7)`),
  category.names = c("scRNAseq","Aurora","Aurora mk2"),
  filename="test.png",
  output=T
  )


## Get HAI data for Aurora and SCE ##
hai.data <- wideFormat(ALFNA16[ , , 
                                # select assays with names that match:
                                # 'HAI'
                                # which returns all data from these two assays at all available days.
                                grepl(names(assays(ALFNA16)), pattern = "HAI") ],
                       # Return ages as well:
                       colDataCols = "Age.Group")


hai.data <- as.data.frame(hai.data)

aurora.hai.data <- hai.data %>% filter(primary %in% aurora.PIDs)
auroraMk2.hai.data <- hai.data %>% filter(primary %in% set2$`Alyssa.ext.set.(day.0.+.7)`)
sce.hai.data <- hai.data %>% filter(primary %in% sce$PID)
replication.hai.data <- hai.data %>% filter(primary %in% sce$PID) %>%
  filter(primary %in% aurora.PIDs | primary %in% set2$`Alyssa.ext.set.(day.0.+.7)`)

## Get ELISPOT data for Aurora and SCE ##
elispot.data <- wideFormat(ALFNA16[ , , 
                                # select assays with names that match:
                                # 'HAI'
                                # which returns all data from these two assays at all available days.
                                grepl(names(assays(ALFNA16)), pattern = "ELISPOT") ],
                       # Return ages as well:
                       colDataCols = "Age.Group")


elispot.data <- as.data.frame(elispot.data)


elispot.data$Age <- factor(c("22-36yo", "67-86yo")[elispot.data$Age.Group])

aurora.elispot.data <- elispot.data %>% filter(primary %in% aurora.PIDs)
auroraMk2.elispot.data <- elispot.data %>% filter(primary %in% set2$`Alyssa.ext.set.(day.0.+.7)`)
sce.elispot.data <- elispot.data %>% filter(primary %in% sce$PID)
replication.elispot.data <- elispot.data %>% filter(primary %in% sce$PID) %>%
  filter(primary %in% aurora.PIDs | primary %in% set2$`Alyssa.ext.set.(day.0.+.7)`)


## Flag sharing ##

aurora.elispot.data <- aurora.elispot.data %>%
  mutate(shared = ifelse(primary %in% sce.elispot.data$primary, "yes", "no"))

sce.elispot.data <- sce.elispot.data %>%
  mutate(shared = ifelse(primary %in% aurora.elispot.data$primary, "yes", "no"))


aurora.hai.data <- aurora.hai.data %>%
  mutate(shared = ifelse(primary %in% sce.hai.data$primary, "yes", "no"))

sce.hai.data <- sce.hai.data %>%
  mutate(shared = ifelse(primary %in% aurora.hai.data$primary, "yes", "no"))

## Plots ##
## HAI first ##
cowplot::plot_grid(
ggplot(data = 
         sce.hai.data %>%
         dplyr::select(primary, Age.Group, shared,
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = variable, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL, group=NULL), col="grey", trim=F) +
  geom_line(position = position_dodge(width=0.4), col = "lightgrey", lwd=0.5) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="day", y="A/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=T, comparisons = list(c(1,2), c(1,3), c(2,3))) +
  scale_y_continuous(trans = "log2", limits = c(1, 15000)) +
  
  facet_grid(.~Age) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="scRNAseq"),

ggplot(data = 
         sce.hai.data %>%
         dplyr::select(primary, Age.Group, 
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age)),
       aes(x = Age, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL), col="grey", trim=F) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="age", y="\nA/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=F, comparisons=list(c(1,2))) +
  scale_y_continuous(trans = "log2", limits = c(1, 30000)) +
  
  facet_grid(.~variable) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="scRNAseq"),

## aurora  ##
ggplot(data = 
         aurora.hai.data %>%
         dplyr::select(primary, Age.Group, shared,
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = variable, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL, group=NULL), col="grey", trim=F) +
  geom_line(position = position_dodge(width=0.4), col = "lightgrey", lwd=0.5) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="day", y="A/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=T, comparisons = list(c(1,2), c(1,3), c(2,3))) +
  scale_y_continuous(trans = "log2", limits = c(1, 15000)) +
  
  facet_grid(.~Age) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="aurora"),

ggplot(data = 
         aurora.hai.data %>%
         dplyr::select(primary, Age.Group, shared,
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = Age, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL), col="grey", trim=F) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="age", y="\nA/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=F, comparisons=list(c(1,2))) +
  scale_y_continuous(trans = "log2", limits = c(1, 30000)) +
  
  facet_grid(.~variable) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="aurora"),


## aurora  ##
ggplot(data = 
         auroraMk2.hai.data %>%
         dplyr::select(primary, Age.Group, 
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = variable, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL, group=NULL), col="grey", trim=F) +
  geom_line(position = position_dodge(width=0.4), col = "lightgrey", lwd=0.5) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="day", y="A/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=T, comparisons = list(c(1,2), c(1,3), c(2,3))) +
  scale_y_continuous(trans = "log2", limits = c(1, 15000)) +
  
  facet_grid(.~Age) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="auroraMk2"),

ggplot(data = 
         auroraMk2.hai.data %>%
         dplyr::select(primary, Age.Group,
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = Age, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL), col="grey", trim=F) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="age", y="\nA/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=F, comparisons=list(c(1,2))) +
  scale_y_continuous(trans = "log2", limits = c(1, 30000)) +
  
  facet_grid(.~variable) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="auroraMk2"),



## aurora replication data (i.e. PIDs NOT in SCE) ##
ggplot(data = 
         replication.hai.data %>%
         dplyr::select(primary, Age.Group, 
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = variable, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL, group=NULL), col="grey", trim=F) +
  geom_line(position = position_dodge(width=0.4), col = "lightgrey", lwd=0.5) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="day", y="A/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=T, comparisons = list(c(1,2), c(1,3), c(2,3))) +
  scale_y_continuous(trans = "log2", limits = c(1, 15000)) +
  
  facet_grid(.~Age) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="aurora replication"),

ggplot(data = 
         replication.hai.data %>%
         dplyr::select(primary, Age.Group,
                       contains("HAI_d"), -contains("log2")) %>%
         mutate(Age.Group=str_replace_all(Age.Group, c("young"="18-22yo", "old"="67-86yo"))) %>%
         mutate("Age" = factor(paste0(Age.Group))) %>%
         dplyr::select(-Age.Group) %>%
         set_colnames(
           gsub(colnames(.), 
                pattern = "HAI_d|_titre", replacement = "")) %>%
         reshape2::melt( id_vars = c(primary, Age, shared)),
       aes(x = Age, y=value,fill = Age, group=primary)) +
  geom_violin(aes(fill=NULL), col="grey", trim=F) +
  geom_point(shape =21, position = position_dodge(width=0.4)) +
  scale_fill_manual(values = c("white", "darkgrey")) +
  labs(x ="age", y="\nA/Cal09 HAI titer") +
  stat_compare_means(method="wilcox", paired=F, comparisons=list(c(1,2))) +
  scale_y_continuous(trans = "log2", limits = c(1, 30000)) +
  
  facet_grid(.~variable) +
  theme_pubr(legend = "none") +
  theme(strip.background = element_rect(fill = "NA")) +
  labs(title="aurora replication"),


##ELISPOTs

ggboxplot(sce.elispot.data, x = "Age", y = "ELISPOT_d7_HA..per.1e6.PBMC",
          ylab = "HA ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,600) +
  xlab(NULL)  + rotate_x_text()  +
  labs(title="scRNAseq"),
ggboxplot(sce.elispot.data, x = "Age", y = "ELISPOT_d7_Fluvax..per.1e6.PBMC",
          ylab = "Fluvax ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,3000) +
  xlab(NULL) + rotate_x_text() +
  labs(title="scRNAseq"),

ggboxplot(aurora.elispot.data, x = "Age", y = "ELISPOT_d7_HA..per.1e6.PBMC",
          ylab = "HA ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,600) +
  xlab(NULL)  + rotate_x_text() +
  labs(title="aurora"),
ggboxplot(aurora.elispot.data, x = "Age", y = "ELISPOT_d7_Fluvax..per.1e6.PBMC",
          ylab = "Fluvax ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,3000) +
  xlab(NULL) + rotate_x_text() +
  labs(title="aurora"),

ggboxplot(aurora.elispot.data, x = "Age", y = "ELISPOT_d7_HA..per.1e6.PBMC",
          ylab = "HA ELISPOT d7\nper 10^6 PBMC", facet.by = "shared") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,600) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="aurora, shared with scRNAseq"),

ggboxplot(aurora.elispot.data, x = "Age", y = "ELISPOT_d7_Fluvax..per.1e6.PBMC",
          ylab = "Fluvax ELISPOT d7\nper 10^6 PBMC", facet.by = "shared") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,3000) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="aurora, shared with scRNAseq"),

ggboxplot(auroraMk2.elispot.data, x = "Age", y = "ELISPOT_d7_HA..per.1e6.PBMC",
          ylab = "HA ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,600) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="auroraMk2"),

ggboxplot(auroraMk2.elispot.data, x = "Age", y = "ELISPOT_d7_Fluvax..per.1e6.PBMC",
          ylab = "Fluvax ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,3000) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="auroraMk2"),


ggboxplot(replication.elispot.data, x = "Age", y = "ELISPOT_d7_HA..per.1e6.PBMC",
          ylab = "HA ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,600) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="aurora, replication"),

ggboxplot(replication.elispot.data, x = "Age", y = "ELISPOT_d7_Fluvax..per.1e6.PBMC",
          ylab = "Fluvax ELISPOT d7\nper 10^6 PBMC") + stat_compare_means(method ="wilcox", paired=F, comparisons = list(c(1,2))) + 
  ylim(0,3000) +
  xlab(NULL)  + rotate_x_text() + 
  labs(title="aurora, replication"),



ncol=4)
