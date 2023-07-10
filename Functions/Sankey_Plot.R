library(ggalluvial)
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggsci)

# Uses celltype count as input 
data_bar = read_excel("~/Manuscript/Data/Supplementary/Cell_Counts_Sankey.xlsx",
                      sheet = 1,na="") 
data_bar$Celltype <- factor(data_bar$`Cell Types`)
table(data_bar$Celltype)
summary(data_bar)

#baseline
data_bar_baseline <- subset(data_bar, select = c("Celltype", "Baseline"))
summary(data_bar_baseline)
data_bar_baseline = data_bar_baseline %>% 
  group_by(Celltype) %>% 
  ungroup() %>% 
  mutate(data_bar_baseline, percentage = round(Baseline/sum(Baseline),3)) %>%
  arrange(factor(data_bar_baseline$Celltype))
table(data_bar_baseline$percentage)
sum(data_bar_baseline$percentage)
Baseline_pecentage <- data_bar_baseline$percentage

#5weeks
data_bar_5weeks <- subset(data_bar, select = c("Celltype", "5 Weeks"))
summary(data_bar_5weeks)
data_bar_5weeks = data_bar_5weeks %>% 
  group_by(Celltype) %>% 
  ungroup() %>% 
  mutate(data_bar_5weeks, percentage = round(`5 Weeks`/sum(`5 Weeks`),3)) %>%
  arrange(factor(data_bar_5weeks$Celltype))
table(data_bar_5weeks$percentage)
sum(data_bar_5weeks$percentage)
W5_pecentage <- data_bar_5weeks$percentage

#8weeks
data_bar_8weeks <- subset(data_bar, select = c("Celltype", "8 Weeks"))
summary(data_bar_8weeks)
data_bar_8weeks = data_bar_8weeks %>% 
  group_by(Celltype) %>% 
  ungroup()%>% 
  mutate(data_bar_8weeks, percentage = round(`8 Weeks`/sum(`8 Weeks`),3)) %>%
  arrange(factor(data_bar_8weeks$Celltype), 
          desc(data_bar_8weeks))
table(data_bar_8weeks$percentage)
sum(data_bar_8weeks$percentage)
W8_pecentage <- data_bar_8weeks$percentage

df=data.frame(
  Cell_type= levels(data_bar$Celltype),
  Baseline=Baseline_pecentage,
  "5 weeks"= W5_pecentage,
  "8 weeks"= W8_pecentage 
)

summary(df)

### Plot Sankey ###

a <- melt(df)
summary(a)

a$subject <- c(1:9,1:9,1:9)

a <- transform(a, variable = factor(variable, levels = c("Baseline","X5.weeks","X8.weeks")))

summary(a)


ggplot(a, aes(x = variable, stratum = Cell_type, alluvium = subject, y = value , 
            fill = Cell_type,label = Cell_type, position = "fill")) +
   scale_x_discrete(expand = c(.1, .1)) +
   geom_flow(alpha = 0.5, width = 1/3,) +
   theme_classic() + scale_fill_npg() +
   geom_stratum(alpha = 1, width = 1/3) +
   #geom_text(stat = "stratum", size = 3) +
   theme(legend.position = "right") + ylab("Ratio") + xlab("") 

ggplot(a,aes(x = variable, stratum = Cell_type, alluvium = subject,y = value ,
             fill = Cell_type, label = Cell_type)) +
  scale_x_discrete(expand = c(.05, .05)) +
  geom_flow( alpha = 1, width = 0 ) + # curve_type = "linear" "quintic" "cubic" "sine" "sigmoid" "arctangent"
  theme_classic() + scale_fill_npg() +
  theme(legend.position = "right") + ylab("Ratio") + xlab("") 

# ggsave("Sankey plot1.pdf",device = "pdf",
#        path = "./",
#        width = 4, height = 4, units = "in")