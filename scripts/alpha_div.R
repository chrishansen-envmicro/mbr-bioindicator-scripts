#https://rpubs.com/lconteville/713954
#load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(phyloseq)
library(vegan)
library(DESeq2)
library(gridExtra)
library(GGally)
library(wesanderson)
library(ggpubr)
library(data.table)
library(cowplot)
library(mgcv)
library(SRS)
library(rstatix)

#import data from MSI: /panfs/roc/groups/6/hamil689/shared/sen/JM_second_paper/res/16SII
sharedfile <- "16S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "16S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated.csv"

#returns a phyloseq object with taxonomy and read counts where the TAXA are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
#import metadata
map <- read.csv(mapfile)
map$location <- recode(map$location, MH = "Mount Hood") %>% #rename location values to improve figure clarity
  recode(GNP = "Glacier National Park") 

#remove blanks
map.recharge <- filter(map, bio_rep != "blank")
map.recharge <- filter(map.recharge, ncbi_name_16s != "")
map.recharge <- sample_data(map.recharge)
rownames(map.recharge) <- map.recharge$ncbi_name_16s

#merge mothur data with metadata returns returns to create phyloseq object with otus, taxa and metadata are rows
moth_merge <- merge_phyloseq(mothur_data, map.recharge)

#replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#remove non-bacterial mitochondria and chloropast seq. should be removed during amp processing (mother), but this is a double-check
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

#after filtering to only include springs/streams, look at read depth
stream.spring <- c("Stream", "Spring")
data_spst <- subset_samples(data, sample_type %in% stream.spring)

#remove singletons
data_filt <- filter_taxa(data_spst, function (x) {sum(x > 0) > 1}, prune=TRUE)

#develop rarefaction curves and visualize distribution of read lengths: OPTIONAL
cols <-  c("#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F", "#660F45", "#707070",
           "#333333", "black", "#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F", "#660F45", "#707070",
           "#333333", "black", "#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F")

line <- c(rep(1,25), rep(4,25), rep(6,21))

otu_tab <- otu_table(data_filt)
class(otu_tab) <- "matrix" # as.matrix() will do nothing. you get a warning here, but this is what we need to have
otu_tab <- t(otu_tab) # transpose observations to rows
rare <- rarecurve(otu_tab, step=100, ylab="OTU",  label=F, col = cols, lty = line)
#develop rarefaction curves and visualize distribution of read lengths: OPTIONAL
#rarecurve(t(otu_table(data_filt)), step = 100, label = F)
abline(v=8000, lw = 1) #arbitrarily selected per rarecurve

#identify which samples had fewer than 10^3 reads and remove (12)
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (8000))
min(sample_sums(data_filt_clean)) #check min at 11000. revise above to be 90% of this min (11919)
sort(sample_sums(data_filt))

#how many samples remain from each mountain range?
sample_data(data_filt_clean) %>%
  group_by(location) %>%
  summarise(n = n())

#Rarefy the samples without replacement. Rarefaction is used to simulate even number 
#of reads per sample. In this example, the rarefaction depth chosen is the 90% of the 
#minimum sample depth in the dataset (in this case 459 reads per sample).
#rarefy based on lowers read count about 1k. (1015)
#https://micca.readthedocs.io/en/latest/phyloseq.html
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=0.9*min(sample_sums(data_filt_clean)), 
                                       replace=F)
#perform srs to standardize counts: https://rdrr.io/github/jfq3/QsRutils/man/srs_p.html
#data_filt_clean_srs <- srs_p(data_filt_clean)

#provides alpha diversity estimates for each sample
#estimate_richness(data_filt_clean_srs)
estimate_richness(data_filt_rarified)

a_my_comparisons <- list( c("Glacier National Park", "Mount Hood"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) #manually add wilcox values dervied below for aesthetics

#plot results from above estimated alpha diversity. isolating observed and shannon
p <- plot_richness(data_filt_rarified, x="location", 
                   measures=c("Observed", "Shannon"),
                   color = "location") 

p$layers <- p$layers[-1]
p +
  geom_jitter(aes(shape = location), width = 0.2, size = 5, alpha = 0.6, show.legend = FALSE)+
  geom_vline(xintercept = 0, color="snow4") +
  geom_hline(yintercept = 0, color="snow4") +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons,
                     label = "p.signif", symnum.args = symnum.args) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  stat_summary(fun.data = "mean_cl_boot",
               show.legend = FALSE,
               position = position_nudge(x = .3, y = 0),
               fun.args=(conf.int=0.95),
               color = "black",
               size = 0.7) +
  scale_x_discrete(element_blank()) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18))

#GET STATISTICS - test for normality -> observed richness data is near normal!
norm_check <- select(p$data, location, variable, value)
norm_check_shannon <- filter(norm_check, variable == "Shannon")
norm_check_obs <- filter(norm_check, variable == "Observed")
norm_check_shannon_mh <- filter(norm_check_shannon, location == "Mount Hood")
shapiro.test(norm_check_obs$value) #normally distributed
shapiro.test(norm_check_shannon$value) #not normally distributed

# R program to illustrate
# Bartlett’s test

# Using bartlett.test() -> MH and GNP DO NOT have differing variances
bartlett.test(value~location, norm_check_obs)

#derive summary statistics
norm_check %>%
  group_by(location, variable) %>%
  summarise(mean = mean(value))



#let's do a wilcox test: https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/
#except, ttest for obs
norm_check_shannon %>% 
  wilcox_test(value ~ location)

# shortcut for 95% confidence intervals
lm_rich <- lm(value ~ location, data=norm_check_obs)
confint(lm_rich, level=0.95)



#norm_check_obs_mh <- filter(norm_check_obs, location == "Mount Hood")
#norm_check_obs_gnp <- filter(norm_check_obs, location == "Glacier National Park")
#t.test(norm_check_obs_mh$value, norm_check_obs_gnp$value)

#should replicates be combined?
#extract potential biological replicates & either alpha diversity metric
reps_shan <- filter(p$data, variable == "Shannon" & bio_rep > 0)
reps_obs <- filter(p$data, variable == "Observed"  & bio_rep > 0)
#for each biological replicate per sample name, calculate the one-way ANOVA test: http://www.sthda.com/english/wiki/one-way-anova-test-in-r
#USE FOR SHANNON DIVERSITY
reps_shan_anova <- filter(reps_shan, sample_name == "Buzzard Fountain Spring" |
                            sample_name == "Choss Seep" |
                            sample_name == "Grinnell Spring" |
                            sample_name == "James Spring" |
                            sample_name == "Palmer B Spring" |
                            sample_name == "Paradise Park Spring" |
                            sample_name == "Piegan North Spring" |
                            sample_name == "Piegan Spring")
anov <- aov(value ~ sample_name, data = reps_shan_anova)
TukeyHSD(anov)

#USE FOR OBSERVED RICHNESS
reps_obs_anova <- filter(reps_obs, sample_name == "Buzzard Fountain Spring" |
                           sample_name == "Choss Seep" |
                           sample_name == "Grinnell Spring" |
                           sample_name == "James Spring" |
                           sample_name == "Palmer B Spring" |
                           sample_name == "Paradise Park Spring" |
                           sample_name == "Piegan North Spring" |
                           sample_name == "Piegan Spring")

anov <- summary(aov(value ~ sample_name, data = reps_obs_anova))
TukeyHSD(anov)

#I want to visualize how different replicates are
data_filt_rarified_test <- subset_samples(data_filt_rarified, sample_name == "Buzzard Fountain Spring" |
                                            sample_name == "Choss Seep" |
                                            sample_name == "Grinnell Spring" |
                                            sample_name == "James Spring" |
                                            sample_name == "Palmer B Spring" |
                                            sample_name == "Paradise Park Spring" |
                                            sample_name == "Piegan North Spring" |
                                            sample_name == "Piegan Spring")

plot_richness(data_filt_rarified_test, x="location", 
              measures=c("Shannon", "Observed"),
              color = "sample_name") +
  geom_jitter()


#create alpha diversity against continuous variable
palphadt <- data.table(p$data)
palphadt_shan <- palphadt[(variable == "Shannon")]
palphadt_obs <- palphadt[(variable == "Observed")]

# Order by Days - NOT NECESSARY(?)
palphadt_shan <- palphadt_shan[order(fi_best)][(is.na(se))]
palphadt_obs <- palphadt_obs[order(fi_best)][(is.na(se))]

library(ggtext)

# Define the plot
a <- ggplot(data = palphadt_shan, 
            mapping = aes(fi_best, value,
                          color = location, shape = location)) +
  # shape = ReportedAntibioticUsage)) + 
  geom_point(size = 5, alpha = 0.6, show.legend = FALSE) + 
  # geom_path() +
  #geom_point(data = alphadt[(ReportedAntibioticUsage == "Yes")], 
  # size = 8, alpha = 0.35) +
  facet_wrap(~location, ncol = 1) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  labs(y = "Shannon Index", x = "f<sub><i>i") +
  geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_markdown())

b <- ggplot(data = palphadt_obs, 
            mapping = aes(fi_best, value,
                          color = location, shape = location)) +
  # shape = ReportedAntibioticUsage)) + 
  geom_point(size = 5, alpha = 0.6, show.legend = FALSE) + 
  # geom_path() +
  #geom_point(data = alphadt[(ReportedAntibioticUsage == "Yes")], 
  # size = 8, alpha = 0.35) +
  facet_wrap(~location, ncol = 1) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  labs(y = "Observed Richness", x = "f<sub><i>i") +
  geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_markdown())

plot_grid(b, a, labels = "auto", label_size = 20)

#summarise correlation data. manually alter for either alpha diverstiy statistic
#FIRST separate by location
palphadt_shan_mh <- filter(palphadt_shan, location == "Mount Hood")
palphadt_shan_gnp <- filter(palphadt_shan, location == "Glacier National Park")
palphadt_obs_mh <- filter(palphadt_obs, location == "Mount Hood")
palphadt_obs_gnp <- filter(palphadt_obs, location == "Glacier National Park")

library(broom)

cor.test(x = palphadt_shan_mh %>% pull(fi_best),
         y = palphadt_shan_mh %>% pull(value)) %>%
  tidy()

cor.test(x = palphadt_shan_gnp %>% pull(fi_best),
         y = palphadt_shan_gnp %>% pull(value)) %>%
  tidy()

#Compare GNP and MH correlations via permutation
#Begin with SHANNON
obs_cors_shan <- palphadt_shan %>%
  group_by(location) %>%
  summarise(cor = cor(fi_best,value)) %>%
  summarise(diff_cor = diff(cor)) %>%
  pull()

perm_shan <- replicate(10000, simplify = FALSE,
                       expr = palphadt_shan %>% 
                         mutate(perm_location = sample(location, replace = FALSE)) %>%
                         group_by(perm_location) %>%
                         summarise(cor = cor(fi_best,value)) %>%
                         summarise(diff_cor = diff(cor))) %>%
  bind_rows()

summarise(perm_shan, pval = mean(abs(diff_cor) >= abs(obs_cors_shan)))

#Now OBSERVED
obs_cors_obs <- palphadt_obs %>%
  group_by(location) %>%
  summarise(cor = cor(fi_best,value)) %>%
  summarise(diff_cor = diff(cor)) %>%
  pull()

perm_shan <- replicate(10000, simplify = FALSE,
                       expr = palphadt_obs %>% 
                         mutate(perm_location = sample(location, replace = FALSE)) %>%
                         group_by(perm_location) %>%
                         summarise(cor = cor(fi_best,value)) %>%
                         summarise(diff_cor = diff(cor))) %>%
  bind_rows()

summarise(perm_shan, pval = mean(abs(diff_cor) >= abs(obs_cors_shan)))

