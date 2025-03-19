# The data in leaves.csv and roots.csv is in a basic .csv format and needs to
# be converted to long format for plotting/other actions.
# We also need to rename it to make it more UTI relevant

# I need to get your 327 MCQ done today, and for the workshop material I’m going to pare it down:
#  
# 1.⁠ ⁠set up the background of the experiment
# 2.⁠ ⁠⁠give them some “identified this gene, what is it, what does it do?” (BLAST search, pubmed search, that sort of thing - reading files competently )
# 3.⁠ ⁠⁠more background: made a knockout and complement of one of the genes (simple case)
# 4.⁠ ⁠⁠here’s the data, plot it in webR using this code
# 5.⁠ ⁠⁠oh-oh! here’s the actual experimental data, which is riddled with batch effects; get them to plot the data, which looks shite
# 6.⁠ ⁠⁠talk a wee bit about regression modelling/accounting for such effects
# 7.⁠ ⁠⁠get them to run the appropriate regression code
# 8.⁠ ⁠⁠have them plot the results - aren’t they much better


library(tidyverse)
library(lme4)

# (1) Clean leaves.csv
# Load data and anonymise columns
data1 <- readr::read_csv("./leaves.csv", col_types = "dddd") 
colnames(data1) <- c("WT", "KO", "empty", "complement", "batch")

# Add batch number and pivot to long format
data1 <- data1 %>%
#  mutate(complement = complement * 1) %>%  # cheat to make complement work
#  mutate(KO = KO * 0.8) %>%
  pivot_longer(cols=c("WT", "KO", "empty", "complement"),
               names_to="label",
               values_to="CFU") %>%
  arrange(label) %>%
  mutate(batch=as.integer(c(rep(1:2, each=5, times=2), rep(3:4, each=5, times=2)))) %>%
  mutate(logCFU=log10(CFU))

# Write long-form data
# We pretend that this data arises from attachment to catheter tubing
write_csv(data1, "catheter.csv")

# (2) Clean roots.csv
# Load data and anonymise columns
data2 <- readr::read_csv("./roots.csv", col_types = "dddd") 
colnames(data2) <- c("WT", "KO", "empty", "complement", "batch")

# Add batch number and pivot to long format
data2 <- data2 %>%
#  mutate(KO = KO + 0.5) %>%
  pivot_longer(cols=c("WT", "KO", "empty", "complement"),
               names_to="label",
               values_to="CFU") %>%
  arrange(label) %>%
  mutate(batch=as.integer(c(rep(1:4, each=5, times=2), rep(5:8, each=5, times=2)))) %>%
  mutate(logCFU=log10(CFU))

# Write long-form data
# We pretend that this data arises from attachment to human tissue culture
# in microfluidics
write_csv(data2, "tissue.csv")

# Test plots of each dataset
p1 <- ggplot(data1, aes(x=factor(label, level=c("WT", "KO", "empty", "complement")),
                        y=logCFU)) +
  geom_boxplot() +
  geom_jitter(width=0.2, aes(colour=factor(batch))) +
  xlab("experiment")
p2 <- ggplot(data2, aes(x=factor(label, level=c("WT", "KO", "empty", "complement")),
                        y=logCFU)) +
  geom_boxplot() +
  geom_jitter(width=0.2, aes(colour=factor(batch))) +
  xlab("experiment")

p1
p2

# Test fit of model to account for batch effects, etc.
data1.mod <- data1 %>%
  mutate(KO = as.integer(label == "KO")) %>%
  mutate(empty = as.integer(label == "empty")) %>%
  mutate(complement = as.integer(label == "complement")) %>%
  mutate(label=factor(label))
  
write_csv(data1.mod, "catheter.csv")

model1a <- lm(logCFU ~ KO + empty + complement, data=data1.mod)
model1b <- lmer(logCFU ~ KO + empty + complement + (1 | batch), data=data1.mod)

data2.mod <- data2 %>%
  mutate(KO = as.integer(label == "KO")) %>%
  mutate(empty = as.integer(label == "empty")) %>%
  mutate(complement = as.integer(label == "complement")) %>%
  mutate(label=factor(label))

write_csv(data2.mod, "tissue.csv")

model2a <- lm(logCFU ~ KO + empty + complement, data=data2.mod)
model2b <- lmer(logCFU ~ KO + empty + complement + (1 | batch), data=data2.mod)