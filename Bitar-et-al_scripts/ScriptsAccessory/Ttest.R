
#!/usr/bin/env Rscript
# module load R/4.0.2

# Use arguments:
args<-commandArgs(TRUE)

library("dplyr")
library("ggpubr")
library(PairedData)

t <- read.table(args[1], sep="\t", header=FALSE)

ref <- t$V1
other <- t$V2

my_data <- data.frame(group = rep(c(args[2], args[3]), each = 8), tsi = c(ref, other))

group_by(my_data, group) %>% summarise(count = n(), mean = mean(tsi, na.rm = TRUE), sd = sd(tsi, na.rm = TRUE))

tiff("TSIboxplots.tiff", units="in", width=5, height=5, res=300)
ggboxplot(my_data, x = "group", y = "tsi", color = "group", palette = c("#00AFBB", "#E7B800"), order = c(args[2], args[3]), ylab = "TSI", xlab = "Groups")
dev.off()

ref <- subset(my_data, group == args[2], tsi, drop = TRUE)
other <- subset(my_data, group == args[3], tsi, drop = TRUE)

tiff("TSIlineplots.tiff", units="in", width=5, height=5, res=300)
pd <- paired(ref, other)
plot(pd, type = "profile") + theme_bw()
dev.off()

# compute the difference
d <- with(my_data, tsi[group == args[2]] - tsi[group == args[3]])

# Shapiro-Wilk normality test for the differences
shapiro.test(d)

# If p-value is greater than the significance level 0.05 implying that the distribution of the differences (d) are not significantly different from normal distribution. In other words, we can assume the normality.

# Run t-test:
res1 <- t.test(tsi ~ group, data = my_data, paired = TRUE)
res2 <- t.test(tsi ~ group, data = my_data, paired = TRUE, alternative = "less")
res3 <- t.test(tsi ~ group, data = my_data, paired = TRUE, alternative = "greater")

res1
res2
res3

