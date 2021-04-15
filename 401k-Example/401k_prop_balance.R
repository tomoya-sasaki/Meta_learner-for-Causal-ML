
cov_select = c("age","educ","marr","inc","hown","twoearn")
d_overlap <- data.frame(data[,cov_select], D = as.factor(data$d), IPW = ifelse(data$d == 1, 1 / pseudo_all[,5], 1 / (1 - pseudo_all[,5])))
colnames(d_overlap) <- c("Age","Education","Married","Income","Homeowner","Two-earner","D","IPW")
plot.ipw <- melt(d_overlap,id.vars = c("D","IPW"))


# Define some colors
cbp <- c("#000000", "#E69F00", "#0072B2", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(plot.ipw, aes(x = value,weight=IPW, fill = D)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  facet_wrap( ~ variable, ncol = 2,scales = "free")+
  theme_cowplot() +
  scale_fill_manual(values = cbp)
