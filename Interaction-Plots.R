sampledata <- read.csv("C:/Users/pierc/Downloads/SHBG-testdata.csv")
library(dplyr)
library(ggplot2)
library(ggsignif)
library(stats)
library(rlang)
library(utils)
library(plotrix)
library(ggpubr)
library(rstatix)

#i - tablename
#j - SHBG
#k - rsxxxxxx
#x - exposure (Veg2)
#y - minor allele
#z - major allele
#a - number of minor alleles (Not yet implemented fully)
makeplots <- function(i,j,k,x,y,z)
{
  #Multiple Allele Prompting
  
  #if(a == 2)
  #{
  # userinput2 <- as_string(readline(prompt = "Minor Allele 2? "))
  #}
  
  #if(a == 3)
  #{
  #  userinput2 <- as_string(readline(prompt = "Minor Allele 2? "))
  #  userinput3 <- as_string(readline(prompt = "Minor Allele 3? "))
  #}
  
  #if(a == 4)
  #{
  #  userinput2 <- as_string(readline(prompt = "Minor Allele 2? "))
  #  userinput3 <- as_string(readline(prompt = "Minor Allele 3? "))
  #  userinput4 <- as_string(readline(prompt = "Minor Allele 4? "))
  #}
  
  #Put variables in quotes if they need them to play nicely with other functions
  i <- enquo(i)
  j <- enquo(j)
  k <- enquo(k)
  x <- enquo(x)
  y <- enquo(y)
  z <- enquo(z)
  
  
  #Creates a dataframe with mean and sd of SHBG, grouped by allele genotype
  testdataframe <- summarise( group_by(sampledata, !!k), mean_SHBG = mean(!!j), se_SHBG = std.error(!!j))
  
  #Split data by exposure and genotype then apply mean and sd functions
  exposure0 <- summarise( group_by(filter(sampledata, !!x == 0), !!k), mean_SHBG = mean(!!j), se_SHBG = std.error(!!j))
  exposure1 <- summarise( group_by(filter(sampledata, !!x == 1), !!k), mean_SHBG = mean(!!j), se_SHBG = std.error(!!j))
  exposure <- rbind(exposure0, exposure1)
  
  #First Graph Bound Calculations
  bound_lower <- min(testdataframe$mean_SHBG)-2*max(testdataframe$se_SHBG)
  bound_lower[bound_lower < 0] <- 0
  bound_upper <- max(testdataframe$mean_SHBG)+3*max(testdataframe$se_SHBG)
  bound_diff <- bound_upper - bound_lower

  signif01 <- max(testdataframe$mean_SHBG) + max(testdataframe$se_SHBG) + 0.09*bound_diff
  signif02 <- max(testdataframe$mean_SHBG) + max(testdataframe$se_SHBG) + 0.18*bound_diff
  signif03 <- max(testdataframe$mean_SHBG) + max(testdataframe$se_SHBG) + 0.27*bound_diff

  #SNP Allele Naming
  homo_minor <- paste(quo_name(y), quo_name(y), sep = "")
  het <- paste(quo_name(y), quo_name(z), sep = "")
  homo_major <- paste(quo_name(z), quo_name(z), sep = "")
  
  #2nd Plot Re-leveling
  sampledatachar <- mutate(sampledata, Genotype = factor(pull(sampledata, !!k), levels = c(0,1,2), labels = c(homo_major, het, homo_minor)))
  sampledatachar$Veg2 <- factor(pull(sampledatachar, !!x), levels = c("0", "1"))
  sampledatachar <- mutate(sampledatachar, exposure_new = pull(sampledatachar,!!x))
  sampledatachar <- mutate(sampledatachar, dependent_new = pull(sampledatachar,!!j))
  
  #Statistical Testing
  stat.test <- add_significance(adjust_pvalue(t_test(group_by(sampledatachar, !!x), dependent_new ~ Genotype),method = "bonferroni"),"p.adj")
  stat.test2 <- add_significance(adjust_pvalue(t_test(sampledatachar, dependent_new ~ exposure_new), method = "bonferroni"), "p.adj")
  stat.test3 <- add_significance(adjust_pvalue(t_test(sampledatachar, dependent_new ~ Genotype), method = "bonferroni"), "p.adj")
  
  
  #Graphing First Plot
  plotoutput <- ggbarplot(sampledatachar, x = "Genotype", y = quo_name(j), fill = "Genotype", palette = c("#FF0000", "#00A08A", "#F2AD00"), add = "mean_se", position = position_dodge(0.8))
    stat.test3 <- add_xy_position(stat.test3, fun = "mean_se", x = "Genotype", dodge = 0.8, step.increase = 0.4)
  plotoutput <- plotoutput + stat_pvalue_manual(stat.test3, label = "p.adj.signif", tip.length = 0.0001, size = 6, y.position = c(signif01, signif02, signif03))
  plotoutput <- ggpar(plotoutput, legend = "right") + coord_cartesian(ylim = c(bound_lower, bound_upper))
                       
  #Set 2nd Plot Bounds
  bound_lower2 <- min(exposure$mean_SHBG)-2*max(exposure$se_SHBG)
  bound_lower2[bound_lower2 < 0] <- 0
  bound_upper2 <- max(exposure$mean_SHBG)+6*max(exposure$se_SHBG)
  bound_diff2 <- bound_upper2 - bound_lower2
  
  #Pairwise Bracket Locations
  signif1 <- max(exposure$mean_SHBG[1:3]) + max(exposure$se_SHBG[1:3]) + 0.09*bound_diff2
  signif2 <- max(exposure$mean_SHBG[1:3]) + max(exposure$se_SHBG[1:3]) + 0.18*bound_diff2
  signif3 <- max(exposure$mean_SHBG[1:3]) + max(exposure$se_SHBG[1:3]) + 0.27*bound_diff2
  signif4 <- max(exposure$mean_SHBG[4:6]) + max(exposure$se_SHBG[4:6]) + 0.09*bound_diff2
  signif5 <- max(exposure$mean_SHBG[4:6]) + max(exposure$se_SHBG[4:6]) + 0.18*bound_diff2
  signif6 <- max(exposure$mean_SHBG[4:6]) + max(exposure$se_SHBG[4:6]) + 0.27*bound_diff2
  
  #Second Plot
  testplot <- ggbarplot(sampledatachar, x = quo_name(x), y = quo_name(j), add = "mean_se", fill = "Genotype", position = position_dodge(width = 0.8), palette = c("#FF0000", "#00A08A", "#F2AD00"))
  stat.test <- add_xy_position(stat.test, fun = "mean_se", x = quo_name(x), dodge = 0.8, step.increase = 0.15)
  stat.test2 <- add_xy_position(stat.test2, fun = "mean_se", x = quo_name(x))
  testplot <- testplot + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = bound_diff2/10000, size = 6, y.position = c(signif1, signif2, signif3, signif4, signif5, signif6)) +
                         stat_pvalue_manual(stat.test2, label = "p.adj.signif", tip.length = bound_diff2/10000, y.position = bound_upper2-0.025*bound_diff2, size = 6)
                           #bound_upper2-(max(exposure$mean_SHBG)+max(exposure$se_SHBG))))
  testplot <- ggpar(testplot, legend = "right") + coord_cartesian(ylim = c(bound_lower2, bound_upper2))
  
  #Arrange two plots together
  plot_final <- ggarrange(plotoutput, testplot, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
  plot_final <- annotate_figure(plot_final, top = text_grob(quo_name(i)))
  print(plot_final)
  
  #Export final plot as .png
  ggexport(plot_final, filename = "C:/Users/pierc/OneDrive/Desktop/SHBG-Plot-Output.png")
  
}
