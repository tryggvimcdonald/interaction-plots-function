library(stats)
library(rstatix)
library(utils)
library(tidyverse)
library(plotrix)
library(ggsignif)
library(ggpubr)


interaction_plots <- function(tablename, 
                              pheno, 
                              SNP, 
                              exposure, 
                              nonref_allele, 
                              ref_allele, 
                              num_plots = 3, 
                              export = TRUE)
{
  #tablename - name of the data to be used for the plots
  #pheno - phenotype being measured
  #SNP - SNP name (e.g. rs161896_A)
  #exposure - exposure (e.g. Vegetarianism)
  #nonref_allele - non-reference allele. 0 is interpreted as homozygous non-reference
  #ref_allele - reference allele. 2 is interpreted as homozygous reference
  #num_plots - specify number of plot panels to display, can be 1, 2, or 3. Default is 3.
  #export - FALSE prevents plot output to file
  
  
  #Put variables in quotes if they need them to play nicely with other functions
  pheno <- enquo(pheno)
  SNP <- enquo(SNP)
  exposure <- enquo(exposure)
  nonref_allele <- enquo(nonref_allele)
  ref_allele <- enquo(ref_allele)
  
  #Import data from user
  sampledatachar <- tablename
  
  #Creates a dataframe with mean and sd of SHBG, grouped by allele genotype
  testdataframe <- sampledatachar %>%
                     group_by(!!SNP) %>%
                       summarise(Means = mean(!!pheno), SE = std.error(!!pheno))
  
  #Split data by exposure and genotype then apply mean and sd functions
  exposure_df <- sampledatachar %>%
                    group_by(!!exposure, !!SNP) %>%
                      summarise(Means = mean(!!pheno), SE = std.error(!!pheno), .groups = "keep")
 
  #Print the mean and se values grouped by SNP and exposure
  cat("\n\nMean and standard error for SNP values, grouped by exposure:\n\n")
  cat(format(exposure_df)[c(-1L,-2L,-4L)], sep = "\n")
  
  #First Graph Bound Calculations
  bound_lower <- min(testdataframe$Means)-2*max(testdataframe$SE)
  bound_lower[bound_lower < 0] <- 0
  bound_upper <- max(testdataframe$Means)+3*max(testdataframe$SE)
  bound_diff <- bound_upper - bound_lower
  
  signif01 <- max(testdataframe$Means) + 
              max(testdataframe$SE) + 
              0.09*bound_diff
  
  signif02 <- max(testdataframe$Means) + 
              max(testdataframe$SE) + 
              0.18*bound_diff
  
  signif03 <- max(testdataframe$Means) + 
              max(testdataframe$SE) + 
              0.27*bound_diff
  
  #SNP Allele Naming
  homo_nonref <- paste(quo_name(nonref_allele),
                      quo_name(nonref_allele),
                      sep = "")
  
  het <- paste(quo_name(ref_allele),
               quo_name(nonref_allele),
               sep = "")
  
  homo_ref <- paste(quo_name(ref_allele),
                      quo_name(ref_allele),
                      sep = "")
  
  #2nd Plot Re-leveling
  sampledatachar <- mutate(sampledatachar, Genotype = factor(pull(sampledatachar, !!SNP),
                                                         levels = c(0,1,2),
                                                         labels = c(homo_nonref, het, homo_ref)))
  
  sampledatachar$Veg2 <- factor(pull(sampledatachar, !!exposure),
                                levels = c("0", "1"))
  
  sampledatachar <- mutate(sampledatachar,
                           exposure_new = pull(sampledatachar,!!exposure))
  
  sampledatachar <- mutate(sampledatachar,
                           dependent_new = pull(sampledatachar,!!pheno))
  
  #Statistical Testing
  stat.test <- add_significance(adjust_pvalue(
    t_test(group_by(sampledatachar, !!exposure),
           dependent_new ~ Genotype),
    method = "bonferroni"),
    "p.adj")

  stat.test2 <- add_significance(adjust_pvalue(
    t_test(sampledatachar, 
           dependent_new ~ exposure_new), 
    method = "bonferroni"),
    "p.adj")
  
  stat.test3 <- add_significance(adjust_pvalue(
    t_test(sampledatachar, 
           dependent_new ~ Genotype), 
    method = "bonferroni"), 
    "p.adj")
  
  #Statistics Outputs in Console
  session_options <- options(pillar.neg = FALSE)
  
  stat_output <- select(stat.test3, group1, group2, statistic, p:p.adj.signif)
  cat("\n\nUnstratified SNP statistics and p-values:\n\n")
  cat(format(stat_output)[c(-1L,-3L)], sep = "\n")
  
  stat_output2 <- select(stat.test, 1, group1, group2, statistic, p:p.adj.signif)
  cat("\n\nStatistics and p-values stratified by ",quo_name(exposure),":", "\n\n", sep = "")
  cat(format(stat_output2)[c(-1L,-3L)], sep = "\n")
  
  #Graphing First Plot
  plotoutput <- ggbarplot(sampledatachar, 
                          x = "Genotype", 
                          y = quo_name(pheno), 
                          fill = "Genotype", 
                          palette = c("#FF0000", "#00A08A", "#F2AD00"), 
                          add = "mean_se", 
                          position = position_dodge(0.8))
  stat.test3 <- add_xy_position(stat.test3, 
                                fun = "mean_se", 
                                x = "Genotype", 
                                dodge = 0.8, 
                                step.increase = 0.4)
  plotoutput <- plotoutput + stat_pvalue_manual(stat.test3, 
                                                label = "p.adj.signif", 
                                                tip.length = 0.0001, 
                                                size = 4, 
                                                y.position = c(signif01, signif02, signif03))
  plotoutput <- ggpar(plotoutput, legend = "right") + 
    coord_cartesian(ylim = c(bound_lower, bound_upper))
  
  #Set 2nd Plot Bounds
  bound_lower2 <- min(exposure_df$Means)-2*max(exposure_df$SE)
  bound_lower2[bound_lower2 < 0] <- 0
  bound_upper2 <- max(exposure_df$Means)+6*max(exposure_df$SE)
  bound_diff2 <- bound_upper2 - bound_lower2
  
  #Significance Tip Lengths
  bar_diff_prop <- (max(exposure_df$Means)-min(exposure_df$Means))/bound_diff2  
  
  #Pairwise Bracket Locations
  signif1 <- max(exposure_df$Means[1:3]) + 
             max(exposure_df$SE[1:3]) + 
             0.09*bound_diff2
  
  signif2 <- max(exposure_df$Means[1:3]) + 
             max(exposure_df$SE[1:3]) + 
             0.18*bound_diff2
  
  signif3 <- max(exposure_df$Means[1:3]) + 
             max(exposure_df$SE[1:3]) + 
             0.27*bound_diff2
  
  signif4 <- max(exposure_df$Means[4:6]) + 
             max(exposure_df$SE[4:6]) + 
             0.09*bound_diff2
  
  signif5 <- max(exposure_df$Means[4:6]) + 
             max(exposure_df$SE[4:6]) + 
             0.18*bound_diff2
  
  signif6 <- max(exposure_df$Means[4:6]) + 
             max(exposure_df$SE[4:6]) + 
             0.27*bound_diff2
  
  #Second Plot
  testplot <- ggbarplot(sampledatachar,
                        x = quo_name(exposure),
                        y = quo_name(pheno),
                        add = "mean_se", 
                        fill = "Genotype", 
                        position = position_dodge(width = 0.8),
                        palette = c("#FF0000", "#00A08A", "#F2AD00"))
  
  stat.test <- add_xy_position(stat.test,
                               fun = "mean_se",
                               x = quo_name(exposure),
                               dodge = 0.8,
                               step.increase = 0.15)
  
  stat.test2 <- add_xy_position(stat.test2,
                                fun = "mean_se",
                                x = quo_name(exposure))
  
  testplot <- testplot + stat_pvalue_manual(stat.test,
                                            label = "p.adj.signif",
                                            tip.length = bar_diff_prop*(bar_diff_prop/100)-(bar_diff_prop/1000),
                                            size = 4,
                                            y.position = c(signif1, signif2, signif3, signif4, signif5, signif6)) +
    stat_pvalue_manual(stat.test2,
                       label = "p.adj.signif",
                       tip.length = bar_diff_prop*(bar_diff_prop/100)-(bar_diff_prop/1000),
                       y.position = bound_upper2-0.025*bound_diff2,
                       size = 4)
  
  testplot <- ggpar(testplot, legend = "right") + 
    coord_cartesian(ylim = c(bound_lower2, bound_upper2))
  
  #Arrange plots
  if(num_plots == 1)
  {
    plot_final <- plotoutput
  }
  
  if(num_plots == 2)
  {
    plot_final <- testplot
  }
  
  if(num_plots != 1 && num_plots != 2)
  {
    plot_final <- ggarrange(plotoutput,
                            testplot,
                            nrow = 1,
                            ncol = 2,
                            common.legend = TRUE,
                            legend = "right")
  }
  
  plot_final <- annotate_figure(plot_final,
                                top = text_grob(c(paste(quo_name(pheno),
                                                      "by genotype, stratified by",
                                                      quo_name(exposure)))))
  print(plot_final)
  
  #Export final plot as .png
  if(export != FALSE)
  {
    ggexport(plot_final, filename = paste(getwd(),"/Interaction-Plot-Output.png", sep = ""))
  }
  
}
