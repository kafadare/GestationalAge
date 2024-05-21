GA_boxplot <- function(df, value){
  ggp_out <- ggplot(df, aes(x = preterm, y = {{value}})) +
    geom_boxplot(aes(fill = preterm))}

group_smooth <- function(df, group, y_value, x_value, se = FALSE, method = "loess"){
  ggplot(df, aes(x = {{x_value}}, y = {{y_value}}, color = {{group}})) + geom_point(alpha = 0.2) + geom_smooth(method = method, se = se) +  labs(title = paste0("By ", deparse(substitute(group)))) + theme_minimal()
}
group_smooth_onlyLine <- function(df, group, y_value, x_value, se = FALSE, method = "loess"){
  ggplot(df, aes(x = {{x_value}}, y = {{y_value}}, color = {{group}})) + geom_smooth(method = method, se = se) +  labs(title = paste0("By ", deparse(substitute(group)))) + theme_minimal()
}


group_median <- function(df, group, value, CI = TRUE, coef = 0, outlier = FALSE, notch = TRUE) {
  if(outlier == FALSE){
    outlier_set = NA
  }else if (outlier == TRUE){
    outlier_set = 1
  }
  if (CI == TRUE){
  # Calculate median and confidence interval
  medians <- df %>%
      group_by({{group}}) %>%
      summarise(median = median({{value}}),
                lower_ci = median({{value}}) - 1.58 * IQR({{value}}) / sqrt(n()),
                upper_ci = median({{value}}) + 1.58 * IQR({{value}}) / sqrt(n()))
  # Plot with y-axis limits set to the confidence interval for the median
  ggplot(medians, aes(x = as.factor({{group}}), ymin = lower_ci, lower = median, middle = median, upper = median, ymax = upper_ci)) +
    geom_boxplot(stat = "identity", width = 0.5, notch = notch, outlier.shape = outlier_set) +
    scale_y_continuous(limits = c(min(medians$lower_ci*1.1), max(medians$upper_ci*1.1))) +
    theme(axis.title.x = element_blank(), legend.position = "none") 
  } else if(CI == FALSE) {
    #plot regular boxplot.
  ggplot(df, aes(x = {{group}}, y = {{value}})) +
    geom_boxplot(coef = coef, notch = notch, outlier.shape = outlier_set) + 
    theme(axis.title.x = element_blank(), legend.position = "none") +  
    scale_y_continuous(limits = quantile(df[,{{value}}], c(0.25*1.1, 0.75*1.1)))
  }
}

minQC_by_region_dens <- function(df, title = "QC Distribution by Region", xlab = "Sythseg Auto QC"){
  p_dens <- ggplot(df) +
    geom_density(aes(x = general.white.matter_qc, fill = "WM",color = "WM"), alpha = 0.4) +
    geom_density(aes(x = general.grey.matter_qc, fill = "GM", color = "GM"), alpha = 0.4) +
    geom_density(aes(x = general.csf_qc, fill = "CSF", color = "CSF"), alpha = 0.4) +
    geom_density(aes(x = cerebellum_qc, fill = "Cerebellum", color = "Cerebellum"), alpha = 0.4) +
    geom_density(aes(x = brainstem_qc, fill = "Brainstem", color = "Brainstem"), alpha = 0.4) +
    geom_density(aes(x = thalamus_qc, fill = "Thalamus", color = "Thalamus"), alpha = 0.4) +
    geom_density(aes(x = putamen.pallidum_qc, fill = "Put/Pal", color = "Put/Pal"), alpha = 0.4) +
    geom_density(aes(x = hippocampus.amygdala_qc, fill = "Hip/Amy",color = "Hip/Amy"), alpha = 0.4) +
    geom_density(aes(x = minQC, fill = "minQC"), alpha = 0) +
    scale_fill_manual(values = c("WM" = "red", "GM" = "green", "CSF" = "blue", "Cerebellum" = "yellow", "Brainstem" = "orange", "Thalamus" = "purple", "Put/Pal" = "pink", "Hip/Amy"= "cyan", "minQC" = "black")) +
    scale_color_manual(values = c("WM" = "red", "GM" = "green", "CSF" = "blue", "Cerebellum" = "yellow", "Brainstem" = "orange", "Thalamus" = "purple", "Put/Pal" = "pink", "Hip/Amy"= "cyan", "minQC" = "black"), guide = guide_none()) +
    labs(title = title, y = "Frequency", x = xlab) }

minQC_by_region_hist <- function(df, title = "QC Distribution by Region", xlab = "Sythseg Auto QC"){
  p_hist <- ggplot(df) +
    geom_histogram(aes(x = general.white.matter_qc, fill = "WM"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = general.grey.matter_qc, fill = "GM"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = general.csf_qc, fill = "CSF"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = cerebellum_qc, fill = "Cerebellum"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = brainstem_qc, fill = "Brainstem"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = thalamus_qc, fill = "Thalamus"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = putamen.pallidum_qc, fill = "Put/Pal"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = hippocampus.amygdala_qc, fill = "Hip/Amy"), alpha = 0.4, position = position_dodge()) +
    geom_histogram(aes(x = minQC, fill = "minQC"), alpha = 0) +
    scale_fill_manual(values = c("WM" = "red", "GM" = "green", "CSF" = "blue", "Cerebellum" = "yellow", "Brainstem" = "orange", "Thalamus" = "purple", "Put/Pal" = "pink", "Hip/Amy"= "cyan", "minQC" = "black")) +
    labs(title = title, y = "Frequency", x = xlab) }


#if output of centile curves is better organized, could potentially edit this function to be customizable to whichever "grouping" (sex, GA, interaction) I want to customize
growthChart_plot <- function(p, df, centileCurves, xlab = "Age at scan (years)", title = NULL, tickMarks = NULL, tickLabels = NULL, print = FALSE, by.sex = FALSE, by.GAbins= FALSE, by.preterm = FALSE){
  #define title if NULL
  if (is.null(title)){
    title <- paste0("Sample Growth Chart for ", p)
  }
  if (is.null(tickMarks)){
    tickMarks <- c()
    for (year in c(0, 1, 2, 5, 10, 20)){ # years
      tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10)) #currently has the standard 280 adjustment, but is this reasonable? Perhaps for the scale yes.
    }
    if (is.null(tickLabels)){
    tickLabels <- c("Birth", "1", "2", "5", "10", "20")
    }
  }
# Plot the original data and the set of centile curves on a figure
  if(by.sex == FALSE & by.preterm == FALSE & by.GAbins == FALSE){
    plot <- ggplot() +
      geom_point(aes(x=df$logAge, df[, p]), alpha=0.1) +
      geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$median), alpha=0.4) +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$median), alpha=0.6) +
      geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$median), alpha=0.8) +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$median)) +
      geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$median), alpha=0.8) +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$median), alpha=0.6) +
      geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$median), alpha=0.4) +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(tickMarks[[1]], max(centileCurves[[1]]$logAge))) +
      labs(title=title) +
      xlab(xlab) +
      ylab(paste0(p)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 18))
      if(print == TRUE){
        print(plot)
      }
    return(plot)
  } else if (by.sex == TRUE) {
    plot_sex <- ggplot() +
      geom_point(aes(x=df$logAge, df[, p], color = df$sex), alpha=0.5) +
      scale_color_manual(values = c("female" = "red", "male" = "blue")) +
      ##male curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$male), alpha=0.4, color = "blue") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$male), alpha=0.6, color = "blue") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$male), alpha=0.8, color = "blue") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$male), color = "blue") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$male), alpha=0.8, color = "blue") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$male), alpha=0.6, color = "blue") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$male), alpha=0.4, color = "blue") +
      ##female curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$female), alpha=0.4, color = "red") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$female), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$female), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$female), color = "red") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$female), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$female), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$female), alpha=0.4, color = "red") +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(tickMarks[[1]], max(centileCurves[[1]]$logAge))) +
      labs(title=paste0(title, "By Sex")) +
      xlab(xlab) +
      ylab(paste0(p)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 18))
      if(print == TRUE){
        print(plot_sex)
      }
    return(plot_sex)
  }
  else if (by.GAbins == TRUE) {
    plot_GAbins <- ggplot() +
      geom_point(aes(x=df$logAge, df[, p], color = df$GAbins_recode), alpha=0.25) + 
      scale_color_manual(values = c("VPM" = "red", "LPM" = "blue", "Term" = "green")) +
      ##VPM curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$VPM), alpha=0.4, color = "red") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$VPM), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$VPM), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$VPM), color = "red") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$VPM), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$VPM), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$VPM), alpha=0.4, color = "red") +
      ##LPM curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$LPM), alpha=0.4, color = "blue") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$LPM), alpha=0.6, color = "blue") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$LPM), alpha=0.8, color = "blue") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$LPM), color = "blue") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$LPM), alpha=0.8, color = "blue") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$LPM), alpha=0.6, color = "blue") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$LPM), alpha=0.4, color = "blue") +
      ##Term curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$Term), alpha=0.4, color = "green") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$Term), alpha=0.6, color = "green") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$Term), alpha=0.8, color = "green") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$Term), color = "green") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$Term), alpha=0.8, color = "green") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$Term), alpha=0.6, color = "green") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$Term), alpha=0.4, color = "green") +
      
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(tickMarks[[1]], max(centileCurves[[1]]$logAge))) +
      labs(title=paste0(title, "By Preterm Status")) +
      xlab(xlab) +
      ylab(paste0(p)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 18))
    if(print == TRUE){
      print(plot_GAbins)
    }
    return(plot_GAbins)
  }
  else if (by.preterm == TRUE) {
    plot_preterm <- ggplot() +
      geom_point(aes(x=df$logAge, df[, p], color = df$preterm), alpha=0.25) + 
      scale_color_manual(values = c("Preterm" = "red", "Term" = "green")) +
      ##Pterm curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$Pterm), alpha=0.4, color = "red") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$Pterm), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$Pterm), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$Pterm), color = "red") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$Pterm), alpha=0.8, color = "red") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$Pterm), alpha=0.6, color = "red") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$Pterm), alpha=0.4, color = "red") +
      ##Term curves
      #geom_line(aes(x=centileCurves[[1]]$logAge, y=centileCurves[[1]]$Term), alpha=0.4, color = "green") +
      geom_line(aes(x=centileCurves[[2]]$logAge, y=centileCurves[[2]]$Term), alpha=0.6, color = "green") +
      #geom_line(aes(x=centileCurves[[3]]$logAge, y=centileCurves[[3]]$Term), alpha=0.8, color = "green") +
      geom_line(aes(x=centileCurves[[4]]$logAge, y=centileCurves[[4]]$Term), color = "green") +
      #geom_line(aes(x=centileCurves[[5]]$logAge, y=centileCurves[[5]]$Term), alpha=0.8, color = "green") +
      geom_line(aes(x=centileCurves[[6]]$logAge, y=centileCurves[[6]]$Term), alpha=0.6, color = "green") +
      #geom_line(aes(x=centileCurves[[7]]$logAge, y=centileCurves[[7]]$Term), alpha=0.4, color = "green") +
      
      scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                         limits=c(tickMarks[[1]], max(centileCurves[[1]]$logAge))) +
      labs(title=paste0(title, "By Preterm Status")) +
      xlab(xlab) +
      ylab(paste0(p)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 18))
    if(print == TRUE){
      print(plot_preterm)
    }
    return(plot_preterm)
  }
}