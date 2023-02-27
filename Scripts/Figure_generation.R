rm(list = ls())
color_code <- c("#6fb5a0","#ec9778","#89a0cc","#c8b1d6","#Bb330c")
packages_needed <- c("survival","dplyr","ggplot2","survminer","cowplot","gridExtra","pROC","swimplot","tidyverse","reshape2","survcomp","epiR","eulerr","data.table","caret")
lapply(packages_needed, require, character.only = TRUE)
##### P1 (7 days postsurgical) P2 (6 months postsurgical)

load("CRC_input_data.Rdata")

theme_plot <- theme(axis.line = element_line(colour = "black", size = 1.25),
                    plot.title=element_text(hjust=0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    legend.title = element_text(size = 13, color = "black", face = "bold"),
                    legend.text = element_text(size = 13, color = "black", face = "bold"),
                    axis.text.x = element_text(size = 13, color = "black", face = "bold"),
                    axis.text.y = element_text(size = 13, color = "black", face = "bold"),
                    axis.title.x = element_text(size = 13, color = "black", face = "bold"),
                    axis.title.y = element_text(size = 13, color = "black", face = "bold"))

print_stats <- function(df){
  t <- format(df * 100, digits = 3)
  results <- paste0(t[[1]], "% (", t[[2]], "-", t[[3]], "%)")
  return(results)
}

#####data processing and table generation

#### P1
P1_df <- read.csv(file = "P1_output.csv", stringsAsFactors = FALSE)
P1_df$PD_event <- as.numeric(as.logical(P1_df$PD_event))
P1_model <- inner_join(P1_df, clinical_info[,c("SampleID","Patient ID")],c("SampleID"="SampleID"))
P1_rocobj <- roc(P1_model$PD_event, P1_model$Cox_pred, direction = "<")
P1_threshold <- coords(P1_rocobj, "local maximas", transpose = FALSE)
P1_threshold <- P1_threshold[P1_threshold$specificity >= 0.9,]
P1_threshold <- P1_threshold[order(P1_threshold$sensitivity, decreasing = TRUE),]
P1_Index <- P1_threshold[1,]$threshold
P1_high_risk <- P1_model[P1_model$Cox_pred > P1_Index, ]
P1_high_risk_id <- P1_high_risk$SampleID
P1_model$risk <- ifelse(P1_model$SampleID %in% P1_high_risk_id, "High", "Low")
P1_model$Time <- P1_model$DFS - P1_model$DFS_adjusted
P1_model <- inner_join(P1_model, MRD_sample_point,c("Patient ID" = "Patient_ID", "Time"="Time"))
P1_model$joined <- 1
P1_model[which(P1_model$risk == "Low" & P1_model$Event == "Negative"),]$joined = 0
P1_model$ctDNA = 1
P1_model[P1_model$Event == "Negative",]$ctDNA = 0

#####for table stats
table(P1_model[,c("PD_event","risk")])
table(P1_model[,c("PD_event","joined")])


##### frag model
tmp_test <- P1_model[,c("risk","PD_event")]
tmp_test[tmp_test$risk == "Low",]$risk = 0
tmp_test[tmp_test$risk == "High",]$risk = 1
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]

##### ctDNA

tmp_test <- P1_model[,c("ctDNA","PD_event")]
colnames(tmp_test) <- c("risk","PD_event")
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]


##### combined
tmp_test <- P1_model[,c("joined","PD_event")]
colnames(tmp_test) <- c("risk","PD_event")
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]



#### P2
P2_df <- read.csv(file = "P2_output.csv", stringsAsFactors = FALSE)
P2_df$PD_event <- as.numeric(as.logical(P2_df$PD_event))
P2_model <- inner_join(P2_df, clinical_info[,c("SampleID","Patient ID")],c("SampleID"="SampleID"))
P2_rocobj <- roc(P2_model$PD_event, P2_model$Cox_pred, direction = "<")
P2_threshold <- coords(P2_rocobj, "local maximas", transpose = FALSE)
P2_threshold <- P2_threshold[P2_threshold$specificity >= 0.9,]
P2_threshold <- P2_threshold[order(P2_threshold$sensitivity, decreasing = TRUE),]
P2_Index <- P2_threshold[1,]$threshold
P2_high_risk <- P2_model[P2_model$Cox_pred > P2_Index, ]
P2_high_risk_id <- P2_high_risk$SampleID
P2_model$risk <- ifelse(P2_model$SampleID %in% P2_high_risk_id, "High", "Low")
P2_model$Time <- P2_model$DFS - P2_model$DFS_adjusted
P2_model <- inner_join(P2_model, MRD_sample_point,c("Patient ID" = "Patient_ID", "Time"="Time"))
P2_model$joined <- 1
P2_model[which(P2_model$risk == "Low" & P2_model$Event == "Negative"),]$joined = 0
P2_model$ctDNA = 1
P2_model[P2_model$Event == "Negative",]$ctDNA = 0

#####for table stats
table(P2_model[,c("PD_event","risk")])
table(P2_model[,c("PD_event","joined")])


##### frag model
tmp_test <- P2_model[,c("risk","PD_event")]
tmp_test[tmp_test$risk == "Low",]$risk = 0
tmp_test[tmp_test$risk == "High",]$risk = 1
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]


##### ctDNA

tmp_test <- P2_model[,c("ctDNA","PD_event")]
colnames(tmp_test) <- c("risk","PD_event")
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]


##### combined
tmp_test <- P2_model[,c("joined","PD_event")]
colnames(tmp_test) <- c("risk","PD_event")
tmp_test$risk = factor(tmp_test$risk, levels = c(1,0))
tmp_test$PD_event = factor(tmp_test$PD_event, levels = c(1,0))

table(tmp_test)
tmp_epi <- epi.tests(dat = table(tmp_test[,c("risk","PD_event")]), conf.level = 0.95)
print_stats(tmp_epi$detail$se)
print_stats(tmp_epi$detail$sp)
print_stats(tmp_epi$detail$pv.pos)
print_stats(tmp_epi$detail$pv.neg)
print_stats(tmp_epi$detail$diag.ac)

tmp_cm <- confusionMatrix(tmp_test$risk, tmp_test$PD_event, positive = "1", mode = "everything")
tmp_cm[["byClass"]]["Precision"]
tmp_cm[["byClass"]]["Recall"]



#####Figure 2 (HRs manually added)

## Fragment

P1_fit <- survfit(formula = Surv(DFS, PD_event)~ risk, data = P1_model)
P1_KM <- ggsurvplot(P1_fit, data = P1_model, risk.table = T, pval = T,
                    break.time.by = 365, surv.median.line = 'hv', 
                    font.title = c(14, "bold", "black"),
                    conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2, title = "7 days postsurgical",
                    # font.main = c(16, "bold", "black"),
                    # font.x = c(14, "bold", "black"),
                    # font.y = c(14, "bold", "black"),
                    # font.tickslab = c(14, "bold", "black"),
                    censor.size=4, size = 1.25,
                    ggtheme = theme_plot,
                    tables.theme = theme_plot,
                    axis.line = element_line(colour = 'black', size = 2),
                    palette = c("red", "blue"),
                    legend.labs = c("High risk","Low risk"),
                    legend.title = "",
                    risk.table.title = "",
                    risk.table.y.text.col = TRUE,
                    risk.table.y.text = FALSE,
                    tables.col = "strata"
) 

P2_fit <- survfit(formula = Surv(DFS, PD_event)~ risk, data = P2_model)
P2_KM <- ggsurvplot(P2_fit, data = P2_model, risk.table = T, pval = T,
                    break.time.by = 365, surv.median.line = 'hv', 
                    font.title = c(14, "bold", "black"),
                    conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2, title = "6 months postsurgical",
                    # font.main = c(16, "bold", "black"),
                    # font.x = c(14, "bold", "black"),
                    # font.y = c(14, "bold", "black"),
                    # font.tickslab = c(14, "bold", "black"),
                    censor.size=4, size = 1.25,
                    ggtheme = theme_plot,
                    tables.theme = theme_plot,
                    axis.line = element_line(colour = 'black', size = 2),
                    palette = c("red", "blue"),
                    legend.labs = c("High risk","Low risk"),
                    legend.title = "",
                    risk.table.title = "",
                    risk.table.y.text.col = TRUE,
                    risk.table.y.text = FALSE,
                    tables.col = "strata"
) 

####
Figure_3_top_plot  <- ggarrange(plotlist=list(P1_KM$plot, P2_KM$plot),
                                ncol = 2, nrow = 1, labels = c("A","B"))

Figure_3_top_risk <- ggarrange(plotlist=list(P1_KM$table, P2_KM$table),
                               ncol = 2, nrow = 1)

Figure_3_top_title <- ggdraw() + draw_label(paste0("Fragmentation"), fontface='bold')

Figure_3_top_panel <- plot_grid(Figure_3_top_title, Figure_3_top_plot, Figure_3_top_risk,
                                ncol=1,
                                rel_heights=c(0.1, 1, 0.5))

####ctDNA
P1_model$ctDNA = factor(P1_model$ctDNA, levels = c(1,0))
P1_mut_fit <- survfit(formula = Surv(DFS, PD_event)~ ctDNA, data = P1_model)
P1_mut_KM <- ggsurvplot(P1_mut_fit, data = P1_model, risk.table = T, pval = T,
                        break.time.by = 365, surv.median.line = 'hv', 
                        font.title = c(14, "bold", "black"),
                        conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2,
                        title = "7 days postsurgical",
                        # font.main = c(16, "bold", "black"),
                        # font.x = c(14, "bold", "black"),
                        # font.y = c(14, "bold", "black"),
                        # font.tickslab = c(14, "bold", "black"),
                        censor.size=4, size = 1.25,
                        ggtheme = theme_plot,
                        tables.theme = theme_plot,
                        axis.line = element_line(colour = 'black', size = 2),
                        palette = c("red", "blue"),
                        legend.labs = c("Pos","Neg"),
                        legend.title = "",
                        risk.table.title = "",
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE,
                        tables.col = "strata"
) 

P2_model$ctDNA = factor(P2_model$ctDNA, levels = c(1,0))
P2_mut_fit <- survfit(formula = Surv(DFS, PD_event)~ ctDNA, data = P2_model)
P2_mut_KM <- ggsurvplot(P2_mut_fit, data = P2_model, risk.table = T, pval = T,
                        break.time.by = 365, surv.median.line = 'hv', 
                        font.title = c(14, "bold", "black"),
                        conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2, 
                        title = "6 months postsurgical",
                        # font.main = c(16, "bold", "black"),
                        # font.x = c(14, "bold", "black"),
                        # font.y = c(14, "bold", "black"),
                        # font.tickslab = c(14, "bold", "black"),
                        censor.size=4, size = 1.25,
                        ggtheme = theme_plot,
                        tables.theme = theme_plot,
                        axis.line = element_line(colour = 'black', size = 2),
                        palette = c("red", "blue"),
                        legend.labs = c("Pos","Neg"),
                        legend.title = "",
                        risk.table.title = "",
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE,
                        tables.col = "strata"
) 

Figure_3_mid_plot  <- ggarrange(plotlist=list(P1_mut_KM$plot,P2_mut_KM$plot),
                                ncol = 2, nrow = 1, labels = c("C","D"))

Figure_3_mid_risk <- ggarrange(plotlist=list(P1_mut_KM$table,P2_mut_KM$table),
                               ncol = 2, nrow = 1)

Figure_3_mid_title <- ggdraw() + draw_label(paste0("Mutation"), fontface='bold')

Figure_3_mid_panel <- plot_grid(Figure_3_mid_title, Figure_3_mid_plot, Figure_3_mid_risk,
                                ncol=1,
                                rel_heights=c(0.1, 1, 0.5))

####joined
P1_model$joined = factor(P1_model$joined, levels = c(1,0))
P1_ctDNA_fit <- survfit(formula = Surv(DFS, PD_event)~ joined, data = P1_model)
P1_ctDNA_KM <- ggsurvplot(P1_ctDNA_fit, data = P1_model, risk.table = T, pval = T,
                          break.time.by = 365, surv.median.line = 'hv', 
                          font.title = c(14, "bold", "black"),
                          conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2,
                          title = "7 days postsurgical",
                          # font.main = c(16, "bold", "black"),
                          # font.x = c(14, "bold", "black"),
                          # font.y = c(14, "bold", "black"),
                          # font.tickslab = c(14, "bold", "black"),
                          censor.size=4, size = 1.25,
                          ggtheme = theme_plot,
                          tables.theme = theme_plot,
                          axis.line = element_line(colour = 'black', size = 2),
                          palette = c("red", "blue"),
                          legend.labs = c("High risk/Pos","Low risk&Neg"),
                          legend.title = "",
                          risk.table.title = "",
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE,
                          tables.col = "strata"
) 

##Panel B, P def 
P2_model$joined = factor(P2_model$joined, levels = c(1,0))
P2_ctDNA_fit <- survfit(formula = Surv(DFS, PD_event)~ joined, data = P2_model)
P2_ctDNA_KM <- ggsurvplot(P2_ctDNA_fit, data = P2_model, risk.table = T, pval = T,
                          break.time.by = 365, surv.median.line = 'hv', 
                          font.title = c(14, "bold", "black"),
                          conf.int = F, ylab = 'PFS probability', xlab = 'Time (Days)', cex.axis = 2, 
                          title = "6 months postsurgical",
                          # font.main = c(16, "bold", "black"),
                          # font.x = c(14, "bold", "black"),
                          # font.y = c(14, "bold", "black"),
                          # font.tickslab = c(14, "bold", "black"),
                          censor.size=4, size = 1.25,
                          ggtheme = theme_plot,
                          tables.theme = theme_plot,
                          axis.line = element_line(colour = 'black', size = 2),
                          palette = c("red", "blue"),
                          legend.labs = c("High risk/Pos","Low risk&Neg"),
                          legend.title = "",
                          risk.table.title = "",
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE,
                          tables.col = "strata"
) 

Figure_3_bottom_plot  <- ggarrange(plotlist=list(P1_ctDNA_KM$plot,P2_ctDNA_KM$plot),
                                   ncol = 2, nrow = 1, labels = c("E","F"))

Figure_3_bottom_risk <- ggarrange(plotlist=list(P1_ctDNA_KM$table,P2_ctDNA_KM$table),
                                  ncol = 2, nrow = 1)

Figure_3_bottom_title <- ggdraw() + draw_label(paste0("Fragmentation + Mutation"), fontface='bold')

Figure_3_bottom_panel <- plot_grid(Figure_3_bottom_title, Figure_3_bottom_plot, Figure_3_bottom_risk,
                                   ncol=1,
                                   rel_heights=c(0.1, 1, 0.5))


ggsave(filename = "Figure_3.pdf",
       plot_grid(Figure_3_top_panel, Figure_3_mid_panel, Figure_3_bottom_panel, nrow = 3),
       dpi = 300,width = 9,height = 18)

##### Figure 4, swimmer plot

swim_df <- rbind(P1_model,P2_model)
colnames(swim_df) <- gsub("Patient ID","Patient_ID",colnames(swim_df))
swim_df$Event = swim_df$risk
swim_df$Patient_ID <- as.character(swim_df$Patient_ID)

######## swimmer plot
fig4_swimmer_plot <- function(x, y,z) {
  swimmer_plot <- ggplot(x) +
    
    geom_line(
      data = y,
      aes(x = follow_up_day, y = Patient_ID, group = Patient_ID),
      size = 0.5,
      color = "grey55"
    ) +
    geom_segment(
      data = y,
      aes(
        x = AT_start_day,
        xend = AT_end_day,
        y = Patient_ID,
        yend = Patient_ID,
        color = AT_type
      ),
      alpha = 0.75,
      size = 2.5
    ) +
    scale_color_manual(
      values = c(
        Chemo = "#80B1D3",
        RT = "#FB8072",
        "Chemo+RT" = "#BEBADA",
        "Chemo+Targeted" = "#8DD3C7"
      )
    ) +
    
    geom_point(
      aes(
        x = as.numeric(Time),
        y = Patient_ID,
        shape = Event,
        fill = Event
      ),
      size = 3,
      color = "grey33"
    ) +
    scale_shape_manual(values = c(
      High = 21,
      Relapse = 4,
      Low = 21,
      Decease = 6
    )) +
    scale_fill_manual(
      values = c(
        High = "firebrick3",
        Low = "steelblue3",
        Relapse = "black",
        Decease = "black"
      )
    ) +
    
    geom_segment(
      data = y,
      aes(
        x = AT_start_day_only,
        xend = AT_start_day_only + 5,
        y = Patient_ID,
        yend = Patient_ID
      ),
      color = "green",
      size = 2.5
    ) +
    scale_x_continuous(breaks = c(0, 180, 360, 540, 720, 900, 1080)) +
    xlab("Time (days)") +
    ylab("Patient ID") +
    theme_classic() +
    theme(legend.position = z, legend.box = "vertical") +
    guides(fill = "none") + 
    
    geom_vline(xintercept = 0,
               linetype = 2,
               alpha = 0.5)
  
  swimmer_plot
  return(swimmer_plot)
}


tmp <- inner_join(swim_df[,c("Patient_ID","Time","Event")],MRD_sample_point[,c(1:2,4:6)],c("Patient_ID","Time"))
tmp_extra <- MRD_sample_point[which(MRD_sample_point$Event == "Decease"),]
tmp <- rbind(tmp,tmp_extra)

tmp_extra <- P1_model[P1_model$PD_event == 1,c("Patient ID","DFS")]
colnames(tmp_extra) <- c("Patient_ID","Time")
tmp_extra$Event <- "Relapse"
tmp_extra$AT_start_day <- NA
tmp_extra$AT_end_day <- NA
tmp_extra$AT_type <- NA

tmp <- rbind(tmp,tmp_extra)

tmp <- tmp[tmp$Patient_ID %in% P1_model[order(P1_model$PD_event),]$`Patient ID`,]


differ_1 <- setdiff(MRD_followup[!is.na(MRD_followup$AT_start_day_only),]$Patient_ID,
                    P1_model[P1_model$Therapy == 1,]$`Patient ID`)

differ_2 <- setdiff(P1_model[P1_model$Therapy == 1,]$`Patient ID`,
                    MRD_followup[!is.na(MRD_followup$AT_start_day_only),]$Patient_ID)


tmp[tmp$Patient_ID %in% differ_1,c(4:6)] <- NA


figure_4_panel_A <- fig4_swimmer_plot(tmp[tmp$Patient_ID %in% P1_model[P1_model$PD_event == 1,]$`Patient ID`,], 
                                      MRD_followup[MRD_followup$Patient_ID %in% P1_model[P1_model$PD_event == 1,]$`Patient ID`,],
                                      "top")

figure_4_panel_B <- fig4_swimmer_plot(tmp[tmp$Patient_ID %in% P1_model[P1_model$PD_event == 0,]$`Patient ID`,],
                                      MRD_followup[MRD_followup$Patient_ID %in% P1_model[P1_model$PD_event == 0,]$`Patient ID`,],
                                      "NA")

figure_4_swimmer <- plot_grid(figure_4_panel_A, figure_4_panel_B, ncol = 1, rel_heights = c(1,2), labels = c("A","B"))

ggsave(filename = "Figure_4.pdf",
       figure_4_swimmer,dpi = 300,width = 8,height = 14)



##### Figure_6

##### Figure 6, sens/spec bar plot + 2 Venn diagram

### Panel A, sens/specs bar plot
stats_df <- data.frame(Value = numeric(), 
                       Lower = numeric(), 
                       Upper = numeric(), 
                       type = character(), 
                       label = character(),
                       Group = character())

tmp <- P1_model
tmp$Predict <- 0
tmp[tmp$risk == "High",]$Predict <- 1
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days","Fragmentation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days","Fragmentation")
stats_df <- rbind(stats_df, tmp_df)

tmp <- P2_model
tmp$Predict <- 0
tmp[tmp$risk == "High",]$Predict <- 1
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","6 months","Fragmentation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","6 months","Fragmentation")
stats_df <- rbind(stats_df, tmp_df)

#mutation only
tmp <- P1_model
tmp$Predict <- tmp$ctDNA
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days","Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days","Mutation")
stats_df <- rbind(stats_df, tmp_df)

tmp <- P2_model
tmp$Predict <- tmp$ctDNA
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","6 months","Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","6 months","Mutation")
stats_df <- rbind(stats_df, tmp_df)


#with ctDNA
tmp <- P1_model
tmp$Predict <- tmp$joined
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days","Fragmentation + Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days","Fragmentation + Mutation")
stats_df <- rbind(stats_df, tmp_df)

tmp <- P2_model
tmp$Predict <- tmp$joined
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","6 months","Fragmentation + Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","6 months","Fragmentation + Mutation")
stats_df <- rbind(stats_df, tmp_df)

stats_df$Value <- as.numeric(stats_df$Value)
stats_df$Lower <- as.numeric(stats_df$Lower)
stats_df$Upper <- as.numeric(stats_df$Upper)
stats_df$Group <- factor(stats_df$Group, levels = c("Fragmentation","Mutation","Fragmentation + Mutation"))
stats_df$label <- factor(stats_df$label, levels = c("7 days","6 months"))


sens_plot <- ggplot(data=stats_df[stats_df$label != "Combined timepoints",], aes(x=label, y=Value,fill=type)) +
  geom_bar(colour="black",stat="identity", position = position_dodge(0.7), width = 0.5)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,position = position_dodge(0.7)) + 
  ylab("Sensitivity/Sepcificity") + xlab("Postsurgical") +
  facet_wrap(~ Group,ncol=3) + ylim(0,1) + 
  theme(axis.line = element_line(colour = "black", size = 1.25),
        strip.text.x = element_text(size = 13, colour = "black", face = "bold"),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 13, color = "black", face = "bold"),
        legend.text = element_text(size = 13, color = "black", face = "bold"),
        axis.text.x = element_text(size = 13, color = "black", face = "bold"),
        axis.text.y = element_text(size = 13, color = "black", face = "bold"),
        axis.title.x = element_text(size = 13, color = "black", face = "bold"),
        axis.title.y = element_text(size = 13, color = "black", face = "bold"))

tmp <- P1_model
P1_venn_fit <- euler(c(Fragmenation = length(tmp[tmp$risk == "High",]$SampleID) - length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID),
                       Mutation = length(tmp[tmp$ctDNA == 1,]$SampleID) - length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID),
                       "Fragmenation&Mutation" = length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID)))
P1_venn_plot <- plot(P1_venn_fit, main = "7 days postsurgical",
                     fills = list(fill=c("steelblue4", "darkgoldenrod1"), alpha = 0.5),
                     edges = FALSE, 
                     # labels = list(col = "black", face = "bold", font = 12),
                     legend = TRUE, 
                     adjust_labels =  FALSE,
                     quantities= TRUE)

tmp <- P2_model
P2_venn_fit <- euler(c(Fragmenation = length(tmp[tmp$risk == "High",]$SampleID) - length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID),
                       Mutation = length(tmp[tmp$ctDNA == 1,]$SampleID) - length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID),
                       "Fragmenation&Mutation" = length(tmp[tmp$risk == "High" & tmp$ctDNA == 1,]$SampleID)))
P2_venn_plot <- plot(P2_venn_fit, main = "6 months postsurgical",
                     fills = list(fill=c("steelblue4", "darkgoldenrod1"), alpha = 0.5),
                     edges = FALSE, 
                     # labels = list(col = "black", face = "bold", font = 12),
                     legend = TRUE, 
                     adjust_labels =  FALSE,
                     quantities= TRUE)

Figure_6_bottom_plot <- ggarrange(plotlist=list(P1_venn_plot, P2_venn_plot),
                                  ncol = 2, nrow = 1, labels = c("A","B"))

ggsave(filename = "Figure_6.pdf",
       plot_grid(Figure_6_bottom_plot, sens_plot, ncol=1, rel_heights=c(0.75, 1), labels = c("","C")),
       dpi = 300,width = 12,height = 6)





##### Figure 2
###### Figure 2 panel A ROC, panel B boxplot

tmp_P1 <- P1_model
tmp_P1$Group = "Progressed"
tmp_P1[tmp_P1$PD_event == 0,]$Group = "Not progressed"

figure_2_panel_C <- ggplot(tmp_P1,aes(x=Group,y=Cox_pred, fill = Group)) + geom_boxplot(outlier.shape = NA) +
  xlab("7 days postsurgical") + ylab("Risk score") + theme_bw() +  ylim(c(-20, 80)) + 
  scale_fill_manual(values = color_code[1:2]) +
  geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
  geom_hline(yintercept=P1_Index, linetype="dashed") +
  theme(legend.position = c(0.25,0.85),text=element_text(face='bold'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


tmp_P2 <- P2_model
tmp_P2$Group = "Progressed"
tmp_P2[tmp_P2$PD_event == 0,]$Group = "Not progressed"

figure_2_panel_D <- ggplot(tmp_P2,aes(x=Group,y=Cox_pred, fill = Group)) + geom_boxplot(outlier.shape = NA) +
  xlab("6 months postsurgical") + ylab("Risk score") + theme_bw() + ylim(c(-20, 80)) +
  scale_fill_manual(values = color_code[1:2]) +
  geom_point(aes(fill = Group),size = 1.5, shape = 21, position = position_jitterdodge()) +
  geom_hline(yintercept=P2_Index, linetype="dashed") +
  theme(legend.position = c(0.25,0.85),text=element_text(face='bold'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggrocs <- function(rocs, breaks = seq(0,1,0.1), legendTitel = "Legend", fontsize = 12, line_alaph = 1, line_size = 1, color = color_code[2:5]) {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {
    require(plyr)

    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        names = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringAsFactors = T
      )
    })
    
    aucAvg <- mean(sapply(rocs, "[[", "auc"))
    
    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = names)) +
      geom_line(alpha = line_alaph, size = line_size) +
      #scale_fill_brewer(palette="Dark2") + 
      scale_color_manual(values= color) + 
      geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5, size = 0.75, colour = "black") + 
      geom_step() +
      scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks) + 
      scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks) +
      theme_bw() +
      #theme_grey() + 
      #theme_linedraw() + 
      #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.text=element_text(size=fontsize),
            legend.title = element_text(size = fontsize)) + 
      coord_equal() + 
      theme(legend.position = c(0.65, 0.2)) + 
      #annotate("text", x = 0.1, y = 0.1, vjust = 0, label = paste("AUC =",sprintf("%.3f",aucAvg))) +
      guides(colour = guide_legend(legendTitel), size = fontsize) +
      #theme(axis.ticks = element_line(color = "grey80"))
      theme(axis.title = element_text(size = fontsize), axis.text = element_text(size = fontsize))
    #guides(colour = guide_legend(legendTitel)) 
    
    rocPlot
  }
}

roc.plot.tmp <- list()

t <- format(roc(PD_event~Cox_pred,tmp_P1,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.tmp[[paste0("7 days postsurgical, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~Cox_pred,tmp_P1,levels=c('0','1'),
                                                                                    percent=F,smooth=F,ci=T)

t <- format(roc(PD_event~Cox_pred,tmp_P2,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.tmp[[paste0("6 months postsurgical, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~Cox_pred,tmp_P2,levels=c('0','1'),
                                                                                      percent=F,smooth=F,ci=T)


figure_2_panel_A <- ggrocs(rocs = roc.plot.tmp, breaks = seq(0,1,0.2), legendTitel = "LOOCV\n AUC (95% CI)")

ggsave(filename = "Figure_2.pdf",
       plot_grid(figure_2_panel_A, figure_2_panel_C, figure_2_panel_D, ncol= 3, labels = c("A","B","C"), rel_widths = c(6,4,4)),
       dpi = 300,width = 14,height = 6)

#####Figure F5
P1_correlation_df <- inner_join(P1_model,maxVAF_P1,c("Patient ID"="Patient_ID"))
P1_correlation_df$PD_event <- factor(P1_correlation_df$PD_event, levels = c(1,0), labels = c("Yes","No"))

F5_panel_A <- ggplot(P1_correlation_df, aes(x=maxVAF, y=Cox_pred, 
                                            shape = PD_event,
                                            color = risk
)) + 
  geom_point()+
  geom_smooth(method="lm", color="black", se = TRUE, aes(group = 1))+
  labs(title="7 days postsurgical",
       x="max VAF", y = "Risk score")+
  theme(
    # legend.position = c(0.85,0.85),
    text=element_text(face='bold'),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  stat_cor(method = "pearson", aes(group = 1))

P2_correlation_df <- inner_join(P2_model,maxVAF_P2,c("Patient ID"="Patient_ID"))
P2_correlation_df$PD_event <- factor(P2_correlation_df$PD_event, levels = c(1,0), labels = c("Yes","No"))

F5_panel_B <- ggplot(P2_correlation_df, aes(x=maxVAF, y=Cox_pred, 
                                            shape = PD_event,
                                            color = risk
)) + 
  geom_point()+
  geom_smooth(method="lm", color="black", se = TRUE, aes(group = 1))+
  labs(title="6 months postsurgical",
       x="max VAF", y = "Risk score")+
  theme(
    # legend.position = c(0.85,0.85),
    text=element_text(face='bold'),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  stat_cor(method = "pearson",aes(group = 1))

P1_feature_data <- read.csv(file = "P1_feature.csv", stringsAsFactors = FALSE)

#####P1_ratio
P1_feature_ratio <- P1_feature_data[,-c(2:9)]
rownames(P1_feature_ratio) <- P1_feature_ratio$SampleID
P1_feature_ratio <- P1_feature_ratio[,-1]
P1_shortA <- P1_feature_ratio[,c(grep("shortA",colnames(P1_feature_ratio)))]
P1_longA <- P1_feature_ratio[,c(grep("longA",colnames(P1_feature_ratio)))]
P1_ratioA <- P1_shortA/P1_longA
P1_shortB <- P1_feature_ratio[,c(grep("shortB",colnames(P1_feature_ratio)))]
P1_longB <- P1_feature_ratio[,c(grep("longB",colnames(P1_feature_ratio)))]
P1_ratioB <- P1_shortB/P1_longB
P1_feature_ratio <- cbind(P1_ratioA, P1_ratioB)

P1_feature_ratio_zscore <- P1_feature_ratio

P1_feature_ratio$SampleID <- row.names(P1_feature_ratio)
P1_feature_ratio <- inner_join(P1_model[,c("SampleID","risk","PD_event")],P1_feature_ratio,c("SampleID"="SampleID"))
P1_feature_ratio <- melt(P1_feature_ratio, id.vars = c("SampleID","risk","PD_event"))
P1_feature_ratio$PD_event <- factor(P1_feature_ratio$PD_event,  levels = c(1,0), labels = c("Yes","No"))

fragcnv_bin = read.csv("fragbin_cnvbin_intersect.bed",header=F,sep="\t")
fragcnv_bin = fragcnv_bin[,c(4:9)]
colnames(fragcnv_bin) = c("CNV_bin","chrom","start","end","arm","frag_bin")

P1_feature_ratio$bin <- as.numeric(gsub("[a-zA-Z]","",P1_feature_ratio$variable))
P1_feature_ratio <- left_join(P1_feature_ratio, fragcnv_bin[,c("frag_bin","arm")],c("bin"="frag_bin"))

P1_plot_F5 <- P1_feature_ratio[which(P1_feature_ratio$variable %like% "shortA" ),]
P1_plot_F5 <- P1_plot_F5[rowSums(is.na(P1_plot_F5)) == 0,]
P1_plot_F5$arm <- factor(P1_plot_F5$arm,
                         levels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q",
                                    "10p","10q","11p" ,"11q", "12p", "12q", "13q","14q", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q",
                                    "20p","20q","21q","22q"))
P1_plot_F5$value <- 100 * P1_plot_F5$value

P1_plot_F5_new <- unique(P1_plot_F5) %>%
  dplyr::group_by(PD_event, variable, arm, bin) %>%
  dplyr::summarise(meanValue = mean(value), medianValue = median(value))

F5_panel_C <- ggplot(
  P1_plot_F5_new,
  aes(x=bin, y=medianValue, colour = PD_event)) +
  geom_line(alpha = 0.4, lwd = 1 ) +
  ylab("Median short fragment ratio (%)") + xlab(paste0("Peak 1")) +
  ylim(12,16) +
  labs(colour = "Progressed") + 
  theme_bw() +
  theme(#legend.position = c(0.85,0.2),
    strip.background=element_blank(),
    strip.text.x=element_text(size=11, angle=90),
    strip.text.y=element_text(size=15, angle=0),
    text=element_text(face='bold',size = 12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size=12),
    panel.border = element_blank(),
    legend.key=element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.spacing.x = unit(0.5, "mm"),
    plot.caption = element_text(hjust=0.5, size=rel(1.2)),
    panel.grid.minor = element_blank()) +
  facet_grid(. ~ arm, scales="free_x", space="free_x")

ggsave(filename = "Figure_F5.pdf",
       plot_grid(plot_grid(F5_panel_A, F5_panel_B, ncol = 2, nrow = 1,labels = c("A","B")), F5_panel_C, nrow = 2, labels = c("","C")),
       dpi = 300,width = 15,height = 9)


##### Figure S1

##### combined_model

combined_model <- rbind(P1_model,P2_model)
combined_model$status <- 0
combined_model[combined_model$risk == "High",]$status <- 1
combined_model$joined <- as.numeric(as.character(combined_model$joined))
combined_model$ctDNA <- as.numeric(as.character(combined_model$ctDNA))

combined_model <- combined_model %>%
  group_by(`Patient ID`) %>%
  dplyr::mutate(maxPFS_alter = max(DFS_adjusted)) %>%
  dplyr::mutate(maxStatus = max(status)) %>%
  dplyr::mutate(maxJoined = max(joined)) %>%
  dplyr::mutate(maxCtDNA = max(ctDNA)) %>%
  ungroup() %>% filter(maxPFS_alter == DFS_adjusted)

table(combined_model[,c("PD_event","maxJoined")])
table(combined_model[,c("PD_event","maxCtDNA")])


stats_df <- data.frame(Value = numeric(), 
                       Lower = numeric(), 
                       Upper = numeric(), 
                       type = character(), 
                       label = character(),
                       Group = character())

tmp <- combined_model
tmp$Predict <- 0
tmp[tmp$maxStatus == 1,]$Predict <- 1
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days + 6 months","Fragmentation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days + 6 months","Fragmentation")
stats_df <- rbind(stats_df, tmp_df)

#mutation only
tmp <- combined_model
tmp$Predict <- tmp$maxCtDNA
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days + 6 months","Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days + 6 months","Mutation")
stats_df <- rbind(stats_df, tmp_df)

#with ctDNA
tmp <- combined_model
tmp$Predict <- tmp$maxJoined
tmp$Predict <- factor(tmp$Predict,levels = c(1,0))
tmp$PD_event <- factor(tmp$PD_event, levels = c(1,0))
tmp_stats <- epi.tests(dat = table(tmp[,c("Predict","PD_event")]), conf.level = 0.95)
tmp_df <- stats_df[0,]
tmp_df[1,] <- c(tmp_stats$detail$se$est,tmp_stats$detail$se$lower, tmp_stats$detail$se$upper, "Sensitivity","7 days + 6 months","Fragmentation + Mutation")
tmp_df[2,] <- c(tmp_stats$detail$sp$est,tmp_stats$detail$sp$lower, tmp_stats$detail$sp$upper, "Specificity","7 days + 6 months","Fragmentation + Mutation")
stats_df <- rbind(stats_df, tmp_df)

stats_df$Value <- as.numeric(stats_df$Value)
stats_df$Lower <- as.numeric(stats_df$Lower)
stats_df$Upper <- as.numeric(stats_df$Upper)
stats_df$Group <- factor(stats_df$Group, levels = c("Fragmentation","Mutation","Fragmentation + Mutation"))
# stats_df$label <- factor(stats_df$label, levels = c("7 days","6 months"))

sens_plot_S1 <- ggplot(data=stats_df[stats_df$label != "Combined timepoints",], aes(x=label, y=Value,fill=type)) +
  geom_bar(colour="black",stat="identity", position = position_dodge(0.7), width = 0.5)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,position = position_dodge(0.7)) + 
  ylab("Sensitivity/Sepcificity") + xlab("Postsurgical") +
  facet_wrap(~ Group,ncol=3) + ylim(0,1) + 
  theme(axis.line = element_line(colour = "black", size = 1.25),
        strip.text.x = element_text(size = 13, colour = "black", face = "bold"),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 13, color = "black", face = "bold"),
        legend.text = element_text(size = 13, color = "black", face = "bold"),
        axis.text.x = element_text(size = 13, color = "black", face = "bold"),
        axis.text.y = element_text(size = 13, color = "black", face = "bold"),
        axis.title.x = element_text(size = 13, color = "black", face = "bold"),
        axis.title.y = element_text(size = 13, color = "black", face = "bold"))


tmp <- combined_model
combined_venn_fit <- euler(c(Fragmenation = length(tmp[tmp$maxStatus == 1,]$SampleID) - length(tmp[tmp$maxStatus == 1 & tmp$maxCtDNA == 1,]$SampleID),
                             Mutation = length(tmp[tmp$maxCtDNA == 1,]$SampleID) - length(tmp[tmp$maxStatus == 1 & tmp$maxCtDNA == 1,]$SampleID),
                             "Fragmenation&Mutation" = length(tmp[tmp$maxStatus == 1 & tmp$maxCtDNA == 1,]$SampleID)))
combined_venn_plot <- plot(combined_venn_fit, main = "7 days + 6 months postsurgical",
                           fills = list(fill=c("steelblue4", "darkgoldenrod1"), alpha = 0.5),
                           edges = FALSE, 
                           # labels = list(col = "black", face = "bold", font = 12),
                           legend = TRUE, 
                           adjust_labels =  FALSE,
                           quantities= TRUE)


ggsave(filename = "Figure_S1.pdf",
       plot_grid(combined_venn_plot, sens_plot_S1, nrow = 1, ncol = 2, rel_widths = c(1,2), scale = c(1,0.9), labels = c("A","B")),
       dpi = 300,width = 15,height = 6)


##### Table 1, multivariate and univeriate plots
##P1
tmp <- inner_join(P1_model, clinical_info[,c("SampleID","Age","Sex","SMOKING","StageTnm","Therapy")],c("SampleID"="SampleID"))
tmp$ModelPrediction = "Low risk"
tmp[tmp$risk == "High",]$ModelPrediction = "High risk"
tmp$ModelPrediction <- factor(tmp$ModelPrediction, levels = c("Low risk","High risk"))
tmp$SMOKING <- as.character(tmp$SMOKING)
tmp[tmp$SMOKING == 0,]$SMOKING <- "No"
tmp[tmp$SMOKING == 1,]$SMOKING <- "Yes"
tmp$ctDNA <- as.character(tmp$ctDNA)
tmp[tmp$ctDNA == 1,]$ctDNA = "Pos"
tmp[tmp$ctDNA == 0,]$ctDNA = "Neg"
tmp$combined <- "high"
tmp[tmp$joined == 0,]$combined <- "low"
tmp$combined <- factor(tmp$combined, levels = c("low","high"))

P1_uni_fit <- coxph(Surv(DFS, PD_event) ~ ModelPrediction, data =  tmp)
P1_uni_plot <- ggforest(P1_uni_fit, tmp, fontsize = 1, main = "Hazard ratio\n7 days postsurgical")
P1_multi_fit <- coxph(Surv(DFS, PD_event) ~ ModelPrediction + Age + Sex + SMOKING + Therapy + StageTnm + ctDNA, 
                      data =  tmp)
P1_multi_plot <- ggforest(P1_multi_fit, tmp, fontsize = 1, main = "Hazard ratio\n7 days postsurgical")

P1_uni_fit_Model <- coxph(Surv(DFS, PD_event) ~ ModelPrediction, data =  tmp)
P1_uni_fit_Age <- coxph(Surv(DFS, PD_event) ~ Age, data =  tmp)
P1_uni_fit_Sex <- coxph(Surv(DFS, PD_event) ~ Sex, data =  tmp)
P1_uni_fit_SMOKING <- coxph(Surv(DFS, PD_event) ~ SMOKING, data =  tmp)
P1_uni_fit_Therapy <- coxph(Surv(DFS, PD_event) ~ Therapy, data =  tmp)
P1_uni_fit_StageTnm <- coxph(Surv(DFS, PD_event) ~ StageTnm, data =  tmp)
P1_uni_fit_ctDNA <- coxph(Surv(DFS, PD_event) ~ ctDNA, data =  tmp)
P1_uni_fit_joined <- coxph(Surv(DFS, PD_event) ~ combined, data =  tmp)

summary(P1_uni_fit_Model)
summary(P1_uni_fit_Age)
summary(P1_uni_fit_Sex)
summary(P1_uni_fit_SMOKING)
summary(P1_uni_fit_Therapy)
summary(P1_uni_fit_StageTnm)
summary(P1_uni_fit_ctDNA)
summary(P1_uni_fit_joined)
summary(P1_multi_fit)


##P2
tmp <- inner_join(P2_model, clinical_info[,c("SampleID","Age","Sex","SMOKING","StageTnm","Therapy")],c("SampleID"="SampleID"))
tmp$ModelPrediction = "Low risk"
tmp[tmp$risk == "High",]$ModelPrediction = "High risk"
tmp$ModelPrediction <- factor(tmp$ModelPrediction, levels = c("Low risk","High risk"))
tmp$SMOKING <- as.character(tmp$SMOKING)
tmp[tmp$SMOKING == 0,]$SMOKING <- "No"
tmp[tmp$SMOKING == 1,]$SMOKING <- "Yes"
tmp$ctDNA <- as.character(tmp$ctDNA)
tmp[tmp$ctDNA == 1,]$ctDNA = "Pos"
tmp[tmp$ctDNA == 0,]$ctDNA = "Neg"
tmp$combined <- "high"
tmp[tmp$joined == 0,]$combined <- "low"
tmp$combined <- factor(tmp$combined, levels = c("low","high"))

P2_uni_fit <- coxph(Surv(DFS, PD_event) ~ ModelPrediction, data =  tmp)
P2_uni_plot <- ggforest(P2_uni_fit, tmp, fontsize = 1, main = "Hazard ratio\n6 months postsurgical")
P2_multi_fit <- coxph(Surv(DFS, PD_event) ~ ModelPrediction + Age + Sex + SMOKING + Therapy + StageTnm + ctDNA, 
                      data =  tmp)
P2_multi_plot <- ggforest(P2_multi_fit, tmp, fontsize = 1, main = "Hazard ratio\n6 months postsurgical")
P2_uni_fit_Model <- coxph(Surv(DFS, PD_event) ~ ModelPrediction, data =  tmp)
P2_uni_fit_Age <- coxph(Surv(DFS, PD_event) ~ Age, data =  tmp)
P2_uni_fit_Sex <- coxph(Surv(DFS, PD_event) ~ Sex, data =  tmp)
P2_uni_fit_SMOKING <- coxph(Surv(DFS, PD_event) ~ SMOKING, data =  tmp)
P2_uni_fit_Therapy <- coxph(Surv(DFS, PD_event) ~ Therapy, data =  tmp)
P2_uni_fit_StageTnm <- coxph(Surv(DFS, PD_event) ~ StageTnm, data =  tmp)
P2_uni_fit_ctDNA <- coxph(Surv(DFS, PD_event) ~ ctDNA, data =  tmp)
P2_uni_fit_joined <- coxph(Surv(DFS, PD_event) ~ combined, data =  tmp)

summary(P2_uni_fit_Model)
summary(P2_uni_fit_Age)
summary(P2_uni_fit_Sex)
summary(P2_uni_fit_SMOKING)
summary(P2_uni_fit_Therapy)
summary(P2_uni_fit_StageTnm)
summary(P2_uni_fit_ctDNA)
summary(P2_uni_fit_joined)
summary(P2_multi_fit)


##### Different Validation methods

P1_5_fold <- read.csv(file = "P1_5_Repeat_5_fold_AUC.txt", stringsAsFactors = FALSE)
P1_10_fold <- read.csv(file = "P1_10_Repeat_10_fold_AUC.txt", stringsAsFactors = FALSE)
P1_50_split <- read.csv(file = "P1_50_times_60vs40_random_split_AUCs.txt", stringsAsFactors = FALSE)

P1_5_fold$Method <- "5 fold"
P1_10_fold$Method <- "10 fold"
P1_50_split$Method <- "50X random split"

P1_validation <- rbind(P1_5_fold[,c("AUC","Method")],P1_10_fold[,c("AUC","Method")],P1_50_split[,c("AUC","Method")])
P1_validation$Timepoint <- "TP1"

P2_5_fold <- read.csv(file = "P2_5_Repeat_5_fold_AUC.txt", stringsAsFactors = FALSE)
P2_10_fold <- read.csv(file = "P2_10_Repeat_10_fold_AUC.txt", stringsAsFactors = FALSE)
P2_50_split <- read.csv(file = "P2_50_times_60vs40_random_split_AUC.txt", stringsAsFactors = FALSE)

P2_5_fold$Method <- "5 fold"
P2_10_fold$Method <- "10 fold"
P2_50_split$Method <- "50X random split"

P2_validation <- rbind(P2_5_fold[,c("AUC","Method")],P2_10_fold[,c("AUC","Method")],P2_50_split[,c("AUC","Method")])
P2_validation$Timepoint <- "TP2"

Merged_validation <- rbind(P1_validation,P2_validation)
Merged_validation$Method <- factor(Merged_validation$Method, levels = c("5 fold","10 fold","50X random split"))
Merged_validation$Timepoint <- factor(Merged_validation$Timepoint, levels = c("TP1","TP2"))

figure_S2 <- ggplot(Merged_validation,aes(x=Method,y=AUC, fill = Method)) + geom_boxplot(outlier.shape = NA) +
  xlab("Validation methods") + ylab("AUCs") + theme_bw() +  ylim(c(0, 1)) + facet_grid(. ~ Timepoint) + 
  scale_fill_manual(values = color_code[1:3]) +
  geom_point(aes(fill = Method),size = 1.5, shape = 21, position = position_jitterdodge()) +
  theme(legend.position = c(0.85,0.25),text=element_text(face='bold'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "Figure_S2.pdf",
       figure_S2,
       dpi = 300,width = 8,height = 6)


summary(P1_5_fold[P1_5_fold$Run != "Average",]$AUC)
P1_5_fold[P1_5_fold$Run == "Average",]$AUC

summary(P2_5_fold[P2_5_fold$Run != "Average",]$AUC)
P2_5_fold[P2_5_fold$Run == "Average",]$AUC

summary(P1_10_fold[P1_10_fold$Run != "Average",]$AUC)
P1_10_fold[P1_10_fold$Run == "Average",]$AUC

summary(P2_10_fold[P2_10_fold$Run != "Average",]$AUC)
P2_10_fold[P2_10_fold$Run == "Average",]$AUC


summary(P1_50_split$AUC)
summary(P2_50_split$AUC)



##### Other algorithm LOOCV results

load("Other_model_output.Rdata")


roc.plot.other_P1 <- list()

t <- format(roc(PD_event~Coxnnet_pred,P1_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P1[[paste0("COX-NNET, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~Coxnnet_pred,P1_other_model_output,levels=c('0','1'),
                                                                                    percent=F,smooth=F,ci=T)

t <- format(roc(PD_event~GBM_pred,P1_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P1[[paste0("GBM, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~GBM_pred,P1_other_model_output,levels=c('0','1'),
                                                                         percent=F,smooth=F,ci=T)


t <- format(roc(PD_event~RSF_pred,P1_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P1[[paste0("RSF, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~RSF_pred,P1_other_model_output,levels=c('0','1'),
                                                                    percent=F,smooth=F,ci=T)

t <- format(roc(PD_event~SVM_pred,P1_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P1[[paste0("SVM, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~SVM_pred,P1_other_model_output,levels=c('0','1'),
                                                                    percent=F,smooth=F,ci=T)



figure_S3_panel_P1 <- ggrocs(rocs = roc.plot.other_P1, breaks = seq(0,1,0.2), legendTitel = "7 days postsurgical\nLOOCV\n AUC (95% CI)")


roc.plot.other_P2 <- list()

t <- format(roc(PD_event~Coxnnet_pred,P2_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P2[[paste0("COX-NNET, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~Coxnnet_pred,P2_other_model_output,levels=c('0','1'),
                                                                              percent=F,smooth=F,ci=T)

t <- format(roc(PD_event~GBM_pred,P2_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P2[[paste0("GBM, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~GBM_pred,P2_other_model_output,levels=c('0','1'),
                                                                         percent=F,smooth=F,ci=T)


t <- format(roc(PD_event~RSF_pred,P2_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P2[[paste0("RSF, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~RSF_pred,P2_other_model_output,levels=c('0','1'),
                                                                         percent=F,smooth=F,ci=T)

t <- format(roc(PD_event~SVM_pred,P2_other_model_output,levels=c('0','1'),
                percent=F,smooth=F,ci=T)$ci, digits = 3)
roc.plot.other_P2[[paste0("SVM, ", t[2]," (",t[1],"-",t[3],")")]] <- roc(PD_event~SVM_pred,P2_other_model_output,levels=c('0','1'),
                                                                         percent=F,smooth=F,ci=T)



figure_S3_panel_P2 <- ggrocs(rocs = roc.plot.other_P2, breaks = seq(0,1,0.2), legendTitel = "6 months postsurgical\nLOOCV\n AUC (95% CI)")


ggsave(filename = "Figure_S3.pdf",
       plot_grid(figure_S3_panel_P1, figure_S3_panel_P2, ncol= 2, labels = c("A","B")),
       dpi = 300,width = 12,height = 6)

##### RFE 

P1_RFE <- read.csv(file = "P1_RFE_LOOCV_AUC.csv", stringsAsFactors = FALSE)
P2_RFE <- read.csv(file = "P2_RFE_LOOCV_AUC.csv", stringsAsFactors = FALSE)

P1_RFE$Feature_number <- P1_RFE$Feature_number + 1
P2_RFE$Feature_number <- P2_RFE$Feature_number + 1

P1_length <- length(P1_RFE)
P2_length <- length(P2_RFE)

P1_RFE_plot <- ggplot(data = P1_RFE, aes(x = Feature_number, y = LOOCV.AUC, group = 1)) + geom_line() + geom_point(color = color_code[1]) + 
  ggtitle("7 days postsurgical") + ylim(0,1) + 
  xlab("Top features") + ylab("LOOCV AUC") + 
  theme_bw() +  
  theme(legend.position = c(0.85,0.25),text=element_text(face='bold'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


P2_RFE_plot <- ggplot(data = P2_RFE, aes(x = Feature_number, y = LOOCV_AUC, group = 1)) + geom_line() + geom_point(color = color_code[2]) + 
  ggtitle("6 months postsurgical") + ylim(0,1) + 
  xlab("Top features") + ylab("LOOCV AUC") + 
  theme_bw() + 
  theme(legend.position = c(0.85,0.25),text=element_text(face='bold'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggsave(filename = "Figure_S4.pdf",
       plot_grid(P1_RFE_plot, P2_RFE_plot, ncol= 2, labels = c("A","B")),
       dpi = 300,width = 12,height = 4)



