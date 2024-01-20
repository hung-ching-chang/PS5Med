library(ggplot2)
library(ggpattern)
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/simulation/May2023final/simulation_Oct_conti_exposure/power/")
# arrange multiple figures in one plot

library(ggpubr)
for(cor.setting in c(0,0.5,0.9)){
  for(num.ts in c(5,10,30)){
    # xmax <- ifelse(num.ts <=5, 0.5, 0.2)
    xmax <- 0.2
    file.name <- paste0("cor",cor.setting, "ts",num.ts)
    
    load(paste0("power_N500p1000_", file.name, ".RData"))
    print(paste0(file.name, ": ", length(power.result)))
    
    #### visualization 1 - power #######
    effect.a.setting <- seq(0, 0.2, 0.01)
    effect.a.setting <- effect.a.setting[-1]
    power.df <- data.frame(matrix(ncol = 4, nrow = 20))
    colnames(power.df) <- c("PS5", "HP2016", "HIMA", "HILMA")
    for(i in 1:20){
      power.df[i,] <- power.result[[i]]$power[1:4]
    }
    power.plot <- data.frame(alpha.effect = rep(effect.a.setting, 4),
                             methods = rep(c("PS5", "H&P", "HIMA", "HILMA"), 
                                           each = length(effect.a.setting)),
                             power = unlist(power.df))
    power.plot$methods <- factor(power.plot$methods, levels = c("H&P", "HILMA",
                                                                "HIMA", "PS5"))
    p <- ggplot(power.plot, aes(x=alpha.effect, y=power)) +
      geom_line(aes(color=methods), linewidth = 2) +
      scale_color_manual(values = c("#eecc16", "#0000a7", "#008176", "#c1272d")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, xmax)) +
      xlab("") +
      ylab("") +
      theme(legend.position="top",
            legend.title = element_text(size=20),
            legend.text=element_text(size=20),
            legend.key.height= unit(2, 'cm'),
            legend.key.width= unit(2, 'cm'),
            axis.title = element_text(size = 20),
            axis.text=element_text(size=18))
    if(num.ts == 5){p <- p + ylab("Power")}
    if(cor.setting == 0.9){p <- p + xlab("Alpha")}
    assign(paste0("power_", file.name), p)

    #### visualization 2 - ME bias #######
    effect.a.setting <- seq(0, 0.4, 0.01)
    effect.a.setting <- effect.a.setting[-1]
    ME.bias.df <- data.frame(matrix(ncol = 4, nrow = 40))
    colnames(ME.bias.df) <- c("PS5", "H&P", "HIMA", "HILMA")
    for(i in 1:40){
      ME.bias.df[i,] <- power.result[[i]]$ME.bias
    }
    # for visualization, fix the maximum value as 100%
    ME.bias.df[ME.bias.df > 1] <- 1
    ME.bias.plot <- data.frame(alpha.effect = rep(effect.a.setting, 4),
                               methods = rep(c("PS5", "H&P", "HIMA", "HILMA"), 
                                             each = length(effect.a.setting)),
                               ME.bias = unlist(ME.bias.df))
    ME.bias.plot$methods <- factor(ME.bias.plot$methods,
                                   levels = c("H&P", "HILMA",
                                              "HIMA", "PS5"))
    p <- ggplot(ME.bias.plot, aes(x=alpha.effect, y=ME.bias)) +
      geom_line(aes(color=methods), linewidth = 2) +
      scale_color_manual(values = c("#eecc16", "#0000a7", "#008176", "#c1272d")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 0.4)) +
      xlab("")  + 
      ylab("") +
      theme(legend.position="top",
            legend.title = element_text(size=20),
            legend.text=element_text(size=20),
            legend.key.height= unit(2, 'cm'),
            legend.key.width= unit(2, 'cm'),
            axis.title.x = element_text(size = 20),
            axis.text=element_text(size=18))
    if(num.ts == 5){p <- p + ylab("Percent relative bias")}
    if(cor.setting == 0.9){p <- p + xlab("Alpha")}
    assign(paste0("MEbias_", file.name), p)
    
    #### visualization 3 - Sensitivity #######
    max.top.k <- ifelse(num.ts <= 10, 20, 50)
    for(effect.a in c(0.1,0.2,0.3)){
      idx <- which(effect.a.setting == effect.a)
      PS5.pval <- power.result[[idx]]$PS5.contri.pval
      HIMA.pval <- power.result[[idx]]$HIMA.contri.pval
      sen.PS5 <- sen.HIMA <- rep(NA, max.top.k)
      for(top.k in 1:max.top.k){
        # select top k mediators
        PS5.topk <- apply(PS5.pval, 1,
                          function(x) rank(x, ties.method = "last") <= top.k)
        HIMA.topk <- apply(HIMA.pval, 1,
                           function(x) rank(x, ties.method = "last") <= top.k)
        # sensitivity
        sen.PS5[top.k] <- mean(PS5.topk[1:num.ts,])
        sen.HIMA[top.k] <- mean(HIMA.topk[1:num.ts,])
      }
      sen.df <- data.frame(top.k = rep(1:max.top.k, 2),
                           sensitvity = c(sen.PS5, sen.HIMA),
                           methods = rep(c("PS5","HIMA"), each =max.top.k))
      sen.df$methods <- factor(sen.df$methods, 
                               levels = c("HIMA", "PS5"))
      p <- ggplot(sen.df, aes(x=top.k, y=sensitvity)) +
        geom_line(aes(color=methods, linetype=methods), linewidth = 2) +
        scale_color_manual(values = c("#008176", "#c1272d")) +
        scale_linetype_manual(values=c("dashed", "solid"))+
        scale_y_continuous(limits = c(0, 1)) +
        xlab("")  + 
        ylab("") +
        theme(legend.position="top",
              legend.title = element_text(size=20),
              legend.text=element_text(size=20),
              legend.key.height= unit(2, 'cm'),
              legend.key.width= unit(2, 'cm'),
              axis.title = element_text(size = 20),
              axis.text=element_text(size=18))
      if(effect.a == 0.1){p <- p + ylab("Sensitivity")}
      if(cor.setting == 0.9){p <- p + xlab("Top k mediators")}
      assign(paste0("sen_alpha", effect.a, "_", file.name), p)
    }
  }
}
ggarrange(power_cor0ts5, power_cor0ts10, power_cor0ts30,
          power_cor0.5ts5, power_cor0.5ts10, power_cor0.5ts30,
          power_cor0.9ts5, power_cor0.9ts10, power_cor0.9ts30,
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "bottom")
ggsave("power_all_truncated.png", width = 24, height = 12, dpi = 1200)
ggarrange(MEbias_cor0ts5, MEbias_cor0ts10, MEbias_cor0ts30,
          MEbias_cor0.5ts5, MEbias_cor0.5ts10, MEbias_cor0.5ts30,
          MEbias_cor0.9ts5, MEbias_cor0.9ts10, MEbias_cor0.9ts30,
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "bottom")
ggsave("MEbias_all.png", width = 24, height = 12, dpi = 1200)

ggarrange(sen_alpha0.1_cor0ts5, sen_alpha0.2_cor0ts5, sen_alpha0.3_cor0ts5,
          sen_alpha0.1_cor0.5ts5, sen_alpha0.2_cor0.5ts5, sen_alpha0.3_cor0.5ts5,
          sen_alpha0.1_cor0.9ts5, sen_alpha0.2_cor0.9ts5, sen_alpha0.3_cor0.9ts5,
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "bottom")
ggsave("sen_ts5_final.png", width = 24, height = 12, dpi = 1200)

ggarrange(sen_alpha0.1_cor0ts10, sen_alpha0.2_cor0ts10, sen_alpha0.3_cor0ts10,
          sen_alpha0.1_cor0.5ts10, sen_alpha0.2_cor0.5ts10, sen_alpha0.3_cor0.5ts10,
          sen_alpha0.1_cor0.9ts10, sen_alpha0.2_cor0.9ts10, sen_alpha0.3_cor0.9ts10,
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "bottom")
ggsave("sen_ts10_final.png", width = 24, height = 12, dpi = 1200)

ggarrange(sen_alpha0.1_cor0ts30, sen_alpha0.2_cor0ts30, sen_alpha0.3_cor0ts30,
          sen_alpha0.1_cor0.5ts30, sen_alpha0.2_cor0.5ts30, sen_alpha0.3_cor0.5ts30,
          sen_alpha0.1_cor0.9ts30, sen_alpha0.2_cor0.9ts30, sen_alpha0.3_cor0.9ts30,
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "bottom")
ggsave("sen_ts30_final.png", width = 24, height = 12, dpi = 1200)



#### supplementary ####
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/simulation/May2023final/simulation_Oct_conti_exposure/gamma/")

for(num.ts in c(1,5,10,30,50)){
  alpha.max <- ifelse(num.ts<=10, 0.4, 0.2)
  file.name <- paste0("ts",num.ts)
  load(paste0("gamma_N500p1000_", file.name, ".RData"))
  effect.a.setting <- seq(0, alpha.max, length.out = 21)
  effect.a.setting <- effect.a.setting[-1]
  #### visualization #######
  power.plot <- data.frame(alpha.effect = rep(effect.a.setting, 2),
                           methods = rep(c("PS5.L1", "PS5.L2"), 
                                         each = length(effect.a.setting)),
                           power = unlist(gamma.result))
  power.plot$methods <- factor(power.plot$methods, levels = c("PS5.L1", "PS5.L2"))
  p <- ggplot(power.plot, aes(x=alpha.effect, y=power)) +
    geom_line(aes(color=methods), linewidth = 2) +
    scale_color_manual(values = c("black", "#c1272d")) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, xmax)) +
    xlab("Alpha") +
    ylab("Power") +
    theme(legend.position="top",
          legend.title = element_text(size=20),
          legend.text=element_text(size=20),
          legend.key.height= unit(2, 'cm'),
          legend.key.width= unit(2, 'cm'),
          axis.title = element_text(size = 20),
          axis.text=element_text(size=14))
  
  tiff(paste0("gamma_", file.name), 
      width = 1440, height = 960, res = 200)
  print(p)
  dev.off()
}





#####################
######## END ########
#####################



####       single plot       #######
#### visualization 1 - power #######
power.df <- data.frame(matrix(ncol = 4, nrow = 20))
colnames(power.df) <- c("PS5", "HP2016", "HIMA", "HILMA")
for(i in 1:length(power.result)){
  power.df[i,] <- power.result[[i]]$power
}
power.plot <- data.frame(alpha.effect = rep(effect.a.setting, 4),
                         methods = rep(c("PS5", "HP2016", "HIMA.minp", "HILMA"), 
                                       each = length(effect.a.setting)),
                         power = unlist(power.df))
power.plot$methods <- factor(power.plot$methods, levels = c("PS5", "HP2016", 
                                                            "HIMA.minp", "HILMA"))
ggplot(power.plot, aes(x=alpha.effect, y=power)) +
  geom_line(aes(color=methods)) +
  scale_color_manual(values = c("#c1272d", "#008176", "#eecc16", "#0000a7")) +
  #scale_linetype_manual(values=c("solid", "dotted", "solid", "dotted", "solid", "solid", "solid")) +, "#0cf371", "#000000"
  scale_y_continuous(limits = c(0, 1)) +
  xlab("alpha") +
  theme(legend.position="top")
ggsave(paste0("power_plot/", file.name, ".png"),
       width = 6, height = 4, dpi = 1200)

#### visualization 2 - ME bias #######
ME.bias.df <- data.frame(matrix(0, ncol = 3, nrow = 20))
colnames(ME.bias.df) <- c("PS5", "HIMA", "HILMA")
for(i in 1:length(power.result)){
  ME.bias.df[i,] <- power.result[[i]]$ME.bias[c("PS5", "HIMA", "HILMA")]
}
# for visualization, fix the maximum value as 100%
ME.bias.df[ME.bias.df > 1] <- 1
ME.bias.plot <- data.frame(global.ME = rep(effect.a.setting, 3),
                           methods = rep(c("PS5", "HIMA", "HILMA"), 
                                         each = length(effect.a.setting)),
                           ME.bias = unlist(ME.bias.df))
ME.bias.plot$methods <- factor(ME.bias.plot$methods,
                               levels = c("PS5", "HIMA", "HILMA"))
ggplot(ME.bias.plot, aes(x=global.ME, y=ME.bias, fill=methods)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("#c1272d", "#eecc16", "#0000a7")) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Alpha")  + 
  ylab("Percent relative bias") 
ggsave(paste0("ME_bias/", file.name, ".png"),
       width = 6, height = 4, dpi = 1200)

#### visualization 3 - sensitivity #######
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/simulation/May2023final/simulation/power/")
num.ts <- 1
cor <- 0.5
file.name <- paste0("cor",cor, "ts",num.ts)
load(paste0("power_N500p1000_", file.name, ".RData"))
threshold.setting <- seq(0, 1, 0.05)
threshold.setting <- threshold.setting[-1]
# show the result of alpha 0.25
if(num.ts <= 10){
  effect.a.setting <- seq(0, 0.5, 0.025)
  effect.a.setting <- effect.a.setting[-1]
}else{
  effect.a.setting <- seq(0, 0.15, 0.01)
  effect.a.setting <- effect.a.setting[-1]
}
sen.df <- data.frame(matrix(ncol = 4, nrow = 0))
for(j in 1:length(effect.a.setting)){
  PS5.contri.pval <- unlist(power.result[[j]]$PS5.contri.pval[,1:num.ts])
  HIMA.contri.pval <- unlist(power.result[[j]]$HIMA.contri.pval[,1:num.ts])
  for(i in 1:length(threshold.setting)){
    sen.df <- rbind(sen.df,
                    c(pval.threhold = threshold.setting[i],
                      effect.a = effect.a.setting[j],
                      sen.ps5 = mean(PS5.contri.pval < threshold.setting[i]),
                      sen.hima = mean(HIMA.contri.pval < threshold.setting[i])))
  }
}
colnames(sen.df) <- c("pval.threhold", "Alpha", "PS5", "HIMA")
sen.df$diff <- sen.df$PS5 - sen.df$HIMA
sen.df$higher.performance <- ifelse(sen.df$diff >= 0, "PS5", "HIMA")
table(sen.df$higher.performance)

sen.df$higher.performance <- factor(sen.df$higher.performance, 
                                    levels = c("HIMA", "PS5"))
ggplot(sen.df, aes(x=Alpha, y=pval.threhold, color=higher.performance, size=abs(diff))) +
  geom_point() +
  scale_color_manual(values = c("#008176", "#c1272d")) +
  scale_x_continuous(limits = c(0, max(effect.a.setting))) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Alpha") +
  ylab("p-value threshold") +
  theme(legend.position="top")
#ggsave(paste0("sen_plot/", file.name, ".png"),
#       width = 6, height = 4, dpi = 1200)



############ end ###################
# MSE plot
# MSE.plot <- data.frame(alpha.effect = rep(effect.a.setting, 2),
#                        methods = rep(c("PECMAN", "HIMA"), 
#                                      each = length(effect.a.setting)),
#                        MSE = unlist(TS.MSE.table))
# MSE.plot$methods <- factor(MSE.plot$methods, levels = c("PECMAN", "HIMA"))
# ggplot(MSE.plot, aes(x=alpha.effect, y=MSE, fill=methods)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   scale_fill_manual(values = c("#045275", "#F0746E")) +
#   xlab("alpha")
# ggsave("tsMSE.png", width = 8, height = 4, dpi = 1200)
# 
# # True positive rate
# Youden.plot <- data.frame(alpha.effect = rep(effect.a.setting, 2),
#                           methods = rep(c("PECMAN", "HIMA"), 
#                                         each = length(effect.a.setting)),
#                           Youden = unlist(Youden.table))
# Youden.plot$methods <- factor(Youden.plot$methods, levels = c("PECMAN", "HIMA"))
# ggplot(Youden.plot, aes(x=alpha.effect, y=Youden, fill=methods, pattern=methods)) +
#   geom_bar_pattern(stat="identity", position=position_dodge(),
#                    pattern_color = "white",
#                    pattern_angle = 45,
#                    pattern_density = 0.01,
#                    pattern_spacing = 0.005,
#                    pattern_alpha = 0.5) +
#   scale_pattern_manual(values = c("none", "none")) +
#   scale_fill_manual(values = c("#0000a7", "#eecc16")) +
#   scale_y_continuous(limits = c(0, 1))
# ggsave("true positive rate.png", width = 8, height = 4, dpi = 1200)







######################################################################
########################## true signal MSE ########################### 
######################################################################
# true.ME <- Alpha.x* Beta.m
# signal <- true.ME != 0
# HIMA.me <- PECMAN.me <- rep(0, p)
# if(!is.null(HIMA.result)){
#   HIMA.selected.mediator <- sapply(rownames(HIMA.result), function(x) which(colnames(mediator) == x))
#   HIMA.me[HIMA.selected.mediator] <- HIMA.result$`alpha*beta`
# }
# if(sum(signal) == 0){
#   PECMAN.ts <- HIMA.ts <- 0
# }else if(any(is.na(PECMAN.result$estimation))){
#   PECMAN.ts <- 0
#   HIMA.ts <- sum((HIMA.me[signal] - true.ME[signal])^2)/length(signal)
# }else{
#   PECMAN.selected.mediator <- sapply(rownames(PECMAN.result$estimation), function(x) which(colnames(mediator) == x))
#   PECMAN.me[PECMAN.selected.mediator] <- PECMAN.result$estimation[,3]
#   PECMAN.ts <- sum((PECMAN.me[signal] - true.ME[signal])^2)/length(signal)
#   HIMA.ts <- sum((HIMA.me[signal] - true.ME[signal])^2)/length(signal)
# }
# TS.MSE[rep.time,] <- c(PECMAN.ts, HIMA.ts)

######################################################################
####################### Sensitivity (TP/TP+FN) ####################### 
######################################################################
# HIMA.pos <- PECMAN.pos <- rep(F, p)
# if(!is.null(HIMA.result)){
#   HIMA.sig.p <- HIMA.result[HIMA.result$Bonferroni.p < 0.05,]
#   HIMA.sig.idx <- sapply(rownames(HIMA.sig.p), function(x) which(colnames(mediator) == x))
#   if(length(HIMA.sig.idx) > 0){HIMA.pos[HIMA.sig.idx] <- T}
# }
# if(sum(signal) == 0){
#   PECMAN.sen <- HIMA.sen <- 0
#   PECMAN.spe <- sum(!PECMAN.pos & !signal)/sum(!signal)
#   HIMA.spe <- sum(!HIMA.pos & !signal)/sum(!signal)
# }else{
#   if(!any(is.na(PECMAN.result$estimation))) PECMAN.pos[PECMAN.selected.mediator] <- (PECMAN.result$estimation[,5] < 0.05)
#   PECMAN.sen <- sum(PECMAN.pos & signal)/sum(signal)
#   HIMA.sen <- sum(HIMA.pos & signal)/sum(signal)
#   PECMAN.spe <- sum(!PECMAN.pos & !signal)/sum(!signal)
#   HIMA.spe <- sum(!HIMA.pos & !signal)/sum(!signal)
# }
# PECMAN.Youden <- PECMAN.sen + PECMAN.spe - 1
# HIMA.Youden <- HIMA.sen + HIMA.spe - 1
# #SEN[rep.time,] <- c(PECMAN.sen, HIMA.sen)
# Youden[rep.time,] <- c(PECMAN.Youden, HIMA.Youden)

for(i in 1:multi.split){
  seed <- 1000+i
  PS5.result <- PS5.mediation(M = mediator, 
                              X = exposure.mat, 
                              Y = diz, 
                              C = conf.mat,
                              n.draw = 10000,
                              seed = seed,
                              dim.reduction = T)
  PS5.global.me[i] <- PS5.result$global.me
  selected.mediators <- which(paste0("M", c(1:p)) %in% row.names(PS5.result$estimation))
  PS5.med.contri.pval[i,selected.mediators] <- PS5.result$estimation$p.value.Bonf
}
