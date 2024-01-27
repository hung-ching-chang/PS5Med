library(ggplot2)
library(ggpattern)
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
    p <- ggplot(power.plot, aes(x = alpha.effect, y = power)) +
      geom_line(aes(color = methods), linewidth = 2) +
      scale_color_manual(values = c("#eecc16", "#0000a7", "#008176", "#c1272d")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, xmax)) +
      xlab("") +
      ylab("") +
      theme(legend.position = "top",
            legend.title = element_text(size=20),
            legend.text = element_text(size=20),
            legend.key.height = unit(2, 'cm'),
            legend.key.width = unit(2, 'cm'),
            axis.title = element_text(size = 20),
            axis.text = element_text(size=18))
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

