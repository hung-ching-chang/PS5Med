library(Matrix)
library(mvtnorm)
library(glmnet)
library(HIMA)
library(ggplot2)
library(cowplot)
library(parallel)
## COPDGene dataset
## What next?
## (1) try smoking as exposure
########################################
### PCA (up to 80%)
########################################
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/")
load("imaging_581_128_new.RData")
pca.extraction <- function(x, pc.num = 5){
  pca.result <- prcomp(x,
                       center = TRUE,
                       scale. = TRUE)
  cum.var <- cumsum(pca.result$sdev^2 / sum(pca.result$sdev^2))
  return(pca.result$x[,1:pc.num])
}
imaging.all.pc <- matrix(NA, nrow = nrow(imaging.all), ncol = 581*5)
for(i in 1:581){
  imaging.all.pc[,5*(i-1)+(1:5)] <- pca.extraction(imaging.all[,(i-1)*128+(1:128)],
                                                   pc.num = 5)
}

#save(imaging.all.pc, file = "context_matter_patch_pc5.RData")

########################################
### 2022/01/12 (PRS -> imaging -> distance walked)
########################################
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/PRS/")
subtype.df <- read.csv("COPDGene_Subtype_Smoker.csv")  # exposure
PRS <- subtype.df[,c("sid" , "prs_composite")]
#top5 <- PRS$prs_composite > quantile(PRS$prs_composite, prob = 0.7, na.rm = T)
#bottom5 <- PRS$prs_composite < quantile(PRS$prs_composite, prob = 0.3, na.rm = T)
#PRS$prs_composite[top5] <- 1
#PRS$prs_composite[bottom5] <- 0
#PRS <- PRS[(top5 | bottom5),]
load("context_matter_patch_pc10.RData")  # mediators
image.id.df <- read.csv("context_matter_id.csv", header = F)
image.id <- apply(image.id.df, 1, function(x) strsplit(x, "_")[[1]][1])
#longitudinal.df <- read.csv("COPDGene_P1P2P3_Flat_SM_NS_Mar20.txt",
#                            sep = '\t', na.strings = c("",".","NA"))  # outcome
clinic.df <- read.csv("clinical_surv_df.txt", sep = " ")
FEV1pp <- clinic.df[,c("sid", "FEV1pp_utah")]
#dist.walked <- longitudinal.df[,c("sid", "distwalked_P1")]
confounding.df <- read.csv("covariate_df.txt", sep = " ") #confounder
confounding.df <- confounding.df[,c("IID", "ATS_PackYears", "gender", "Height_CM", "Age_Enroll")]
# match id
match.id <- Reduce(intersect, list(PRS$sid, image.id, FEV1pp$sid, confounding.df$IID))
PRS <- PRS[sapply(match.id, function(x) which(PRS$sid == x)),2]
imaging.dat <- imaging.all.pc[sapply(match.id, function(x) which(image.id == x)),]
FEV1pp <- FEV1pp[sapply(match.id, function(x) which(FEV1pp$sid == x)),2]
confounding.all <- confounding.df[sapply(match.id, function(x) which(confounding.df$IID == x)), -1]
confounding.all <- as.matrix(confounding.all)
COPD.PRS.dataset <- list(PRS, imaging.dat, FEV1pp, confounding.all)
#save(COPD.PRS.dataset, file = "Dataset_for_journal.RData")

#####################start from HERE!##################
load("Dataset_for_journal.RData")
PRS <- COPD.PRS.dataset[[1]]
imaging.dat <- COPD.PRS.dataset[[2]]
FEV1pp <- COPD.PRS.dataset[[3]]
confounding.all <- COPD.PRS.dataset[[4]]

# remove missing samples
missing.samples <- is.na(PRS) | is.na(FEV1pp)
exposure <- PRS[!missing.samples]
mediators <- imaging.dat[!missing.samples,]
colnames(mediators) <- paste0(rep(paste0("M", 1:581), each = 10), "-", 1:10)
outcome <- FEV1pp[!missing.samples]
confounders <- confounding.all[!missing.samples,]

# total effect
plot(exposure, outcome)
total.effect.model <- lm(outcome ~ exposure + confounders)
summary(total.effect.model)
# multi-split
multi.split <- 500
#PS5.estimation <- PS5.pvalue <- NULL
# multi sample split
str.time <- Sys.time()
PS5.result <- PS5.multi.split(M = mediators, 
                              X = exposure, 
                              Y = outcome, 
                              C = confounders,
                              n.draw = 10000,
                              dim.reduction = T,
                              multi.num = multi.split,
                              cores = detectCores()-1)
Sys.time() - str.time
save(PS5.result, file = "result/PRS_result10k.RData")
global.IE.test <- PS5.result$global.test
global.IE <- PS5.result$global.me
global.IE.prop <- PS5.result$global.me.prop
med.contri <- PS5.result$mediation.contri
cancel.out.prop <- 1 - global.IE.prop[1]/sum(abs(global.IE.prop[2:3]))


########## smoking by packyear #################
clinic.df <- read.csv("clinical_surv_df.txt", sep = " ")
FEV1pp <- clinic.df[,c("sid", "FEV1pp_utah")]   # outcome
#dist.walked <- longitudinal.df[,c("sid", "distwalked_P1")]
confounding.df <- read.csv("covariate_df.txt", sep = " ") #confounder
confounding.df <- confounding.df[,c("IID", "ATS_PackYears", "gender", "Height_CM", "Age_Enroll")]
# match id
match.id <- Reduce(intersect, list(image.id, FEV1pp$sid, confounding.df$IID))
imaging.dat <- imaging.all.pc[sapply(match.id, function(x) which(image.id == x)),]
FEV1pp <- FEV1pp[sapply(match.id, function(x) which(FEV1pp$sid == x)),2]
confounding.all <- confounding.df[sapply(match.id, function(x) which(confounding.df$IID == x)), -1]
confounding.all <- as.matrix(confounding.all)
PackYears <- confounding.all[,1]
confounding.all <- confounding.all[,-1]
COPD.smoking.dataset <- list(PackYears, imaging.dat, FEV1pp, confounding.all)

## prepare input
exposure <- PackYears
mediators <- imaging.dat
colnames(mediators) <- paste0(rep(paste0("M", 1:581), each = 10), "-", 1:10)
outcome <- FEV1pp
confounders <- confounding.all

# multi-split
multi.split <- 500
#PS5.estimation <- PS5.pvalue <- NULL
# multi sample split
str.time <- Sys.time()
PS5.result <- PS5.multi.split(M = mediators, 
                              X = exposure, 
                              Y = outcome, 
                              C = confounders,
                              n.draw = 10000,
                              dim.reduction = T,
                              multi.num = multi.split)
Sys.time() - str.time
save(PS5.result, file = "../PackYears/result/PackYears_result10k.RData")
global.IE.test <- PS5.result$global.test
global.IE <- PS5.result$global.me
global.IE.prop <- PS5.result$global.me.prop
med.contri <- PS5.result$mediation.contri
cancel.out.prop <- 1 - global.IE.prop[1]/sum(abs(global.IE.prop[2:3]))




##########visualization via multiple figures#################
# keep top 50 mediators with higher p-value
top50 <- med.contri[1:50,]
top50$mediatorID <- row.names(top50)
top50$mediatorID <- factor(top50$mediatorID, levels = row.names(top50))
top50$significance <- (top50$pval.BY < 0.05)
bar.prop.df <- top50
bar.prop.plt <- ggplot(bar.prop.df, aes(x = mediatorID, y = preselected.prop, fill = significance)) +
  geom_bar(stat="identity") +
  scale_x_discrete(drop=F) +
  scale_y_continuous(name="Pre-selected proportions", limits = c(0,1)) +
  scale_fill_manual(values=c("grey", "#31a354")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
bar.pvalue.df <- top50
bar.pvalue.plt <- ggplot(bar.pvalue.df, aes(x = mediatorID, y = -log10(pval.BY), fill = significance)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop=F, name="Top 50 mediators' ID") +
  scale_y_continuous(name="-log10 p-value") +
  scale_fill_manual(values=c("grey", "#31a354")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
# mix plots
png(file = "result/each_mediation_contribution.png", 
    height = 1300, width = 2200,res = 300)
cowplot::plot_grid(bar.prop.plt, bar.pvalue.plt,
                   ncol = 1, rel_heights = c(2, 2),
                   align = 'v', axis = 'lr')
dev.off()

bar.global.IE.df <- data.frame(direction = factor(c("positive", "negative"), c("positive", "negative")),
                               abs.contri = abs(global.IE.prop[2:3]),
                               contri = global.IE.prop[2:3])
bar.global.IE.plt <- ggplot(bar.global.IE.df, aes(x = abs.contri, y = direction, fill = direction)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=paste0(round(contri,1), "%")), 
            hjust=-0.1, color="black", size=4) +
  scale_x_continuous(name="contribution (%)", 
                     limits = c(0,max(bar.global.IE.df$abs.contri+15)))
png(file = "result/total_mediation_contribution.png", 
    height = 500, width = 1000,res = 300)
bar.global.IE.plt
dev.off()

####### 3D visualization #########
# plot pre-selected proportion and aggregated p-value on 3D image
#selected.proportion.all <- rep(0, 581)
#for(i in 1:multi.split){
#  patch.level.name <- as.numeric(sapply(row.names(PS5.result$detail.estimation[[i]]), 
#                                        function(x) strsplit(x, split = "M|-")[[1]][2]))
#  selected.proportion.all[patch.level.name] <- selected.proportion.all[patch.level.name]+1
#}
#selected.proportion.all <- selected.proportion.all/multi.split
#write.table(selected.proportion.all, file = "result/3D_visualization/PS5/preselected_proportion.txt",
#            row.names = F, col.names = F)

# result for PRS and PackYears
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/")
load("PRS/final_result/PRS_result10k.RData")
load("PackYears/final_result/PackYears_result10k.RData")
med.contri <- PS5.result$mediation.contri
# aggregated p-value
patch.level.name <- as.numeric(sapply(row.names(med.contri), 
                                      function(x) strsplit(x, split = "M|-")[[1]][2]))
pval.all <- rep(0, 581)
unique.gene.id <- unique(patch.level.name)
for(i in unique.gene.id){
  idx <- which(patch.level.name %in% i)
  pval.all[i] <- -log10(min(med.contri$pval[idx]))
}

write.table(pval.all, file = "PRS/final_result/pvalue_all.txt",
            row.names = F, col.names = F)

#######################
###### Figure #######
#######################
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/")
patch.loc <- read.table("patch_loc.txt", sep = ",")

# result for PRS and PackYears
load("PRS/final_result/PRS_result10k.RData")
PRS.res <- PS5.result$mediation.contri
load("PackYears/final_result/PackYears_result10k.RData")
PY.res <- PS5.result$mediation.contri

# selected patch
PRS.patch <- row.names(PRS.res)[PRS.res$pval < 0.01]
PY.patch <- row.names(PY.res)[PY.res$pval < 0.01]
PRS.patch.idx <- sapply(PRS.patch, function(x) unlist(strsplit(x, split = "-"))[1])
PRS.patch.idx <- sapply(PRS.patch.idx, function(x) unlist(strsplit(x, split = "M"))[2])
PRS.patch.idx <- as.numeric(PRS.patch.idx)
PRS.patch.idx <- unique(PRS.patch.idx)

PY.patch.idx <- sapply(PY.patch, function(x) unlist(strsplit(x, split = "-"))[1])
PY.patch.idx <- sapply(PY.patch.idx, function(x) unlist(strsplit(x, split = "M"))[2])
PY.patch.idx <- as.numeric(PY.patch.idx)
PY.patch.idx <- unique(PY.patch.idx)

length(PRS.patch.idx);length(PY.patch.idx);length(intersect(PRS.patch.idx, PY.patch.idx))

# venn diagram 1
library(ggvenn)
venn.list <- list(PRS = PRS.patch,
                  `Pack Years` = PY.patch)

ggvenn(
  venn.list, 
  fill_color = c("cyan3", "darkgoldenrod2"),
  stroke_size = 0.5, set_name_size = 4
)


# Table
union.patch <- union(PRS.patch, PY.patch)
union.patch.PRS <- PRS.res[union.patch,]
union.patch.PY <- PY.res[union.patch,]
library(metap)
fisher.pval <- sapply(c(1:length(union.patch)), 
                      function(i) sumlog(c(union.patch.PRS$pval[i], 
                                           union.patch.PY$pval[i]))$p)

union.patch.table <- data.frame(PRS.pval = union.patch.PRS$pval,
                                PY.pval = union.patch.PY$pval,
                                PRS.med.contri = union.patch.PRS$mediation.contribution,
                                PY.med.contri = union.patch.PY$mediation.contribution,
                                PRS.contri.prop = union.patch.PRS$contribution.proportion,
                                PY.contri.prop = union.patch.PY$contribution.proportion,
                                Fisher.combined.pval = fisher.pval)
row.names(union.patch.table) <- row.names(union.patch.PRS)
write.csv(union.patch.table, file = "union.patch.table.csv")

# prepare for 3D figure
patch.level.name <- as.numeric(sapply(row.names(union.patch.table), 
                                      function(x) strsplit(x, split = "M|-")[[1]][2]))
union.patch.table$idx <- patch.level.name
pval.all <- rep(0, 581)
unique.patch.id <- unique(patch.level.name)
for(i in unique.patch.id){
  PRS.pval <- min(union.patch.table[union.patch.table$idx == i,"PRS.pval"])
  PY.pval <- min(union.patch.table[union.patch.table$idx == i,"PY.pval"])
  if(PRS.pval < 0.01 & PY.pval < 0.01){
    pval.all[i] <- 3   # both = 3
  }else if(PRS.pval < 0.01){
    pval.all[i] <- 2   #PRS = 2
  }else if(PY.pval < 0.01){
    pval.all[i] <- 1   #PY = 1
  }
}
write.table(pval.all, file = "detected_patch.txt",
            row.names = F, col.names = F)

# histogram - x 
overlap.idx <- intersect(PY.patch.idx, PRS.patch.idx)
overlap.patch.loc <- patch.loc[overlap.idx,]
PY.patch.loc <- patch.loc[PY.patch.idx[!(PY.patch.idx %in% overlap.idx)],]
PRS.patch.loc <- patch.loc[PRS.patch.idx[!(PRS.patch.idx %in% overlap.idx)],]
overlap.patch.loc.v1 <- data.frame(table(factor(overlap.patch.loc$V1, 
                                           levels = sort(unique(patch.loc$V1)))), 
                                   exposure = "Both")
PY.patch.loc.v1 <- data.frame(table(factor(PY.patch.loc$V1, 
                                           levels = sort(unique(patch.loc$V1)))), 
                              exposure = "Pack Years")
PRS.patch.loc.v1 <- data.frame(table(factor(PRS.patch.loc$V1, 
                                            levels = sort(unique(patch.loc$V1)))), 
                               exposure = "PRS")
patch.loc.v1 <- rbind(overlap.patch.loc.v1, PY.patch.loc.v1, PRS.patch.loc.v1)
colnames(patch.loc.v1) <- c("X coordinates", "Number of significant patches",
                            "Exposure")
patch.loc.v1$Exposure <- factor(patch.loc.v1$Exposure,
                                levels = c("PRS", "Pack Years", "Both"))
ggplot(data = patch.loc.v1, aes(x = `X coordinates`, 
                                y = `Number of significant patches`, 
                                fill = Exposure)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("cyan3", "darkgoldenrod2", "red")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10))

# histogram - y
overlap.patch.loc.v2 <- data.frame(table(factor(overlap.patch.loc$V2, 
                                                levels = sort(unique(patch.loc$V2)))), 
                                   exposure = "Both")
PY.patch.loc.v2 <- data.frame(table(factor(PY.patch.loc$V2, 
                                           levels = sort(unique(patch.loc$V2)))), 
                              exposure = "Pack Years")
PRS.patch.loc.v2 <- data.frame(table(factor(PRS.patch.loc$V2, 
                                            levels = sort(unique(patch.loc$V2)))), 
                               exposure = "PRS")
patch.loc.v2 <- rbind(overlap.patch.loc.v2, PY.patch.loc.v2, PRS.patch.loc.v2)
colnames(patch.loc.v2) <- c("Y coordinates", "Number of significant patches",
                            "Exposure")
patch.loc.v2$Exposure <- factor(patch.loc.v2$Exposure,
                                levels = c("PRS", "Pack Years", "Both"))
ggplot(data = patch.loc.v2, aes(x = `Y coordinates`, 
                                y = `Number of significant patches`, 
                                fill = Exposure)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("cyan3", "darkgoldenrod2", "red")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10))

# histogram - z
overlap.patch.loc.v3 <- data.frame(table(factor(overlap.patch.loc$V3, 
                                                levels = sort(unique(patch.loc$V3)))), 
                                   exposure = "PRS & PY")
PY.patch.loc.v3 <- data.frame(table(factor(PY.patch.loc$V3, 
                                           levels = sort(unique(patch.loc$V3)))), 
                              exposure = "PY only")
PRS.patch.loc.v3 <- data.frame(table(factor(PRS.patch.loc$V3, 
                                            levels = sort(unique(patch.loc$V3)))), 
                               exposure = "PRS only")
patch.loc.v3 <- rbind(overlap.patch.loc.v3, PY.patch.loc.v3, PRS.patch.loc.v3)
colnames(patch.loc.v3) <- c("Z coordinates", "Number of significant patches",
                            "Exposure")
patch.loc.v3$Exposure <- factor(patch.loc.v3$Exposure,
                                levels = c("PRS only", "PY only", "PRS & PY"))
ggplot(data = patch.loc.v3, aes(x = `Z coordinates`, 
                                y = `Number of significant patches`, 
                                fill = Exposure)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("cyan3", "darkgoldenrod2", "red")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10))

# fisher's exact test
dat <- data.frame(PY_yes = c(9, 11),
                  PY_no = c(4, 557),
                  row.names = c("PRS_yes", "PRS_no"),
                  stringsAsFactors = FALSE
                  )
fisher.test(dat)






# MDS plot
all.selected.patch <- union(PRS.patch.idx, PY.patch.idx)
patch.loc <- patch.loc[all.selected.patch,]
mds.patch <- cmdscale(dist(patch.loc))
mds.patch.PRS <- mds.patch[row.names(mds.patch) %in% PRS.patch.idx,]
mds.patch.PY <- mds.patch[row.names(mds.patch) %in% PY.patch.idx,]
par(mfrow = c(3,1))
xlim <- range(mds.patch[,1]*1.1)
ylim <- range(mds.patch[,2]*1.1)
plot(mds.patch[,1], mds.patch[,2], xlab = "MDS Dimension 1",
     ylab = "MDS Dimension 2", xlim = xlim, ylim = ylim)
plot(mds.patch.PRS[,1], mds.patch.PRS[,2], xlab = "MDS Dimension 1",
     ylab = "MDS Dimension 2", col = "blue", 
     xlim = xlim, ylim = ylim)
plot(mds.patch.PY[,1], mds.patch.PY[,2], xlab = "MDS Dimension 1",
     ylab = "MDS Dimension 2", col = "red", 
     xlim = xlim, ylim = ylim)

# multi-plot
mds.patch.df <- data.frame(dim1 = mds.patch[,1],
                           dim2 = mds.patch[,2],
                           c1 = ifelse(row.names(mds.patch) %in% PRS.patch.idx,
                                       "1", NA),
                           c2 = ifelse(row.names(mds.patch) %in% PY.patch.idx,
                                       "1", NA))
mds.plt <- ggplot(mds.patch.df, aes(x = dim1, y = dim2, colour = c1)) +
  geom_point()
  geom_bar(stat="identity") +
  scale_x_discrete(drop=F) +
  scale_y_continuous(name="Pre-selected proportions", limits = c(0,1)) +
  scale_fill_manual(values=c("grey", "#31a354")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
bar.pvalue.df <- top50
bar.pvalue.plt <- ggplot(bar.pvalue.df, aes(x = mediatorID, y = -log10(pval.BY), fill = significance)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop=F, name="Top 50 mediators' ID") +
  scale_y_continuous(name="-log10 p-value") +
  scale_fill_manual(values=c("grey", "#31a354")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
# mix plots
cowplot::plot_grid(bar.prop.plt, bar.pvalue.plt,
                   ncol = 1, rel_heights = c(2, 2),
                   align = 'v', axis = 'lr')



# cluster 1 vs cluster 2
cluster <- ifelse(mds.patch[,1] > 0, "1", "2")
c1 <- as.numeric(names(cluster[cluster == "1"]))
c2 <- as.numeric(names(cluster[cluster == "2"]))

pval.all <- read.table(file = "PRS/final_result/pvalue_all.txt")
# keep cluster 1
cluster.pval <- pval.all[,1]
cluster.pval[cluster.pval < 2] <- 0   # only keep pval < 0.01
cluster.pval[c2] <- -cluster.pval[c2]
write.table(cluster.pval, file = "PRS/final_result/pvalue_cluster.txt",
            row.names = F, col.names = F)



######### end ##############

####### other method #########
# HIMA
HIMA.result <- hima(X = exposure,
                    Y = outcome,
                    M = mediators,
                    COV.XM = confounders,
                    COV.MY = confounders,
                    Y.family = "gaussian",
                    M.family = "gaussian",
                    parallel = T,
                    ncore = 5,
                    penalty = "MCP")
HIMA.result <- HIMA.result[order(HIMA.result$BH.FDR),]
head(HIMA.result)
save(HIMA.result, file = "result/HIMA_result.RData")

# aggregated p-value
patch.index <- as.numeric(sapply(row.names(HIMA.result), 
                                 function(x) strsplit(x, split = "M|-")[[1]][2]))
pval.all.HIMA <- rep(0, 581)
pval.all.HIMA[patch.index] <- -log10(HIMA.result$BH.FDR)
write.table(pval.all.HIMA, file = "result/3D_visualization/HIMA/pvalue_HIMA.txt",
            row.names = F, col.names = F)

# HP2016
HP2016.pval <- HP2016(M=mediators,
                      X=exposure,
                      Y=outcome,
                      C=confounders,
                      n.draw=100000,
                      adaptive=0.8)
print(HP2016.pval)































selected.mediators <- sapply(row.names(PS5.result$mediation.contri), function(x) which(x == paste0("M", c(1:p))))
PS5.med.contri.pval[rep.time,selected.mediators] <- PS5.result$mediation.contri$pval.Bonf


for(i in 1:rep.time){
  print(i)
  str.time <- Sys.time()
  result <- PECMAN(M = mediators,
                   X = exposure,
                   Y = outcome,
                   C = confounders,
                   dim.reduction = T,
                   seed = (i+1000),
                   n.draw = 100000)
  print(Sys.time() - str.time)
  print(result$pval)
  print(head(result$estimation))
  PECMAN.estimation[[i]] <- result$estimation
  if(any(result$pval == 0)){result$pval[result$pval == 0] <- 1e-16}
  PECMAN.pvalue[[i]] <- result$pval
}
save(PS5.estimation, file = "result/FEV1_estimation.RData")
save(PS5.pvalue, file = "result/FEV1_pvalue.RData")

# multisplit p-value aggregation
pvalue.matrix <- matrix(unlist(PECMAN.pvalue), ncol = 3, byrow = T)
overall.pvalue <- apply(pvalue.matrix, 2, function(x) multisplit.pval(x, min.gamma = 0.5))
# estimation part
gene.list.all <- unlist(lapply(PECMAN.estimation, function(x) row.names(x)))
gene.list.unique <- unique(gene.list.all)
# summarize results by gene names
gene.selected.proportion <- c()
gene.selected.proportion <- sapply(gene.list.unique, function(x) sum(gene.list.all == x)/rep.time)
gene.selected.proportion <- sort(gene.selected.proportion, decreasing = T)
print(head(gene.selected.proportion))
all.estimation <- do.call(rbind, PECMAN.estimation)
mediation.prop <- NULL
marginal.pvalue.BH <- marginal.pvalue.Bonf <- rep(NA, length(gene.list.unique))
for(i in 1:length(gene.list.unique)){ # top 50 genes ranked by selected proportion
  gene.idx <- which(gene.list.all == names(gene.selected.proportion)[i])
  tmp.gene.res <- all.estimation[gene.idx,]  # result for single gene
  # p-value
  tmp.gene.res$p.value.BH[tmp.gene.res$p.value.BH == 0] <- 1/100000
  tmp.gene.res$p.value.Bonf[tmp.gene.res$p.value.Bonf == 0] <- 1/100000
  tmp.pvalue.BH <- c(tmp.gene.res$p.value.BH, rep(1, rep.time - nrow(tmp.gene.res)))
  tmp.pvalue.Bonf <- c(tmp.gene.res$p.value.Bonf, rep(1, rep.time - nrow(tmp.gene.res)))
  # keep 1/2 informative result to aggregate p-value
  marginal.pvalue.BH[i] <- multisplit.pval(tmp.pvalue.BH, min.gamma = nrow(tmp.gene.res)/rep.time/2)
  marginal.pvalue.Bonf[i] <- multisplit.pval(tmp.pvalue.Bonf, min.gamma = nrow(tmp.gene.res)/rep.time/2)
  #marginal.pvalue.BH[i] <- harmonic.mean(tmp.pvalue.BH)
  #marginal.pvalue.Bonf[i] <- harmonic.mean(tmp.pvalue.Bonf)
  # mediation effect (if p-value > 0.05, then mediation = 0)
  #tmp.mediation.prop <- tmp.gene.res$mediation.proportion*(tmp.gene.res$p.value.BH<0.05)
  
  mediation.prop <- rbind(mediation.prop, data.frame(names(gene.selected.proportion)[i],
                                                     tmp.gene.res$mediation.proportion,
                                                     marginal.pvalue.BH[i]<0.05))
}


# check 451~453 patch
significant.patch <- matrix(0, nrow = rep.time, ncol = 3)
colnames(significant.patch) <- c("M451", "M452","M453")
for(i in 1:rep.time){
  patch.level.name <- sapply(row.names(PECMAN.estimation[[i]]), function(x) strsplit(x, split = "M|-")[[1]][2])
  if("451" %in% patch.level.name) significant.patch[i,1] <- 1
  if("452" %in% patch.level.name) significant.patch[i,2] <- 1
  if("453" %in% patch.level.name) significant.patch[i,3] <- 1
}
table(rowSums(significant.patch))
###############end################
##########visualization via multiple figures#################
# keep top 50 mediators with higher p-value
top50.idx <- order(marginal.pvalue.BH)[1:50]
top50.name <- names(gene.selected.proportion)[top50.idx]
marginal.pvalue.BH <- marginal.pvalue.BH[top50.idx]
gene.selected.proportion <- gene.selected.proportion[top50.idx]
# boxplot for estimation
box.df <- mediation.prop
colnames(box.df) <- c('mediatorID',"Mediation proportion", "Significant")
box.df$mediatorID <- factor(box.df$mediatorID, levels = top50.name)
box.df <- box.df[!is.na(box.df$mediatorID),]
#box.df$geneID <- factor(box.df$geneID, levels = names(gene.selected.proportion)[1:50])
box.plt <- ggplot(box.df, aes(x=mediatorID, y=`Mediation proportion`, fill = Significant)) +
  geom_boxplot() +
  scale_fill_manual(values=c("grey", "#31a354")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
bar.prop.df <- box.df
bar.prop.plt <- ggplot(bar.prop.df, aes(x = mediatorID, fill = Significant)) +
  geom_bar(aes(y=..count../rep.time)) +
  scale_x_discrete(drop=F) +
  scale_y_continuous(name="Pre-selected proportions") +
  scale_fill_manual(values=c("grey", "#31a354")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
bar.pvalue.df <- data.frame(mediatorID = names(gene.selected.proportion),
                            pvalue = -log10(marginal.pvalue.BH),
                            Significant = marginal.pvalue.BH<0.05)
bar.pvalue.df$mediatorID <- factor(bar.pvalue.df$mediatorID, levels = top50.name)
bar.pvalue.plt <- ggplot(bar.pvalue.df, aes(x = mediatorID, y = pvalue, fill = Significant)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop=F, name="Top 50 mediators' ID") +
  scale_y_continuous(name="-log10 p-value") +
  scale_fill_manual(values=c("grey", "#31a354")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
# mix plots
png(file = "result/FEV1_result.png", height = 2200, width = 2200,res = 300)
cowplot::plot_grid(box.plt, bar.prop.plt, bar.pvalue.plt,
                   ncol = 1, rel_heights = c(3, 2, 2),
                   align = 'v', axis = 'lr')
dev.off()















# visualization
p.alpha <- apply(mediators, 2,
                 function(x) summary(lm(x ~ exposure + confounders))$coefficients[2,4])
p.beta  <- apply(mediators, 2,
                 function(x) summary(lm(outcome ~ x + exposure + confounders))$coefficients[2,4])
pvalue.df <- data.frame(geneID = colnames(mediators),
                        alpha = -log10(p.alpha),
                        beta = -log10(p.beta))

scatter.plt <- ggplot(pvalue.df, aes(x=alpha, y=beta)) +
  geom_point() +
  geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  xlab("log p-value of X-M association") +
  ylab("log p-value of M-Y association") +
  scale_x_continuous(limits = c(0, 10))
scatter.plt

# HIMA
HIMA.result <- hima(X = exposure,
                    Y = outcome,
                    M = mediators,
                    COV.XM = confounders,
                    COV.MY = confounders,
                    Y.family = "gaussian",
                    M.family = "gaussian",
                    parallel = T,
                    ncore = 5,
                    penalty = "MCP")
HIMA.result <- HIMA.result[order(HIMA.result$Bonferroni.p),]
print(head(HIMA.result))

# HP2016
HP2016.pval <- HP2016(M=mediators,
                      X=exposure,
                      Y=outcome,
                      C=confounders,
                      n.draw=100000,
                      adaptive=0.8)
print(HP2016.pval)

########################################
### 2022/12/20
########################################
### analyze without multiple splitting
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/")
for(snps.idx in 3:13){
  #imaging.all <- read.csv("context_matter_patch.csv", header = F)
  #load("context_matter_patch_new.RData")
  load("context_matter_patch_pc1.RData")
  image.id.df <- read.csv("context_matter_id.csv", header = F)
  image.id <- apply(image.id.df, 1, function(x) strsplit(x, "_")[[1]][1])
  clinic.df <- read.csv("clinical_surv_df.txt", sep = " ")
  confounding.df <- read.csv("covariate_df.txt", sep = " ")
  snp.df <- read.csv("risk_loci_nature_com.txt", sep = "\t")
  snp.df <- snp.df[,c(1,2,7,9,8,10,15,12,14,16,11,13,17)]   # re-order
  #dir.create(paste0("pc1_result/", colnames(snp.df)[snps.idx]))
  # match id
  match.idx <- intersect(intersect(image.id, clinic.df$sid), snp.df$IID)
  clinic.df <- clinic.df[sapply(match.idx, function(x) which(clinic.df$sid == x)),]
  snp.df <- snp.df[sapply(match.idx, function(x) which(snp.df$IID == x)),]
  confounding.df <- confounding.df[sapply(match.idx, function(x) which(confounding.df$IID == x)),]
  imaging.all <- as.matrix(imaging.all[sapply(match.idx, function(x) which(image.id == x)),])   # mediators
  FEV1pp <- clinic.df$FEV1pp_utah                                                               # outcome
  conf <- as.matrix(confounding.df[,3:18])                                                      # confounders
  snp <- as.matrix(snp.df[,snps.idx])                                                                 # exposure (7~15)
  
  result <- NULL
  result <- PECMAN(M = imaging.all,
                   X = snp,
                   Y = FEV1pp,
                   C = conf,
                   dim.reduction = F,
                   seed = 1001,
                   n.draw = 10000)
  print(result$pval)
  print(head(result$estimation))
  save(result, file = paste0("pc1_result/", colnames(snp.df)[snps.idx], ".RData"))
  
  ### visualization
  p.alpha <- apply(imaging.all, 2,
                   function(x) summary(lm(x ~ snp + conf))$coefficients[2,4])
  p.beta  <- apply(imaging.all, 2,
                   function(x) summary(lm(FEV1pp ~ x + snp + conf))$coefficients[2,4])
  pvalue.df <- data.frame(geneID = colnames(imaging.all),
                          alpha = -log10(p.alpha),
                          beta = -log10(p.beta))
  
  scatter.plt <- ggplot(pvalue.df, aes(x=alpha, y=beta)) +
    geom_point() +
    geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "red") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    xlab("log p-value of X-M association") +
    ylab("log p-value of M-Y association") +
    scale_x_continuous(limits = c(0, 10))
  scatter.plt
}





################################################
################# patch level ##################
################################################
setwd("~/OneDrive/Pitts/COPD/mediation_analysis/real_application/COPD/")
for(snps.idx in 3:13){
  #imaging.all <- read.csv("context_matter_patch.csv", header = F)
  #load("context_matter_patch_new.RData")
  load("context_matter_patch_pc1.RData")
  image.id.df <- read.csv("context_matter_id.csv", header = F)
  image.id <- apply(image.id.df, 1, function(x) strsplit(x, "_")[[1]][1])
  clinic.df <- read.csv("clinical_surv_df.txt", sep = " ")
  confounding.df <- read.csv("covariate_df.txt", sep = " ")
  snp.df <- read.csv("risk_loci_nature_com.txt", sep = "\t")
  snp.df <- snp.df[,c(1,2,7,9,8,10,15,12,14,16,11,13,17)]   # re-order
  #dir.create(paste0("pc1_result/", colnames(snp.df)[snps.idx]))
  # match id
  match.idx <- intersect(intersect(image.id, clinic.df$sid), snp.df$IID)
  clinic.df <- clinic.df[sapply(match.idx, function(x) which(clinic.df$sid == x)),]
  snp.df <- snp.df[sapply(match.idx, function(x) which(snp.df$IID == x)),]
  confounding.df <- confounding.df[sapply(match.idx, function(x) which(confounding.df$IID == x)),]
  imaging.all <- as.matrix(imaging.all[sapply(match.idx, function(x) which(image.id == x)),])   # mediators
  FEV1pp <- clinic.df$FEV1pp_utah                                                               # outcome
  conf <- as.matrix(confounding.df[,3:18])                                                      # confounders
  snp <- as.matrix(snp.df[,snps.idx])                                                                 # exposure (7~15)
  
  result <- NULL
  count <- 1
  for(seed in 1001:1100){
    result[[count]] <- PECMAN(M = imaging.all,
                              X = snp,
                              Y = FEV1pp,
                              C = conf,
                              dim.reduction = T,
                              seed = seed,
                              n.draw = 10000,
                              pre.selection = "MCP")
    print(result[[count]]$pval)
    print(head(result[[count]]$estimation))
    count <- count + 1
  }
  PECMAN.pval <- sapply(result, function(x) x$pval[1])
  pval.df <- data.frame(pval = PECMAN.pval)
  comb.pval <- 1- pcauchy(1/length(PECMAN.pval)*sum(tan(pi*(0.5-PECMAN.pval))))
  comb.pval
  save(result, file = paste0("pc1_result/PECMAN_result_", colnames(snp.df)[snps.idx], ".RData"))
  # patch weight
  patch_weight <- matrix(0, nrow = 581, ncol = 100)
  for(replication.idx in 1:100){
    feature_id <- sapply(rownames(result[[replication.idx]]$estimation),
                         function(x) as.numeric(strsplit(x, "V")[[1]][2]))
    patch_weight[feature_id,replication.idx] <- (result[[replication.idx]]$estimation$mediation.proportion)
  }
  median_effect <- abs(sapply(1:581, function(x) median(patch_weight[x,])))
  write.table(median_effect, file = paste0("pc1_result/patch_weight_", colnames(snp.df)[snps.idx], ".txt"), row.names = F, col.names=F)
}

# scree plot of pre-selected proportion
nonzero <- patch_weight != 0
sort_nonzero <- sort(rowSums(nonzero), decreasing = T)
plot(sort_nonzero,
     col = ifelse(sort_nonzero > 60, 2, 1),
     ylab = "pre-selected proportion")

# boxplot for estimation
library(ggplot2)
patch_weight_selected <- patch_weight[rowSums(nonzero) > 60,]
weight_df <- stack(as.data.frame(t(patch_weight_selected)))
ggplot(weight_df, aes(x=ind, y=values)) +
  geom_boxplot() +
  theme_classic()


for(snps.idx in 3:11){
  #dir.create(paste0("pc1_result/", colnames(snp.df)[snps.idx]))
  #imaging.all <- read.csv("context_matter_patch.csv", header = F)
  #load("context_matter_patch_new.RData")
  load("context_matter_patch_pc1.RData")
  image.id.df <- read.csv("context_matter_id.csv", header = F)
  image.id <- apply(image.id.df, 1, function(x) strsplit(x, "_")[[1]][1])
  clinic.df <- read.csv("clinical_surv_df.txt", sep = " ")
  confounding.df <- read.csv("covariate_df.txt", sep = " ")
  snp.df <- read.csv("risk_loci_nature_com.txt", sep = "\t")
  snp.df <- snp.df[,c(1,2,7,8,13,10,12,14,9,11,15)]   # re-order
  # match id
  match.idx <- intersect(intersect(image.id, clinic.df$sid), snp.df$IID)
  clinic.df <- clinic.df[sapply(match.idx, function(x) which(clinic.df$sid == x)),]
  snp.df <- snp.df[sapply(match.idx, function(x) which(snp.df$IID == x)),]
  confounding.df <- confounding.df[sapply(match.idx, function(x) which(confounding.df$IID == x)),]
  imaging.all <- as.matrix(imaging.all[sapply(match.idx, function(x) which(image.id == x)),])   # mediators
  FEV1pp <- clinic.df$FEV1pp_utah                                                               # outcome
  conf <- as.matrix(confounding.df[,3:18])                                                      # confounders
  snp <- as.matrix(snp.df[,snps.idx])                                                                 # exposure (7~15)
  HIMA.result <- hima(X = snp,
                      Y = FEV1pp,
                      M = imaging.all,
                      COV.XM = conf,
                      COV.MY = conf,
                      Y.family = "gaussian",
                      M.family = "gaussian",
                      parallel = T,
                      ncore = 5,
                      penalty = "MCP")
  HIMA.result <- HIMA.result[order(HIMA.result$Bonferroni.p),]
  print(head(HIMA.result))
  
  HP2016.pval <- HP2016(M=imaging.all,
                        X=snp,
                        Y=FEV1pp,
                        C=conf,
                        n.draw=100000,
                        adaptive=0.8)
  print(HP2016.pval)
}

### visualization
snp.df <- read.csv("risk_loci_nature_com.txt", sep = "\t")
snp.df <- snp.df[,c(1,2,7,9,8,10,15,12,14,16,11,13,17)]   # re-order
median_effect <- matrix(nrow = 581, ncol = 9)
for(snps.idx in 3:11){
  patch_weight <- read.table(paste0("pc1_result/", colnames(snp.df)[snps.idx], "/patch_weight_", colnames(snp.df)[snps.idx], ".txt"))
  median_effect[,(snps.idx-2)] <- patch_weight[,1]
}







###################################################
## calculate average of mediation proportion
sig.feature <- NULL
for(i in 1:100){
  feature_id <- sapply(rownames(result[[i]]$estimation),
                       function(x) as.numeric(strsplit(x, "V")[[1]][2]))
  result[[i]]$estimation$feature.id <- feature_id
  sig.feature <- rbind(sig.feature,
                       result[[i]]$estimation[result[[i]]$estimation$p.value < 0.05,])
}
med.prop.avg <- sapply(unique(sig.feature$feature.id),
                       function(x) mean(sig.feature$mediation.proportion[sig.feature$feature.id == x]))
significant.count <- sapply(unique(sig.feature$feature.id),
                            function(x) sum(sig.feature$feature.id == x))
PECMAN.table <- data.frame(feature.id = unique(sig.feature$feature.id),
                           patch = ceiling(unique(sig.feature$feature.id)/128),
                           `avg mediation proportion` = med.prop.avg,
                           `significant.time` = significant.count)
PECMAN.table <- PECMAN.table[order(PECMAN.table$significant.time, decreasing = T),]
write.csv(PECMAN.table, file = "PECMAN_result_rs12914385.csv", row.names = F)
## visualize by significant time over 100 replications
patch.weight <- rep(0, 581)
patch.weight[PECMAN.table$patch] <- PECMAN.table$significant.time
write.table(patch.weight, file = "patch_weight_rs12914385.txt", row.names = F, col.names = F)

