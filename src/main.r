rm(list = ls())

malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));gc()

graphics.off()

library(data.table)
library(data.table)
library(dplyr)
library(tidyr)
library(dtw)
library(dtwclust)
library(caret)
library(ggdendro)
library(ggalluvial)
library(ggplot2)
library(mice)
library(parallel)
library(tidyr)
library(TTR)
library(randomForest)
library(readstata13)
library(reshape2)
library(tidyr)
library(zoo)
library(scico)
library(imputeTS)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(purrr)
library(dtw)
library(dplyr)
library(matrixStats)
library(cluster)
library(survival) 
library(survminer) 

do.main.alignments = F

data.path0 = "dp0"
data.path1 = "dp1"


source("preproc_alignment_functions.r")

hcc_data = read.dta13(paste(data.path0,"/hcc_data.dta", sep = ""))
ids.ever.hcc = unique(hcc_data$projectid)

if(do.main.alignments == T) source("heading_correct_labdata_allpats_units.r")

if(do.main.alignments == T){
    lab.ts.after.tdf.interpolated.pp.full = lab.ts.after.tdf.interpolated
    ids.ever.hcc.passed = intersect(unique(ids.ever.hcc), names(lab.ts.after.tdf.interpolated))
    set.seed = 1359
    calib_prop <- 0.7  # 70% calibration, 30% test
    ## Shuffle the indices
    shuffled_indices <- sample(ids.ever.hcc.passed)
    ## Compute split point
    N = length(unique(ids.ever.hcc.passed))
    n_calib <- floor(N * calib_prop)
    ## Split indices
    calib_ids.ever.hcc <- shuffled_indices[1:n_calib]
    test_ids.ever.hcc  <- shuffled_indices[(n_calib + 1):N]
    message("in main arrancando con el hcc set")
    
    write.table(calib_ids.ever.hcc,"calib_ids_ever_hcc")
    write.table(test_ids.ever.hcc,"test_ids_ever_hcc")
    
    lab.ts.after.tdf.interpolated.pp = lab.ts.after.tdf.interpolated[!is.element(names(lab.ts.after.tdf.interpolated), calib_ids.ever.hcc)]
}


calib_ids.ever.hcc = read.table("calib_ids_ever_hcc")[,1]
test_ids.ever.hcc = read.table("test_ids_ever_hcc")[,1]
    
file.name.consensus = "consensus_reference_df_allpat"
file.name.alignment = "alignment_al_allpat"
file.name.figure.consensus = "reference_all.png"
file.name.imputation_outcome_example = "imputation_outcome_example.png"
file.name.original_to_aligned = "original_to_aligned.png"
if(do.main.alignments == T) source("generic_alignment_correct_labdata.r")

message("Started alignment of the all set")
#now the ones with hcc diagnosis
if(do.main.alignments == T) { ids.hcc.to.align = intersect(calib_ids.ever.hcc, names(lab.ts.after.tdf.interpolated)) }
if(do.main.alignments == T) { lab.ts.after.tdf.interpolated.pp = lab.ts.after.tdf.interpolated[is.element(names(lab.ts.after.tdf.interpolated), ids.hcc.to.align)]}
file.name.consensus = "consensus_reference_df_hcc"
file.name.alignment = "alignment_al_hcc"
file.name.figure.consensus = "reference_hcc.png"
file.name.imputation_outcome_example = "imputation_outcome_example_hcc.png"
file.name.original_to_aligned = "original_to_aligned_hcc.png"
if(do.main.alignments == T) source("generic_alignment_correct_labdata.r")

consensus.reference.df.hbv.hcc = read.table(file.name.consensus)
consensus.reference.df.hbv.hcc$population = "hbv_hcc"
alignment.al.hbv.hcc.with.ids = read.table(file.name.alignment)
alignment.al.hbv.hcc.with.ids$population = "hbv_hcc"

consensus.reference.df.allpat = read.table("consensus_reference_df_allpat")
consensus.reference.df.allpat$population = "all"

message("1. Align reference time series")
consensus.references.df. =  rbind(consensus.reference.df.hbv.hcc, consensus.reference.df.allpat)
consensus.references.df.$id = 0

warning("I am alignning the actual time series [after inputation and interpolation]")
alignment.al.allpat = read.table("lab_ts_after_tdf_interpolated_ids")

alignment.al.allpat$population = "all"

message("read the unaligned all set for alingment to HCC set")
warning("ojo, make sure the names of the variables are in the same order")

alignments.df =  rbind(alignment.al.hbv.hcc.with.ids, alignment.al.allpat)

x11(title = "cons references");
pp = ggplot(consensus.references.df., aes(x= time, y= value, col= population, lty = population)) +
    facet_wrap(variable~., scale = "free", ncol = 3) +
    geom_line() +
    theme_classic() + scale_color_manual(values = c("all" = "green", "hbv_hcc" = "red"))
print(pp)

png("references.png", width = 5, height = 12.5, res = 300, units = "in");print(pp);dev.off()

alignments.df.long <- alignments.df %>%
    pivot_longer(cols = names(alignments.df)[-grep("time|id|population",names(alignments.df))], names_to = "variable", values_to = "value")

x11(title = "post alignment  - todo")
pp = ggplot(alignments.df.long, aes(x = time, y = value, grouping = as.factor(id), col = as.factor(population))) +
                             geom_line(size = 0.5, alpha = 0.25) + 
    theme_minimal() +   facet_wrap(~variable, ncol = 3, scale = "free")   + scale_y_log10() +
    geom_line(data = consensus.references.df., size = 0.5, col = "black", aes(lty =  as.factor(population)))

print(pp)

png("alignments_independent.png", width = 7, height = 12.5, res = 300, units = "in");print(pp);dev.off()

consensus.reference.df.hbv.hcc.short = dcast(consensus.reference.df.hbv.hcc, time + population ~ variable)
consensus.reference.df.allpat.short = dcast(consensus.reference.df.allpat, time + population ~ variable)

consensus.reference.df.hbv.hcc.short$apri = with(consensus.reference.df.hbv.hcc.short,(AST/50)*100/pla)
consensus.reference.df.allpat.short$apri = with(consensus.reference.df.allpat.short,(AST/50)*100/pla)
consensus.reference.df.hbv.hcc.short$cd4.cd8.ratio = with(consensus.reference.df.hbv.hcc.short, 100*cd4/cd8)
consensus.reference.df.allpat.short$cd4.cd8.ratio = with(consensus.reference.df.allpat.short, 100*cd4/cd8)

alignment.al.hbv.hcc.with.ids = alignment.al.hbv.hcc.with.ids[,c("id", names(consensus.reference.df.hbv.hcc.short)[-grep("apri|cd4.cd8.ratio", names(consensus.reference.df.hbv.hcc.short))]  )]

hcc.hbv.alignment.list = by(alignment.al.hbv.hcc.with.ids, alignment.al.hbv.hcc.with.ids$id, function(x) { x[, -grep("id|population|time", names(x))]} )

##all to hcc consensus reference
 alignment <- dtw(consensus.reference.df.hbv.hcc.short[,-grep(c("time|population|apri|cd4.cd8.ratio"), names(consensus.reference.df.hbv.hcc.short))],
                  consensus.reference.df.allpat.short[,-grep(c("time|population|apri|cd4.cd8.ratio"), names(consensus.reference.df.allpat.short))],
                  keep = T, step.pattern = rabinerJuangStepPattern(6, "c"),
                  dist.method = "Euclidean"
                  )

 aligned <- as.data.frame(matrix(NA, nrow = nrow(consensus.reference.df.hbv.hcc.short) , ncol = ncol(consensus.reference.df.hbv.hcc.short)))
 aligned[alignment$index1, ] <- consensus.reference.df.allpat.short[alignment$index2, ]  # Map alignment
 names(aligned) =  names(consensus.reference.df.allpat.short)## 
 aligned.references = aligned

 aligned.references.long <- aligned.references %>%
     pivot_longer(cols = names(aligned.references)[-grep("time|population",names(aligned.references))], names_to = "variable", values_to = "value")

aligned.references.long$population = "all, aligned to hcc"

x11(title = "cons references +  aligned allpat");
pp= ggplot(consensus.references.df., aes(x= time, y= value, col= population)) +
    facet_wrap(variable~., scale = "free", ncol = 3) +
    geom_line() + 
    scale_y_log10() +
    geom_line(data = aligned.references.long, lty = 2) +
    scale_color_manual(values = c("all" = "green", "hbv_hcc" = "red", "all, aligned to hcc" = "blue")) + theme_classic()
print(pp)

png("references_all_to_hcchbv.png", width = 5, height = 12.5, res = 300, units = "in");print(pp);dev.off()

alignment.al.allpat. = alignment.al.allpat[,c("id",names(consensus.reference.df.hbv.hcc.short)[-grep("apri|cd4.cd8.ratio", names(consensus.reference.df.hbv.hcc.short))])]
alignment.al.allpat.list =  by(alignment.al.allpat., alignment.al.allpat.$id, function(x) x[,-grep("id|population",colnames(alignment.al.allpat.))] )


message("2. About to align hcc - no hcc")
alignment.all.to.hcc.list =  alignment_to_reference(series_list = alignment.al.allpat.list,
                                reference.ts = consensus.reference.df.hbv.hcc.short[,-grep("apri|population|cd4.cd8.ratio", names(consensus.reference.df.hbv.hcc.short))] )

names(alignment.all.to.hcc.list) = unlist(lapply(alignment.all.to.hcc.list, function(x) x$id[1]))

alignment.all.to.hcc.df = do.call(rbind, alignment.all.to.hcc.list)

alignment.all.to.hcc.df.long = alignment.all.to.hcc.df %>%
    pivot_longer(cols = names(consensus.reference.df.hbv.hcc.short)[-grep("time|population|apri|cd4.cd8.ratio",names(consensus.reference.df.hbv.hcc.short))], names_to = "variable", values_to = "value")

##now compute a post-concensus:

alignment.all.to.hcc.df.consensus =  as.data.frame(do.call(rbind, by(alignment.all.to.hcc.df, alignment.all.to.hcc.df$time, function(x) {
    time = x[1,"time"]
    x =  x[, -grep("time|id", colnames(x))]
    out = c(time = time, apply(x,2,median))
    out
})))

alignment.all.to.hcc.df.consensus.long = alignment.all.to.hcc.df.consensus %>%
    pivot_longer(cols = names(alignment.all.to.hcc.df.consensus)[-grep("time|population",names(alignment.all.to.hcc.df.consensus))], names_to = "variable", values_to = "value")
alignment.all.to.hcc.df.consensus.long$id = 1
alignment.all.to.hcc.df.consensus.long$population = "all -al to hcc hbv" 

alignment.all.to.hcc.df.consensus.long[which(alignment.all.to.hcc.df.consensus.long$variable =="rna"),"value"] = 10^alignment.all.to.hcc.df.consensus.long[which(alignment.all.to.hcc.df.consensus.long$variable == "rna"),"value"]

consensus.references.df.revert = consensus.references.df.

consensus.references.df.revert[ which(consensus.references.df.revert$variable == "rna"),"value"] = 10^consensus.references.df.revert[ which(consensus.references.df.revert$variable == "rna"),"value"]

alignment.all.to.hcc.df.long.revert = alignment.all.to.hcc.df.long
alignment.all.to.hcc.df.long.revert[alignment.all.to.hcc.df.long.revert$variable == "rna","value"] = 10^alignment.all.to.hcc.df.long.revert[alignment.all.to.hcc.df.long.revert$variable == "rna","value"]

rm(alignments.combined.init);rm(alignment.all.to.hcc.df.long); malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));gc()

x11(title = "alignment all to hcc  con apri")
pp = ggplot(alignment.all.to.hcc.df.long.revert, aes(x = time, y = value, grouping = as.factor(id))) +
    theme_minimal() +   facet_wrap(~variable, ncol = 6, scale = "free")   + 
    geom_line(data = consensus.references.df.revert, size = 0.5,  aes(lty =  as.factor(population), col = as.factor(population) )) +
    scale_color_manual(values = c("all" = "green", "hbv_hcc" = "red", "all -al to hcc hbv" =  "blue")) +
    geom_line(data = alignment.all.to.hcc.df.consensus.long, size = 0.5, col = "blue", aes(lty =  as.factor(population))) + scale_y_log10() 
                
print(pp)    

png("all_to_consensus_hcchbv.png", width = 10, height = 10, res = 300, units = "in");print(pp);dev.off()

message("Are there regions in time that contribute more to the costs of the alignment")

result = compute_bootstrap_confidence_dtw(aligned_list = hcc.hbv.alignment.list,  consensus = consensus.reference.df.allpat.short[,-grep("time|id|apri|population|cd4.cd8.ratio",names(consensus.reference.df.allpat.short))], n_boot = 1000)

x11(title =  "bootstrap confidence dtw");
par(mfrow = c(2,1))
plot(consensus.reference.df.allpat.short$time, result$observed_contributions,
     cex = 0.5, xlab = "Time", ylab = "Contribution to costs", cex.lab = 2);
lines(consensus.reference.df.allpat.short$time,result$ci_lower, lty = 2)
lines(consensus.reference.df.allpat.short$time,result$ci_upper, lty = 2)
plot(consensus.reference.df.allpat.short$time, result$p_values, cex = 0.5, xlab = "Time", ylab = "p-value")
lines(consensus.reference.df.allpat.short$time, result$p_values, cex = 0.2)

png("bootstrap_confidence_dtw.png", width = 12.5, height = 7.25, res = 300, units = "in")
par(mfrow = c(1,2))
plot(consensus.reference.df.allpat.short$time, result$observed_contributions, cex = 0.5, xlab = "Time", ylab = "Contribution to costs",cex.lab = 1.5);
lines(consensus.reference.df.allpat.short$time,result$ci_lower, lty = 2)
lines(consensus.reference.df.allpat.short$time,result$ci_upper, lty = 2)
plot(consensus.reference.df.allpat.short$time, result$p_values, cex = 0.5, xlab = "Time", ylab = "p-value",cex.lab = 1.5)
lines(consensus.reference.df.allpat.short$time, result$p_values, cex = 0.2)
dev.off()

message("3. BRUT comparison between the negative and positive alignments:")

graphics.off()

hcc.hbv.alignment.list.withcolas   = hcc.hbv.alignment.list

##
no.cancer.ids = setdiff(names(alignment.all.to.hcc.list), ids.ever.hcc)
calib_prop <- 0.7  # 80% calibration, 20% test
## Shuffle the indices
shuffled_indices <- sample(no.cancer.ids)
## Compute split point
N = length(unique(no.cancer.ids))
n_calib <- floor(N * calib_prop)
## Split indices
calib_ids.nohcc <- shuffled_indices[1:n_calib]
test_ids.nohcc  <- shuffled_indices[(n_calib + 1):N]
write.table(calib_ids.nohcc,"calib_ids_nohcc")
write.table(test_ids.nohcc,"test_ids_nohcc")
##

alignment.all.to.hcc.list.ALL = alignment.all.to.hcc.list
alignment.all.to.hcc.list =  alignment.all.to.hcc.list.ALL[which(is.element(names(alignment.all.to.hcc.list), calib_ids.nohcc))]
alignment.all.to.hcc.list = lapply(alignment.all.to.hcc.list, function(x) x[,-grep("time|id", colnames(x))] )
alignment.all.to.hcc.list.withcolas = alignment.all.to.hcc.list

warning("removing time series flat ends before comparing alignments")
hcc.hbv.alignment.list   = lapply(hcc.hbv.alignment.list.withcolas, remove_flat_ends)
alignment.all.to.hcc.list = lapply(alignment.all.to.hcc.list.withcolas, remove_flat_ends)

hcc.hbv.alignment.list =  lapply(hcc.hbv.alignment.list, function(m) {
    m$apri           = with(m,(AST/50)*100/pla)
    m$cd4.cd8.ratio = with(m,100*cd4/cd8)
    m
})

alignment.all.to.hcc.list  =  lapply(alignment.all.to.hcc.list, function(m) {
    m$apri           = with(m,(AST/50)*100/pla)
    m$cd4.cd8.ratio = with(m,100*cd4/cd8)
    m
})

message("computing the differences in level")
aa = compute_median_variance(group1=(hcc.hbv.alignment.list), group2 = (alignment.all.to.hcc.list))
colnames(aa$median_diff) =  colnames(hcc.hbv.alignment.list[[1]])
colnames(aa$variance) =  colnames(hcc.hbv.alignment.list[[1]])
colnames(aa$sign_diff) =  colnames(hcc.hbv.alignment.list[[1]])

hcc.hbv.alignment.list.deriv   =  lapply(hcc.hbv.alignment.list, function(m){ m[is.na(m)]=0 ;cumsum(m[2:nrow(m),] - m[1:(nrow(m)-1),])} )
alignment.al.allpat.list.deriv =  lapply(alignment.all.to.hcc.list, function(m){m[is.na(m)]=0 ; cumsum(m[2:nrow(m),] - m[1:(nrow(m)-1),])} )

message("computing the differences in derivative")
dd =  compute_median_variance(group1=(hcc.hbv.alignment.list.deriv), group2 = (alignment.al.allpat.list.deriv), relative = FALSE)
dd.init = dd
colnames(dd$sign_diff) =  colnames(hcc.hbv.alignment.list[[1]])
colnames(aa$pair.comparison.availability) =  colnames(hcc.hbv.alignment.list[[1]])
aa$pair.comparison.availability$time = consensus.reference.df.hbv.hcc.short$time

median_diff_hcc_all = aa$median_diff
median_diff_hcc_all.init = median_diff_hcc_all
median_diff_hcc_all$time = consensus.reference.df.hbv.hcc.short$time

df_long <- melt(median_diff_hcc_all, id.vars = "time")  ##donde difieren más
colnames(df_long) <- c("time", "variable", "value")

# Create heatmap
x11(title = "dtw discrepancy matrix map - meandiff")
pp = ggplot(df_long, aes(x = time, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  ## Rotate x-axis labels
    labs(title = "", fill = "value")
print(pp)

sign_diff_hcc_all = aa$sign_diff  
sign_diff_hcc_all$time = consensus.reference.df.hbv.hcc.short$time

sign_deriv_diff_hcc_all = as.data.frame(dd$sign_diff)
sign_deriv_diff_hcc_all$time = consensus.reference.df.hbv.hcc.short$time[2:length(consensus.reference.df.hbv.hcc.short$time)]

df_long <- melt(sign_diff_hcc_all, id.vars = "time")
colnames(df_long) <- c("time", "variable", "value")
# Create heatmap
x11(title = "dtw discrepancy matrix map - sign diff")
pp = ggplot(df_long, aes(x = time, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_minimal() +
    theme(
        axis.text.x = element_text(hjust = 1,size=17),
        axis.text.y = element_text(hjust = 1,size=17),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)
          ) + 
    labs(title = "", fill = "value")
print(pp)
png("dtw_signdiff_all.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()

df_long <- melt(sign_deriv_diff_hcc_all, id.vars = "time")
colnames(df_long) <- c("time", "variable", "value")
x11(title = "dtw discrepancy matrix map - deriv sign diff")
pp = ggplot(df_long, aes(x = time, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_minimal() +
    theme(
        axis.text.x = element_text(hjust = 1,size=17),
        axis.text.y = element_text(hjust = 1,size=17),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)
          ) +
    labs(title = "", fill = "value")
print(pp)
png("dtw_signdiff_deriv_all.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()

df_long <- melt(aa$pair.comparison.availability, id.vars = "time")
colnames(df_long) <- c("time", "variable", "value")

x11(title = "dtw_pair_comparison_availability")
pp = ggplot(df_long, aes(x = time, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_minimal() +
    theme(
axis.text.x = element_text(hjust = 1,size=17),
axis.text.y = element_text(hjust = 1,size=17),
axis.title.x = element_text(size=20),
axis.title.y = element_text(size=20)
          ) + 
    labs(title = "", fill = "value")
print(pp)
png("dtw_pair_comparison_availability.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()

message("computing the levels for the comparison matrices")
##selecting the region
sign_diff_hcc_all.level       = compute_dtw_discrepancy_matrix.level(sign_diff_hcc_all, similar= FALSE, threshold = "75%")
sign_diff_hcc_deriv_all.level = compute_dtw_discrepancy_matrix.level(sign_deriv_diff_hcc_all, similar= FALSE, threshold = "75%")
pair.comparability.level =  compute_dtw_discrepancy_matrix.level(aa$pair.comparison.availability, similar= FALSE, threshold = "25%")

cancer.region.1 = as.matrix(1-is.na(sign_diff_hcc_all.level))
cancer.region.2 = as.matrix(rbind(rep(0,ncol(sign_diff_hcc_deriv_all.level)) ,(1-is.na(sign_diff_hcc_deriv_all.level))))
cancer.region.pair.comparability =  as.matrix(1-is.na(pair.comparability.level))

cancer.region = as.data.frame(((cancer.region.1 + cancer.region.2 )>0) + 0)*cancer.region.pair.comparability
cancer.region$time = consensus.reference.df.hbv.hcc.short$time

df_long <- melt(cancer.region, id.vars = "time")
colnames(df_long) <- c("time", "variable", "value")

x11(title = "cancer region")
pp = ggplot(df_long, aes(x = time, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_minimal() +
    theme(
        axis.text.x = element_text(hjust = 1,size=17),
        axis.text.y = element_text(hjust = 1,size=17),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20)
          ) + 
    labs(title = "", fill = "value")
print(pp)

png("canceregion.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()

cancer.pattern            = ifelse(cancer.region[,-grep("time", names(cancer.region))] == 1,1,NA)*consensus.reference.df.hbv.hcc.short[,-grep("time|population", names(consensus.reference.df.hbv.hcc.short))]
complement.cancer.pattern = ifelse(cancer.region[,-grep("time", names(cancer.region))] == 1,NA,1)*consensus.reference.df.hbv.hcc.short[,-grep("time|population", names(consensus.reference.df.hbv.hcc.short))]

cancer.pattern$time = cancer.region$time
complement.cancer.pattern$time = cancer.region$time

df_long_cancer.pattern <- melt(cancer.pattern, id.vars = "time")
colnames(df_long_cancer.pattern) <- c("time", "variable", "value")

df_long_complement.cancer.pattern <- melt(complement.cancer.pattern, id.vars = "time")
colnames(df_long_complement.cancer.pattern) <- c("time", "variable", "value")

x11(title = "cancer pattern")
pp =  ggplot(df_long_cancer.pattern, aes(x = time, y = value)) + geom_line(col = "red") + facet_wrap(variable~., scale = "free", ncol = 3) + theme_classic()
print(pp)

nocancer.pattern = ifelse(cancer.region[,-grep("time", names(cancer.region))] ==1,1,NA)*alignment.all.to.hcc.df.consensus[,-grep("time|population", names(alignment.all.to.hcc.df.consensus))]

complement.nocancer.pattern = ifelse(cancer.region[,-grep("time", names(cancer.region))] ==1,NA,1)*alignment.all.to.hcc.df.consensus[,-grep("time|population", names(alignment.all.to.hcc.df.consensus))]

           nocancer.pattern$time = cancer.region$time
complement.nocancer.pattern$time = cancer.region$time

nocancer.pattern.long =  melt(nocancer.pattern, id.vars = "time")
  colnames(nocancer.pattern.long) <- c("time", "variable", "value")
complement.nocancer.pattern.long =  melt(complement.nocancer.pattern, id.vars = "time")
colnames(complement.nocancer.pattern.long) <- c("time", "variable", "value")

nocancer.pattern.long$origin = "nocancer"
df_long_complement.cancer.pattern$origin = "cancer"
df_long_cancer.pattern$origin = "cancer"
complement.nocancer.pattern.long$origin = "nocancer" 

cancer.nocancer.pattern.long =  rbind(nocancer.pattern.long, df_long_cancer.pattern)
cancer.nocancer.pattern.long.nona = subset(cancer.nocancer.pattern.long, !is.na(value))
nn = names(table(cancer.nocancer.pattern.long.nona$variable))[(table(cancer.nocancer.pattern.long.nona$variable)>=0)] #min two years "patterns"

cancer.nocancer.pattern.long. = subset(cancer.nocancer.pattern.long, is.element(variable,nn))
df_long_complement.cancer.pattern.nona = subset(df_long_complement.cancer.pattern, is.element(variable,nn)) 
complement.nocancer.pattern.long.nona = subset(complement.nocancer.pattern.long, is.element(variable,nn))

cancer.nocancer.pattern.long.[which(cancer.nocancer.pattern.long.$variable == "rna"),"value"] = 10^(cancer.nocancer.pattern.long.[which(cancer.nocancer.pattern.long.$variable == "rna"),"value"])

x11(title = "cancer/no cancer pattern COMPL, vars")
pp =  ggplot(data = cancer.nocancer.pattern.long., aes(x = time, y = value, lty =  as.factor(origin), col = as.factor(origin) )) + geom_line(size = 1) + scale_y_log10()  +
    geom_line(data = df_long_complement.cancer.pattern.nona, aes(x = time, y = value, lty =  as.factor(origin) ), col = "grey", size = 1 ) + scale_y_log10()  +
geom_line(data = complement.nocancer.pattern.long.nona, aes(x = time, y = value, lty =  as.factor(origin) ), col = "grey", size = 1 ) +
    facet_wrap(variable~., scale = "free", ncol = 1) + scale_color_manual(values = c("nocancer" = "green", "cancer" = "red")) + geom_vline(xintercept =  seq(-10,30, by = 1), size = 0.1, col = "grey") +
    theme_classic()+ scale_y_log10() 
print(pp)

png("hcc_pattern_comparison_class_v_sense.png", width = 5, height = 25, res = 300, units = "in")
print(pp)
dev.off()

x11(title = "cancer/no cancer pattern COMPL, vars")
pp =  ggplot(data = cancer.nocancer.pattern.long., aes(x = time, y = value, lty =  as.factor(origin), col = as.factor(origin) )) + geom_line(size = 1) + scale_y_log10()  +
    geom_line(data = df_long_complement.cancer.pattern.nona, aes(x = time, y = value, lty =  as.factor(origin) ), col = "grey", size = 1 ) + scale_y_log10()  +
geom_line(data = complement.nocancer.pattern.long.nona, aes(x = time, y = value, lty =  as.factor(origin) ), col = "grey", size = 1 ) +
    facet_wrap(variable~., scale = "free", ncol = 5) + scale_color_manual(values = c("nocancer" = "green", "cancer" = "red")) + geom_vline(xintercept =  seq(-10,30, by = 1), size = 0.1, col = "grey") +
    theme_classic() + theme(xlim = c(-5, 20)) + scale_y_log10() 
print(pp)

png("hcc_pattern_comparison_class_h_sense.png", width = 12.5, height = 10, res = 300, units = "in")
print(pp)
dev.off()

message("preparing to plot regions")

hcc.hbv.alignment.withtime.list  = lapply(hcc.hbv.alignment.list, function(m) {m$time = consensus.reference.df.hbv.hcc.short$time; m })
hcc.hbv.alignment.withtime.list.noflat  = lapply(hcc.hbv.alignment.withtime.list, remove_flat_ends)
hcc.hbv.alignment.withtime.list.noflat =  lapply(1:length(hcc.hbv.alignment.withtime.list.noflat), function(ii){ x = hcc.hbv.alignment.withtime.list.noflat[[ii]]; x$id = names(hcc.hbv.alignment.withtime.list.noflat)[[ii]];x })
names(hcc.hbv.alignment.withtime.list.noflat) = names(hcc.hbv.alignment.withtime.list)
hcc.hbv.alignment.withtime.list.noflat.df = do.call(rbind, hcc.hbv.alignment.withtime.list.noflat)  
hcc.hbv.alignment.withtime.list.noflat.df$population = "hbv hcc"
hcc.hbv.alignment.withtime.list.noflat.df$apri =  with(hcc.hbv.alignment.withtime.list.noflat.df,(AST/50)*100/pla)
hcc.hbv.alignment.withtime.list.noflat.df$cd4.cd8.ratio =  with(hcc.hbv.alignment.withtime.list.noflat.df, 100*cd4/cd8) 

hcc.hbv.alignment.withtime.list.noflat.df.long = as.data.frame(hcc.hbv.alignment.withtime.list.noflat.df %>%
                                                               pivot_longer(cols = names(hcc.hbv.alignment.withtime.list.noflat.df)[-grep("time|population|id",
                                                                                   names(hcc.hbv.alignment.withtime.list.noflat.df))],
                                                                            names_to = "variable", values_to = "value")
)

##all to hcc
alignment.all.to.hcc.withtime.list.noflat = lapply(alignment.all.to.hcc.list, function(m) {m$time = consensus.reference.df.hbv.hcc.short$time; m }) ## ya no esta flat
alignment.all.to.hcc.withtime.list.noflat =  lapply(1:length(alignment.all.to.hcc.withtime.list.noflat), function(ii){ x = alignment.all.to.hcc.withtime.list.noflat[[ii]]; x$id = names(alignment.all.to.hcc.withtime.list.noflat)[[ii]];x })

names(alignment.all.to.hcc.withtime.list.noflat) = names(alignment.all.to.hcc.withtime.list.noflat)
alignment.all.to.hcc.withtime.list.noflat.df = do.call(rbind, alignment.all.to.hcc.withtime.list.noflat)
alignment.all.to.hcc.withtime.list.noflat.df$population = "all -al to hcc hbv" 
alignment.all.to.hcc.withtime.list.noflat.df$apri =  with(alignment.all.to.hcc.withtime.list.noflat.df,(AST/50)*100/pla)
alignment.all.to.hcc.withtime.list.noflat.df$cd4.cd8.ratio =  with(alignment.all.to.hcc.withtime.list.noflat.df, 100*cd4/cd8)

alignment.all.to.hcc.withtime.list.noflat.df.long =
    as.data.frame(
        alignment.all.to.hcc.withtime.list.noflat.df %>%
pivot_longer(cols = names(alignment.all.to.hcc.withtime.list.noflat.df)[-grep("time|population|id",
                    names(alignment.all.to.hcc.withtime.list.noflat.df))],
             names_to = "variable", values_to = "value")
)

alignments.combined =  rbind(hcc.hbv.alignment.withtime.list.noflat.df.long, alignment.all.to.hcc.withtime.list.noflat.df.long)

alignments.combined[which(alignments.combined$variable == "rna"), "value"] = 10^alignments.combined[which(alignments.combined$variable == "rna"), "value"] 
consensus.reference.df.hbv.hcc.short$rna = 10^consensus.reference.df.hbv.hcc.short$rna
consensus.reference.df.allpat.short$rna = 10^consensus.reference.df.allpat.short$rna

consensus.references.df.short = rbind(consensus.reference.df.hbv.hcc.short, consensus.reference.df.allpat.short)

consensus.references.df.short.long = as.data.frame(consensus.references.df.short%>%
pivot_longer(cols = names(consensus.references.df.short)[-grep("time|population|id",
                    names(consensus.references.df.short))],
             names_to = "variable", values_to = "value"))

consensus.references.df.short.long$id = "all" 


source("plots_trajectories_consensus_cancerregion.r")

##Table numbers:

message("Table numbers: Cirrhosis ever - all set");table(unique(alignments.combined[,c("id","cirrhosis.ever")])$cirrhosis.ever)
message("Table numbers: Cirrhosis ever - HCC calib set");table(unique(subset(alignments.combined, is.element(id,calib_ids.ever.hcc))[,c("id","cirrhosis.ever")])$cirrhosis.ever)
message("Table numbers: Cirrhosis ever - HCC fullset");table(unique(subset(alignments.combined, is.element(id,ids.ever.hcc))[,c("id","cirrhosis.ever")])$cirrhosis.ever)


##source("tumour_reconstruction.r")

alignment.all.to.hcc.list.noflat.ALL = alignment.all.to.hcc.list.ALL
alignment.all.to.hcc.list.noflat.ALL = lapply(alignment.all.to.hcc.list.noflat.ALL, function(m){
    m$apri = with(m, (AST/50)*100/pla)
    m$cd4.cd8.ratio = with(m, 100*cd4/cd8)
    m
})

list_of_series. = lapply(alignment.all.to.hcc.list.noflat.ALL, function(x) x[, -grep("time|id", names(alignment.all.to.hcc.list.noflat.ALL[[1]]))] )
names(list_of_series.) = names(alignment.all.to.hcc.list.noflat.ALL)

cancer.region.nas = cancer.region
cancer.region.nas[,-grep("time", names(cancer.region.nas) )] = ifelse(cancer.region.nas[,-grep("time", names(cancer.region.nas) )] ==1,100,1)

x11(title = "Cancer region NA");image(as.matrix(cancer.region.nas[,-grep("time", colnames(cancer.region.nas))]))

list_of_series.all =  lapply(list_of_series., function(m) m*cancer.region.nas[,-grep("time", names(cancer.region))])

rm(list_of_series.); malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));gc()

list_of_series.all.plus.hcc.norep = list_of_series.all
message("about to compute distance for the big set")
list_of_series.all.plus.hcc.norep= list_of_series.all.plus.hcc.norep[c(1:1235,1238:1910)]
D.dist = compute_dtw_dist.alternative(list_of_series.all.plus.hcc.norep)
D = as.matrix(D.dist)
write.table(as.matrix(D),"dtw_dist_combined")

D = as.matrix(read.table("dtw_dist_combined"))
rownames(D) = names(list_of_series.all.plus.hcc.norep) 
colnames(D) = names(list_of_series.all.plus.hcc.norep) 

D_fw_tri = ifelse(lower.tri(D),1,NA)*D 

ideh = which(is.element(rownames(D_fw_tri),calib_ids.ever.hcc))
idrow.test.set =  which(is.element(rownames(D_fw_tri),c(test_ids.ever.hcc, test_ids.nohcc)))

D_fw_cancer = D_fw_tri[idrow.test.set,ideh]      

D_fw_cancer.min = apply(D_fw_cancer, 1, function(x) min(x, na.rm= T))
x11(title = "Histogram of listed");hist(D_fw_cancer.min,100)

qq = quantile(D_fw_cancer.min, na.rm = T)
message("distances between cancer y no cancer ts"); print(qq)
x11(title = "Distances"); hist(as.numeric(D_fw_cancer[D_fw_cancer<Inf]), 200)

D_fw_cancer.df =  as.data.frame(D_fw_cancer)
names(D_fw_cancer.df) = rownames(D_fw_tri)[ideh]
D_fw_cancer.df$id = rownames(D_fw_cancer)
D_fw_cancer.df$min= apply(D_fw_cancer, 1, function(x) min(x, na.rm= T))
D_fw_cancer.df.comparable = subset(D_fw_cancer.df, min !=Inf)

source("cum_incidence.r")
