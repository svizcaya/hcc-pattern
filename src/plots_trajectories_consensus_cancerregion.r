
df_long_complement.cancer.pattern.nona$population =  "hbv hcc"
df_long_complement.cancer.pattern.nona$cirrhosis.ever = "agg"
complement.nocancer.pattern.long.nona$population = "all -al to hcc hbv"
complement.nocancer.pattern.long.nona$cirrhosis.ever = "agg"
valid = df_long_complement.cancer.pattern.nona
valid$value[!is.na(valid$value)] = 1 
valid = unique( valid[,c("time","value", "variable")] )
valid[is.na(valid$value),"value"] = 0
valid$value = 1- valid$value

highlight_regions. <- valid %>%
  arrange(variable, time) %>%
  group_by(variable) %>%
  mutate(group = cumsum(c(TRUE, diff(value) != 0))) %>%  ## identify runs
  filter(value == 1) %>%
  group_by(variable, group) %>%
  summarise(
    xmin = min(time),
    xmax = max(time),
    ymin = 0,
    ymax = Inf,
    .groups = "drop"
  )


alignments.combined. = do.call(rbind, by(alignments.combined, alignments.combined$variable, function(m) {    
    outlier = quantile(m$value, 1-1e-6, na.rm = T);   
    outlier.indices = which(abs(m$value)>=abs(outlier))
    if(length(outlier.indices) >=0 ) {m =  m[-outlier.indices,]}
    m
} ) )

alignments.combined.init =  alignments.combined
alignments.combined =  alignments.combined.

x11(title = "alignment.all.to.hcc")
pp = ggplot(data = alignments.combined, aes(x = as.numeric(time), y = as.numeric(value), grouping = as.factor(id))) +
                             geom_line(size = 0.5, alpha = 0.5, aes(col = as.factor(population))) + 
    theme_minimal() +   facet_wrap(~variable, ncol = 6, scale = "free") +  scale_y_log10() + 
    geom_line(data = consensus.references.df.short.long, size = 0.5,  aes(lty =  as.factor(population)) ) +
  geom_rect(data = as.data.frame(highlight_regions.), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "red", alpha = 0.2, inherit.aes = FALSE) 
print(pp)    

png("alignment_all_and_hcc_to_hcc_consensus_traj_consensusseries.png", width = 10, height = 15, res = 300, units = "in")
print(pp)
dev.off()

alignments.combined.cirr = subset(alignments.combined, variable == "apri")
cirrhosis.ever = do.call(rbind, by(alignments.combined.cirr, alignments.combined.cirr$id, function(x) data.frame(cirrhosis.ever = sum(x$value>=2, na.rm = T)>0 + 0, id = x$id[1] ) ))
alignments.combined.cirr = merge(alignments.combined.cirr, cirrhosis.ever, by = "id", suffixes = c("",""))
alignments.combined.cirr = alignments.combined.cirr[,unique(names(alignments.combined.cirr))]

pp = ggplot(data = alignments.combined.cirr, aes(x = as.numeric(time), y = as.numeric(value), grouping = as.factor(id))) +
    geom_line(size = 0.25, alpha = 0.8, aes(col = as.factor(population))) + scale_y_log10() + 
    geom_hline(yintercept = c(0.5,1,1.5,2), col = "red", lty = 3, lwd = 1) +
    theme_minimal() + facet_grid(.~cirrhosis.ever) + ##scale_colour_viridis_d()+
    stat_summary(fun = median, geom = "line", aes(group =  as.factor(population), lty = as.factor(population) ), size = 0.75, col = "black") +
  geom_rect(data = subset(as.data.frame(highlight_regions.), variable == "apri"), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "red", alpha = 0.2, inherit.aes = FALSE) 

x11(title = "APRI Cirrhosis");print(pp)
png("APRI_CIRRHOSISever_alignment_all_and_hcc_to_hcc_consensus_traj.png", width = 10, height = 10, res = 300, units = "in")
print(pp)
dev.off()

alignments.combined =  merge(alignments.combined, cirrhosis.ever, by = "id", suffixes = c("",""))


x11(title = "alignment.all.to.hcc no cirr")
pp = ggplot(data = alignments.combined, aes(x = as.numeric(time), y = as.numeric(value), grouping = as.factor(id))) +
                             geom_line(size = 0.25, alpha = 0.8, aes(col = as.factor(population))) + scale_y_log10()+
    theme_minimal() +   facet_wrap(~variable, ncol = 3, scale = "free") + 
    stat_summary(fun = median, geom = "line", aes(group =  as.factor(population), lty = as.factor(population) ), col = "black")+
  geom_rect(data = as.data.frame(highlight_regions.), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "red", alpha = 0.2, inherit.aes = FALSE)  
print(pp)
png("alignment_all_and_hcc_to_hcc_consensus_traj.png", width = 10, height = 15, res = 300, units = "in")
print(pp)
dev.off()

alignments.combined.zoom =  subset(alignments.combined, is.element(variable, c("ALB",
                                                                               "ALT",
                                                                               "apri",
                                                                               "AST",
                                                                               "BIL",
                                                                               "cd4",
                                                                               "cd8",
                                                                               "CHOL",
                                                                               "CRE",
                                                                               "GLUC",
                                                                               "hem",
                                                                               "leu",
                                                                               "pla",
                                                                               "rna",
                                                                               "TRIG",
                                                                               "cd4.cd8.ratio"
                                                                               )))


highlight_regions.. =  subset(as.data.frame(highlight_regions.), is.element(variable, names(table(alignments.combined.zoom$variable))))
highlight_regions..$variable =  as.character(highlight_regions..$variable)

alignments.combined.zoom = subset(alignments.combined.zoom, time >=0 )

## -- plot in the report
pp =   ggplot(alignments.combined.zoom, aes(x = time, y = value, col = as.factor(population)     ))+
    theme_minimal() +
    theme(legend.position =  "bottom",
          axis.text.x = element_text(hjust = 1,size=13),
          axis.text.y = element_text(hjust = 1,size=13),
          strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size = 20)
          ) + labs(x =  "Time (years)")+

    facet_wrap(.~variable, ncol = 3, scale = "free") +
    stat_summary(fun = median, geom = "line", fun.args = list(na.rm = TRUE)) +
    scale_y_log10() +
     geom_rect(data = subset(highlight_regions.,is.element(variable, unique(alignments.combined.zoom$variable))), 
             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
             fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(name = "",values = c("green", "red"), labels = c("without HCC diagnosis","with HCC diagnosis")) 

x11(title = "no traje HCC, no hcc"); print(pp)
png("report_zoom_post_alignment_consensus_hcc_notraj.png", width = 10, height = 12.5, res = 300, units = "in");print(pp);dev.off()

x11(title = "Zoom alignment.all.to.hcc no cirr TRAJ")
pp = ggplot(data = alignments.combined.zoom, aes(x = as.numeric(time), y = as.numeric(value), grouping = as.factor(id))) +
                             geom_line(size = 0.25, alpha = 0.8, aes(col = as.factor(population))) + scale_y_log10()+
    theme_minimal() +   facet_wrap(variable~cirrhosis.ever, ncol = 6, scale = "free") + 
     stat_summary(fun = median, geom = "line", aes(group =  as.factor(population), lty = as.factor(population) ), col = "black") +
  geom_rect(data = subset(as.data.frame(highlight_regions..), is.element(variable, unique(alignments.combined.zoom$variable)) ), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "red", alpha = 0.2, inherit.aes = FALSE)
print(pp)
png("zoom_CIRRHOSIS_post_alignment_consensus_trajectories.png", width = 12.5, height = 15, res = 300, units = "in");print(pp);dev.off()

x11(title = "no traje CIRR HCC"); print(pp)
pp = ggplot(alignments.combined.zoom, aes(x = time, y = value, col = as.factor(population),lty = as.factor(cirrhosis.ever) ))+
    theme_minimal() +  theme(legend.position =  "bottom") + facet_wrap(.~variable, ncol = 3, scale = "free") +
    stat_summary(fun = median, geom = "line") + scale_color_manual(values = c("green", "red")) +scale_y_log10() +
    geom_rect(data = subset(highlight_regions.,is.element(variable, unique(alignments.combined.zoom$variable))), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "red", alpha = 0.1, inherit.aes = FALSE)
png("zoom_CIRRHOSIS_post_alignment_consensus.png", width = 12.5, height = 15, res = 300, units = "in");print(pp);dev.off()

