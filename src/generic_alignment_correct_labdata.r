
alignment.list = consensus_alignment(series_list = lab.ts.after.tdf.interpolated.pp)

alignment.ref = alignment.list$reference
alignment.init.ref = alignment.list$ref.init
alignment.al = alignment.list$aligned

reference.df = melt(as.data.frame(alignment.ref), id.vars = "time")
reference.init.df = melt(as.data.frame(alignment.init.ref), id.vars = "time")

ids = names(alignment.al)
alignment.al.withids =  lapply(1:length(ids), function(ii) {x = alignment.al[[ii]];x$id = ids[[ii]]; x } )

alignment.al.df = do.call(rbind, alignment.al.withids)

write.table(reference.df, file.name.consensus)
write.table(alignment.al.df, file.name.alignment)

x11(title = "cons reference");
pp = ggplot(reference.df, aes(x= time, y= value)) +
    facet_wrap(variable~., scale = "free", ncol = 3) +
    geom_line(col = "green") + geom_point(col = "green", size = 0.01)+
    scale_y_log10()+
    theme_classic() 
print(pp)

png(file.name.figure.consensus, width = 7, height = 12, res = 300, units = "in");print(pp);dev.off()

##plotear el alignment:
x11(title = "pre todo");plot_feature_all.2(ts_list = lab.ts.after.tdf.interpolated.pp, colapsar = T, ref.df = reference.df) 
x11(title = "post todo"); plot_feature_all.2(ts_list = alignment.al, colapsar = T, ref.df = reference.df) 

##looking at a random id:
id.check = names(alignment.al)[sample(1:length(unique(names(alignment.al))),1)]

alignment.al.df.= alignment.al.df
alignment.al.df.$origin = "alignment"

alignment.al.df. = alignment.al.df.[,names(pp.data.imputed.)]
data.pp.to.impute.tetris.wide. = data.pp.to.impute.tetris.wide
data.pp.to.impute.tetris.wide.$origin = "to.impute"
data.pp.all = rbind(pp.data.imputed., data.pp.to.impute.tetris.wide., alignment.al.df.)

data.pp.all.long <- data.pp.all %>%
    pivot_longer(cols = names(data.pp.all)[-grep("time|id|origin",names(data.pp.all))], names_to = "variable", values_to = "value")

reference.df. = reference.df
reference.df.$origin = "reference"

pp = ggplot(as.data.frame(subset(data.pp.all.long, id ==  id.check)), aes(x = time, y = value, col = origin )) +
    facet_wrap(variable~id, scale = "free", ncol = 3) +
    geom_point(alpha = 0.75, size = 0.75) +  ##+ theme_classic()
    geom_line(alpha = 0.5, size = 0.75) +  ##+ theme_classic()
    geom_line(data = as.data.frame(reference.df.) )+
  theme_classic()

x11(title = "just and imputation case");print(pp)

png(file.name.imputation_outcome_example, width = 10, height = 10, res = 300, units = "in")

pp = ggplot(subset(as.data.frame(subset(data.pp.all.long, id ==  id.check)), origin !="imputed"), aes(x = time, y = value, col = origin )) +
    facet_wrap(variable~id, scale = "free", ncol = 3) +
    geom_point(alpha = 0.75, size = 0.75) +  
    geom_line(alpha = 0.5, size = 0.75) + 
    geom_line(data = as.data.frame(reference.df.), size =1 )+scale_y_log10()+
    theme_classic()+
    geom_vline(xintercept = seq(-5,40, by = 1), col = "grey", size = 0.25)

x11(title = "just and imputation case - no imp");print(pp)

png(file.name.original_to_aligned, width = 7, height = 12, res = 300, units = "in");print(pp);dev.off()
