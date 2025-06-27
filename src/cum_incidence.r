message("Entering predictive power")

D_fw_cancer.short = D_fw_cancer.df.comparable[,c("id", "min")]
D_fw_cancer.short = subset(D_fw_cancer.short, min !=Inf)
D_fw_cancer.short$distance = D_fw_cancer.short$min 
D_fw_cancer.short$distance[D_fw_cancer.short$distance == 0] = 1e-16
D_fw_cancer.short$distance = log10(D_fw_cancer.short$distance)
D_fw_cancer.short$proximity = 100*with(D_fw_cancer.short, 1/(1+distance)) 

x11(title = "distance, proximity transform")
par(mfrow = c(2,2))
hist(D_fw_cancer.short$min,100)
hist(D_fw_cancer.short$distance,100)
hist(D_fw_cancer.short$proximity)

data = subset(read.table(paste(data.path1,"/lab_wide.baseline.NOhcc.data", sep = "")))
data$hcc[is.na(data$hcc) ] = 0 
data$before.hcc[is.na(data$before.hcc) ] = 1 
data = subset(data, before.TFV_TAF==0)

lab_wide.baseline.NOhcc.data = subset(data, before.hcc  == 1 )

data.distance.listed =  merge(D_fw_cancer.short, lab_wide.baseline.NOhcc.data, by.id =  "TRUE", all.x = T)

data.distance.listed$apri =  with(data.distance.listed, (AST/50)*100/pla )

data.distance.listed.. = 
unique(data.frame(data.distance.listed %>%
    pivot_longer(cols = names(data.distance.listed)[-grep("time|id|cohort|lab_d|lab_u|enrol_d|sex|transmission_mode|ethnic|hcc_d",names(data.distance.listed))], names_to = "variable", values_to = "value")
    ))

data.distance.listed. = unique(data.distance.listed..[!is.na(data.distance.listed..$value),])

data.distance.listed <- data.distance.listed. %>%
  group_by(id, time,cohort,lab_d,lab_u,enrol_d,sex,transmission_mode,ethnic ,variable) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

data.distance.listed.wide <- as.data.frame(pivot_wider(data.distance.listed, names_from = variable, values_from = value, values_fill = NA_real_))

write.table(data.distance.listed.wide, "data_distance_listed_wide")

data.distance.listed.wide = read.table("data_distance_listed_wide")

points.for.age.f = function(x){ (x>29)*2 + (x>39)*2 + (x>49)*2 + (x>59)*2 + (x>69)*2 + (x>79)*2 }
points.for.pla.f = function(x){ (x<200)*6 + (x<100)*3   }

data.distance.listed.wide$PAGEB =  with(data.distance.listed.wide,
                                        points.for.age.f(age.updated) +
                                        6*ifelse(sex == "Male",1,0) +
                                        points.for.pla.f(pla)
                                        )

data.distance.listed.wide.page = subset(data.distance.listed.wide, !is.na(PAGEB) )

cirrh.ever.data = do.call(rbind, by(data.distance.listed.wide.page, data.distance.listed.wide.page$id, function(x){ data.frame(id = x$id[1],
                                                                                                                               cirrh.ever.count = sum(x$apri >=2, na.rm = T),
                                                                                                                                     cirrh.ever = (sum(x$apri >=2, na.rm = T)>0) + 0 ) } )
                          )

data.distance.listed.wide.page. = do.call(rbind, by(data.distance.listed.wide.page,
                                                    data.distance.listed.wide.page$id,
                                                    function(x) {
                                                        subset(x, time == max(x$time)) }))

data.distance.listed.wide.page. = merge(data.distance.listed.wide.page., cirrh.ever.data, by = "id")
 
aa = data.distance.listed.wide.page.[,c("id","PAGEB", "proximity", "sex","hcc", "apri", "ethnic", "distance", "cirrh.ever.count", "cirrh.ever")]

bb = with(aa, is.element(id, test_ids.ever.hcc))
cc = with(aa, is.element(id, calib_ids.ever.hcc))

aa$test =  bb
aa$calbi = cc

x11(title = " ")
plot(aa$PAGEB, log(aa$proximity))

aa$cirrh =  ifelse(aa$apri>=2,1,0)
aa$hcc[is.na(aa$hcc)] = 0

aa$cirrh.exp = 0
aa$cirrh.exp[aa$cirrh.ever.count>=1] = 1
aa$cirrh.exp[aa$cirrh.ever.count>=4] = 2
aa$cirrh.exp[aa$cirrh.ever.count>=8] = 3

print("table hcc: ");print(table(aa[,c("hcc")]))
print("table cirrhosis: ");print(table(aa[,c("hcc", "cirrh.ever")]))
print("table cirrhosis: ");print(table(aa[,c("hcc", "cirrh.ever","test")]))

x11(title = " ")
pp = ggplot(aa, aes(x = as.factor(hcc), y = as.numeric(proximity), col =  as.factor(hcc), pch = test )) +
    geom_point(size = 3, alpha = 0.5)  + ## coord_trans(y = "log2")
    geom_point(data = subset(aa, test == T), size = 3, col = "blue")  +
    geom_hline(yintercept = c(12.75,13.4), col = "red", lty = 2)  
print(pp)

aa$proximity.group = 0
aa$proximity.group[aa$proximity>=12.75 ] = 1
aa$proximity.group[aa$proximity>=13.4 ] = 2


aa$PAGEB.group = "10" 
aa$PAGEB.group[aa$PAGEB>10] = "10-17"
aa$PAGEB.group[aa$PAGEB>17] = "18"

# Filter only rows with event = 1
events_only <- aa %>%
  filter(hcc == 1) %>%
  group_by(proximity.group, test, PAGEB.group) %>%
      mutate(y = row_number()) %>%  # Assign y position for each event
  ungroup()

## Plot stacked circles
x11(title = "circles")
ggplot(events_only, aes(x = interaction(proximity.group, test), y = y, fill = as.factor(PAGEB.group),
                        )) +
     geom_point(shape =  21, size = 6, stroke = 0.3) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "Stacked Circles for Binary Events",
       x = "Proximity Group / Test", y = "", fill = "Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

##Let's go back to data
lab_wide.baseline.NOhcc.data$time.since.tenofovir.start =  with(lab_wide.baseline.NOhcc.data, time.year - TFV_TAF_start.time.year)
todo =  merge(lab_wide.baseline.NOhcc.data, aa, by = "id", suffixes = c("","",""), all.y = T)

todo.for.survival = do.call(rbind,by(todo, todo$id, function(m){    
    
    data.frame(id =  m$id[1],
               time = max(m$time.since.tenofovir.start),
               event = m$hcc[1],
               proximity.group = m$proximity.group[1],
               proximity =  m$proximity[1],
               PAGEB = m$PAGEB.group[1],
               PAGEB.group = m$PAGEB.group[1],
               test = m$test[1],
               cirrh.ever =  m$cirrh.ever[1]
               )
} ))

fit <- survfit(Surv(time, event) ~ proximity.group, data  =todo.for.survival)
x11(title = "by proximity")
pp = ggsurvplot(fit,
           fun = "event",  # 1 - survival
           conf.int = F,
           risk.table = TRUE,
           ggtheme = theme_minimal(),
           title = "Cumulative Incidence by Group", ylim =  c(0, 0.075))

pp = pp$plot +
    theme(
        legend.position = "right",
        axis.text.x = element_text(hjust = 1,size=17),
        axis.text.y = element_text(hjust = 1,size=17),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size = 20)
          )


print(pp)

png("cuminc_by_proximity.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()


fit <- survfit(Surv(time, event) ~ PAGEB.group, data  =todo.for.survival)
                                        # Plot using ggsurvplot
x11(title = "by PAGE")
pp = ggsurvplot(fit,
           fun = "event",  # 1 - survival
           conf.int = FALSE,
           risk.table = TRUE,
           ggtheme = theme_minimal(),
           title = "Cumulative Incidence by Group", ylim =  c(0, 0.075))


pp = pp$plot +
    theme(
        legend.position = "right",
        axis.text.x = element_text(hjust = 1,size=17),
        axis.text.y = element_text(hjust = 1,size=17),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size = 20)
          )

print(pp)

png("cuminc_by_pageb.png", width = 7, height = 7, res = 300, units = "in");print(pp);dev.off()
