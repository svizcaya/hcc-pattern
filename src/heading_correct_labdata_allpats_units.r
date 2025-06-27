
##source("units_corrector.r");message("done with the units_corrector")
    
lab.read = read.table(paste(data.path1,"/lab",sep = ""))
lab.read$lab_d = as.Date(lab.read$lab_d)


lab = lab.read

## First, spread the data to a wider format with pivot_wider
lab_wide <- lab %>%
  pivot_wider(
    names_from = lab_id,
    values_from = lab_r
   # This handles duplicates by summing
  )

# Output the reshaped data
lab_wide.data = as.data.frame(lab_wide)

#####---
## lab_wide.data es la de lab
fibroscan_validated = read.dta13(paste(data.path0,"/fibroscan_validated.dta", sep = ""))
fibroscan = fibroscan_validated[,c("projectid", "cohort", "date", "median_lsm")]
fibroscan <- fibroscan %>% rename(fibroscan_median_lsm = median_lsm)
fibroscan <- fibroscan %>% rename(lab_d = date)
baseline_validated = read.dta13(paste(data.path0,"/baseline_validated.dta", sep = ""))
baseline = baseline_validated[,c("projectid", "enrol_d", "sex", "transmission_mode", "age_at_enrol", "cohort", "ethnic")]
baseline = subset(baseline, !is.na(enrol_d))
baseline$enrol.time.year = as.numeric(baseline$enrol_d)/365.25 + 1970
rm(lab_wide);rm(lab);malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0))); malloc.trim();gc()
baseline$cohort[baseline$cohort == "Eurosida"] = "Euro" 
fibroscan$cohort[fibroscan$cohort == "Eurosida"] = "Euro"
lab_wide.baseline.data = merge(lab_wide.data, baseline, by = c("projectid","cohort"), suffixes = c("", ""))

lab_wide.baseline.data. = merge(lab_wide.baseline.data, fibroscan, by = c("projectid", "lab_d", "cohort"), all.x =  T)
lab_wide.baseline.data.$time.year = as.numeric(as.Date(lab_wide.baseline.data.$lab_d))/365.25 + 1970 
lab_wide.baseline.data = lab_wide.baseline.data.
lab_wide.baseline.data$time.since.enrol  =  with(lab_wide.baseline.data, time.year- enrol.time.year)
lab_wide.baseline.data =  subset(lab_wide.baseline.data, time.since.enrol>= -5)
lab_wide.baseline.data$age.updated = lab_wide.baseline.data$time.since.enrol + as.numeric(lab_wide.baseline.data$age_at_enrol)

##merge with the SHCS data 
aa = unlist(lapply(names(lab_wide.baseline.data), function(nn){
   ## print(nn)
    bb = lab_wide.baseline.data[!is.na(lab_wide.baseline.data[,nn]),]
    aux = table(bb[,"cohort"])
    aux. = names(aux)
test = (sum(is.element(aux., "SHCS")) ==0)
if(test){out = nn}else{out =NA} 
    out
}))
names.pull.shcs = aa[!is.na(aa)]

rm(aa);malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0))); malloc.trim();gc()

lab12.shcs = read.table(paste(data.path1,"/lab12", sep = ""))[,c("id","labdate",names.pull.shcs)] 
lab12.shcs = lab12.shcs %>% rename(projectid = id)
lab12.shcs = lab12.shcs %>% rename(lab_d = labdate)
lab12.shcs$cohort = "SHCS"

###error levels
lab12.shcs$BIL[as.numeric(lab12.shcs$BIL)>600] = NA; lab12.shcs$BIL[as.numeric(lab12.shcs$BIL)<=0] = NA 
lab12.shcs$ALB[as.numeric(lab12.shcs$ALB)>70] = NA; lab12.shcs$ALB[as.numeric(lab12.shcs$ALB)<15] = NA 
lab12.shcs$ALT[as.numeric(lab12.shcs$ALT)>2000] = NA; lab12.shcs$ALT[as.numeric(lab12.shcs$ALT)<=0] = NA 
lab12.shcs$AST[as.numeric(lab12.shcs$AST)>500] = NA; lab12.shcs$AST[as.numeric(lab12.shcs$AST)<=0] = NA 
lab12.shcs$CHOL[as.numeric(lab12.shcs$CHOL)>20] = NA; lab12.shcs$CHOL[as.numeric(lab12.shcs$CHOL)<1] = NA
lab12.shcs$TRIG[as.numeric(lab12.shcs$TRIG)>40] = NA; lab12.shcs$TRIG[as.numeric(lab12.shcs$TRIG)<=0] = NA
lab12.shcs$CRE[as.numeric(lab12.shcs$CRE)>2000] = NA; lab12.shcs$CRE[as.numeric(lab12.shcs$CRE)<=0] = NA
lab12.shcs$GLUC[as.numeric(lab12.shcs$GLUC)>40] = NA; lab12.shcs$GLUC[as.numeric(lab12.shcs$GLUC)<=0] = NA

block.aux. = lab_wide.baseline.data[, names(lab12.shcs)]
block.aux = subset(block.aux., cohort == "SHCS")[,c("projectid", "lab_d", "cohort")]
block.aux.shcs.pull = merge(block.aux, lab12.shcs, by =  c("projectid", "lab_d", "cohort"), all.x = T, all.y = T)  
rr.block.no.shcs = subset(block.aux., cohort != "SHCS")
rr.block =  rbind(rr.block.no.shcs, block.aux.shcs.pull)
block.aux.complement = lab_wide.baseline.data[,c("projectid","lab_d","cohort",
                                                 names(lab_wide.baseline.data)[!is.element(names(lab_wide.baseline.data),names(lab12.shcs))])]

rm(lab12.shcs); malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0))); malloc.trim();gc()

lab_wide.baseline.data = merge(block.aux.complement, rr.block, by = c("projectid", "lab_d", "cohort"))

rm(rr.block.no.shcs)
rm(rr.block)
rm(block.aux.complement)
rm(block.aux)
rm(block.aux.shcs.pull)
rm(baseline)
try(rm(adh_euro_validated))
rm(fibroscan_validated)

malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0))); malloc.trim();gc()

vars.to.ts =  c(
    "time",
    "aghbe",
    "aghbs",
    "ALB",
    "ALT",
    "AMY",
    "antihcv",
    "antihbc",
    "antihbe",
    "antihbs",
    "antihdv",
    "AST",
    "BIL",
    "cd3",
    "cd3p",
    "cd4",
    "cd4p",
    "cd8",
    "cd8p",
    "CHOL",
    "antiCMV",
    "CRE",
    "GLUC",
    "hem",
    "HBVDNA",
    "HCV_RNA",
    "rna",
    "QUI",
    "leu",
    "lym",
    "pla",
    "UPROT",
    "toxo",
    "TRIG"
    )

lab_wide.baseline.data$rna = as.numeric(lab_wide.baseline.data$rna)
lab_wide.baseline.data$rna[lab_wide.baseline.data$rna<=0] =  1e-16
lab_wide.baseline.data$rna = log10(lab_wide.baseline.data$rna)

aux = lab_wide.baseline.data[,vars.to.ts[-grep("time|age.updated", vars.to.ts)]]

test = apply(aux, 1, function(l) sum(!is.na(l)) ==0 )

setDT(lab_wide.baseline.data)
lab_wide.baseline.data = unique(lab_wide.baseline.data[!test,])

fup.time.all.iqr = quantile(by(lab_wide.baseline.data$time.since.enrol, lab_wide.baseline.data$projectid, function(x) max(x, na.rm = T)))
print(" max fup time all:");print(fup.time.all.iqr)

##-now for the subpopulations
ids.hbv.ever = unique(subset(lab_wide.baseline.data, aghbs == 1)$projectid)

print("# ids ever positive antihbc:"); print(length(ids.hbv.ever))

lab_wide.baseline.data.hbv = subset(lab_wide.baseline.data, is.element(projectid, ids.hbv.ever))
fup.time.hbv.iqr = quantile(by(lab_wide.baseline.data.hbv$time.since.enrol, lab_wide.baseline.data.hbv$projectid, function(x) max(x, na.rm = T)))
print("max fup time ever antihbc+:");print(fup.time.hbv.iqr)

lab_wide.baseline.data = lab_wide.baseline.data %>% rename(id = projectid)
lab_wide.baseline.data$time = lab_wide.baseline.data$time.since.enrol

hcc_data = read.dta13(paste(data.path0,"/hcc_data.dta", sep = ""))

hcc_data$id = hcc_data$projectid
hcc_data$cohort[hcc_data$cohort == "Eurosida"] = "Euro"

hcc_data = hcc_data[order(hcc_data$hcc_d),]
hcc_data$id = hcc_data$projectid
hcc_data.init = hcc_data
hcc_data =  hcc_data[!duplicated(hcc_data$id),]

message("merge 1")
lab_wide.baseline.NOhcc.data = merge(lab_wide.baseline.data, hcc_data, by = c("id","cohort"), suffixes = c("", ""), all.x = T)
lab_wide.baseline.NOhcc.data$hcc.time.year = as.numeric(lab_wide.baseline.NOhcc.data$hcc_d)/365.25 + 1970 
lab_wide.baseline.NOhcc.data = unique(lab_wide.baseline.NOhcc.data)

##fitst hit in "ART_validated" with a TDF medicine
ART_validated = read.dta13(paste(data.path0,"/ART_validated.dta", sep = ""))
ART_validated$start.time.year =  as.numeric(ART_validated$start_date)/365.25 + 1970
ART_validated$stop.time.year =  as.numeric(ART_validated$stop_date)/365.25 + 1970

ART_validated$TFV_TAF =  0
ART_validated$TFV_TAF[grep("TFV|TAF|TDF|Tenofovir|tenofovir|Truvada|Atripla|tdf|taf|tfv", ART_validated$art_name)] = 1

ART_validated.start.tenofovir. =  subset(ART_validated, TFV_TAF ==1)

ART_validated.start.tenofovir = do.call(rbind, by(ART_validated.start.tenofovir., ART_validated.start.tenofovir.$projectid,
                                        function(x) {  x = data.frame(x);
                                            data.frame(projectid = x$projectid,
                                                       cohort =  x$cohort,
                                                       TFV_TAF_start.time.year = min(x$start.time.year),
                                                       TFV_TAF_stop.time.year = max(x$stop.time.year) ## will just ignore in-between
                                                       ) }
                                           ) )   ##will use small in "start" for the one i create


ART_validated.start.tenofovir = ART_validated.start.tenofovir %>% rename(id = projectid)

ART_validated.start.tenofovir$cohort[ART_validated.start.tenofovir$cohort == "Eurosida"] = "Euro"

ART_validated.start.tenofovir = unique(ART_validated.start.tenofovir)

message("merge 2")
lab_wide.baseline.NOhcc.data = merge(lab_wide.baseline.NOhcc.data, ART_validated.start.tenofovir[,c("id", "cohort", "TFV_TAF_start.time.year","TFV_TAF_stop.time.year")], by = c("id","cohort"), suffixes = c("", "") )

lab_wide.baseline.NOhcc.data$before.TFV_TAF = with(lab_wide.baseline.NOhcc.data, ifelse(time.year <= TFV_TAF_start.time.year, 1,0) )

lab_wide.baseline.NOhcc.data$before.hcc = with(lab_wide.baseline.NOhcc.data, ifelse(time.year <= hcc.time.year, 1,0) )

print("# all patients with follow-up after tenofovir start");print(length(unique(lab_wide.baseline.NOhcc.data$id)))
print("Table 1: # all patients with follow-up after tenofovir start and ever had a positive aghbs");print(length(unique(subset(lab_wide.baseline.NOhcc.data, aghbs ==1)$id)))
print("Table 1: # all patients with follow-up after tenofovir start and ever had a positive aghbs and a hcc diagnosis"); print(length(unique(subset(lab_wide.baseline.NOhcc.data, ((hcc ==1) + (aghbs ==1)) ==2  )$id)))

message("upto here all info to be used combined")

lab_wide.baseline.data.hbv.NOhcc = subset(lab_wide.baseline.NOhcc.data, is.element(id, ids.hbv.ever)  )

fup.time.hbv.hcc.iqr = quantile(by(lab_wide.baseline.data.hbv.NOhcc$time.since.enrol, lab_wide.baseline.data.hbv.NOhcc$id, function(x) max(x, na.rm = T)), seq(0,1,0.05))

print("max fup time ever hcc and positive antihbc:");print(as.matrix(fup.time.hbv.hcc.iqr))

ids.atleast.Xfup.NOhcc = do.call(rbind, by(lab_wide.baseline.data.hbv.NOhcc[,c("id","time.since.enrol")],
                                                  lab_wide.baseline.data.hbv.NOhcc$id, function(x) data.frame(id 
  = x$id[1], fup = max(x$time.since.enrol, na.rm = T)) ))$id  #will pracmatically cut at 2y


print("# Table 1: of ids with or without HCC in HBV+ & long enough follow-up");print(length(ids.atleast.Xfup.NOhcc))
##now take only those with enough follow-up
lab_wide.baseline.data.hbv.NOhcc.Xfup = subset(lab_wide.baseline.data.hbv.NOhcc, is.element(id,ids.atleast.Xfup.NOhcc))

message("Table 1: aghbs+ & fup tenof and hcc");print(length( unique(subset(lab_wide.baseline.data.hbv.NOhcc.Xfup, hcc==1)$id)))

##will now check the number of registries per person(id)
registries.per.id = as.matrix(quantile(as.numeric(by(lab_wide.baseline.data.hbv.NOhcc.Xfup$time.since.enrol, lab_wide.baseline.data.hbv.NOhcc.Xfup$id, length)), seq(0,1,0.05) ))
print("distribution of number of registries among ids with NO HCC in HBV+ & long enough follow-up")
print(registries.per.id)

##n)will not exclude anyone, no need to  up to know
print("distribution across cohorts ids with or without HCC in HBV+ & long enough follow-up");table(lab_wide.baseline.data.hbv.NOhcc.Xfup[!duplicated(lab_wide.baseline.data.hbv.NOhcc.Xfup$id),]$cohort)

##checking the interval between records

#organize according to the order of the records
lab_wide.baseline.data.hbv.NOhcc.Xfup = do.call(rbind, by(lab_wide.baseline.data.hbv.NOhcc.Xfup, lab_wide.baseline.data.hbv.NOhcc.Xfup$id, function(x) {x[order(x$time.year),]}))

##naive interval between records:

aux.1 = by(lab_wide.baseline.data.hbv.NOhcc.Xfup$time.year, lab_wide.baseline.data.hbv.NOhcc.Xfup$id, function(x) {
   out =  (x[2:length(x)] - x[1:(length(x)-1)])*365.25

    out
}
)

pp.data = subset(lab_wide.baseline.data.hbv.NOhcc.Xfup, !is.na( aghbs))

aux.2  =  by(pp.data$time.year, pp.data$id, function(x) {
   out =  (x[2:length(x)] - x[1:(length(x)-1)])*365.25
    out
}
)

print("for specific e.g cd4, the interval between measurements distribution");print(quantile(unlist(aux.2), seq(0,1,0.05), na.rm = T))

lab_wide.baseline.data.hbv.NOhcc.Xfup.init = lab_wide.baseline.data.hbv.NOhcc.Xfup 
lab_wide.baseline.data.hbv.NOhcc.Xfup$before.hcc[is.na(lab_wide.baseline.data.hbv.NOhcc.Xfup$before.hcc)] = 1
lab_wide.baseline.data.hbv.NOhcc.Xfup = subset(lab_wide.baseline.data.hbv.NOhcc.Xfup, before.hcc == 1)

data.pp = lab_wide.baseline.data.hbv.NOhcc.Xfup
data.pp$rna[ data.pp$rna == -1] = NA
lab_wide.baseline.NOhcc.data$rna[lab_wide.baseline.NOhcc.data$rna == -1]= NA

write.table(lab_wide.baseline.NOhcc.data, paste(data.path1,"/lab_wide.baseline.NOhcc.data", sep = ""))

data.pp.to.impute.bin = with(unique(data.pp),
                         data.frame(
                             id        = id,
                             time = as.numeric(time),
                             rna =  as.numeric(rna),
                             cd4 =  as.numeric(cd4),
                             cd8 = as.numeric(cd8) ,
                              AST = as.numeric(AST),  
                               CHOL = as.numeric(CHOL), 
                               hem = as.numeric(hem),
                               pla =  as.numeric(pla),
                               TRIG =  as.numeric(TRIG),
                              ALT =  as.numeric(ALT),
                              CRE =  as.numeric(CRE),
                              GLUC = as.numeric(GLUC),
                              leu = as.numeric(leu) ,
                              BIL = as.numeric(BIL),
                               ALB =  as.numeric(ALB),
                             aghbs = as.numeric(aghbs),
                             antihcv = as.numeric(antihcv),
                             antihbc = as.numeric(antihbc),
                             antihbs =  as.numeric(antihbs)                    
                         ))

data.pp.to.impute = data.pp.to.impute.bin[,c(
    "id",
   "time",
    "rna",
    "cd4",
    "cd8",
    "AST",
    "CHOL",
    "hem",
    "pla",
    "TRIG",
    "ALT",
    "CRE",
    "GLUC",
    "leu",
    "BIL",
    "ALB", "aghbs","antihcv","antihbc","antihbs"
    )
]

data.pp.to.impute = unique(data.pp.to.impute)

test.pp =  do.call(rbind, by(data.pp.to.impute, data.pp.to.impute$id, function(m){apply(m,2, function(x) sum(!is.na(x)))}))
ids.2.reg = rownames(test.pp)[apply(test.pp, 1, function(l) min(l)>=0 )]
warning("not applying the ids.2reg condition")

##voy a crear otro test: que tenga al menos 2 records in one of the liver function values: pla, ALT, BIL, ALB, AST
ids.key.reg = rownames(test.pp)[apply(test.pp, 1, function(l) {sum( (l[c("pla", "ALT","BIL", "ALB", "AST")]) !=0 )>= 2  } )]
print("# with at least 2 of the 5 liver enzimes");print(length(ids.key.reg))
print("# with at least 2 of the 5 liver enzimes with HCC");print(length(intersect(ids.key.reg, ids.ever.hcc)))

data.pp.to.impute = subset(data.pp.to.impute, is.element(id, ids.key.reg))
test.horizontal = apply(data.pp.to.impute, 1, function(x) { sum(!is.na(x[-grep("time|id",names(data.pp.to.impute))])) !=0 })
data.pp.to.impute = data.pp.to.impute[test.horizontal,]

rm(ids.2.reg);rm(test.horizontal);malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));malloc.trim();gc()

data.pp.to.impute.tetris = 
unique(data.frame(data.pp.to.impute %>%
    pivot_longer(cols = names(data.pp.to.impute)[-grep("time|id",names(data.pp.to.impute))], names_to = "variable", values_to = "value")
    ))

data.pp.to.impute.tetris = unique(data.pp.to.impute.tetris[!is.na(data.pp.to.impute.tetris$value),])

data.pp.to.impute.tetris <- data.pp.to.impute.tetris %>%
  group_by(id, time, variable) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

data.pp.to.impute.tetris.wide <- as.data.frame(pivot_wider(data.pp.to.impute.tetris, names_from = variable, values_from = value, values_fill = NA_real_))

data.pp.to.impute.tetris.wide = data.pp.to.impute.tetris.wide[,names(data.pp.to.impute)]

 write.table(data.pp.to.impute.tetris.wide, paste(data.path1,"/data_pp_to_impute_all",sep = ""))
data.pp.to.impute.tetris.wide = read.table(paste(data.path1,"/data_pp_to_impute_all", sep = ""))

rm(data.pp.to.impute.tetris)
rm(data.pp.to.impute)

malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));malloc.trim();gc()

##conditions before imputation, para evitar la fatiga
ii.longer.fup.space =  as.character((by(data.pp.to.impute.tetris.wide, data.pp.to.impute.tetris.wide$id,  function(m) {
    test1 = (max(m$time) - min(m$time))>1.5
    test2 = (max(m$time)>0)

    test = ((test2 + test1) == 2)
    
    if(test) m$id[1]
})))
                              
data.pp.to.impute.tetris.wide = subset(data.pp.to.impute.tetris.wide, is.element(id,ii.longer.fup.space))

print("all right before imputing, enough fup"); print(length(unique((data.pp.to.impute.tetris.wide$id))))
print("all right before imputing, enough fup - HCC");print(length(intersect(unique((data.pp.to.impute.tetris.wide$id)), unique(ids.ever.hcc))))

message("ABOUT TO IMPUTE")

pp.data.imputed = as.data.frame(impute_missing_time_dependent.3(data.pp.to.impute.tetris.wide,
                                                                variable_types = c(rep("continuos", length(names(data.pp.to.impute.tetris.wide))-2 -4), rep("binary",4))))
    
write.table(pp.data.imputed, paste(data.path1,"/pp_data_imputed_allpat", sep = ""))
pp.data.imputed = read.table(paste(data.path1,"/pp_data_imputed_allpat",sep = ""))

pp.data.imputed <- pp.data.imputed %>% group_by(id) %>% arrange(time) %>% fill(aghbs, .direction = "down") %>% ungroup()
pp.data.imputed <- pp.data.imputed %>% group_by(id) %>% arrange(time) %>% fill(antihcv, .direction = "down") %>% ungroup()
pp.data.imputed <- pp.data.imputed %>% group_by(id) %>% arrange(time) %>% fill(antihbc, .direction = "down") %>% ungroup()
pp.data.imputed <- pp.data.imputed %>% group_by(id) %>% arrange(time) %>% fill(antihbs, .direction = "down") %>% ungroup()

pp.data.imputed = as.data.frame(pp.data.imputed)

pp.data.imputed$aghbs[is.na(pp.data.imputed$aghbs)] = 0
pp.data.imputed$antihcv[is.na(pp.data.imputed$antihcv)] = 0
pp.data.imputed$antihbc[is.na(pp.data.imputed$antihbc)] = 0
pp.data.imputed$antihbs[is.na(pp.data.imputed$antihbs)] = 0

write.table(pp.data.imputed, paste(data.path1,"/pp_data_imputed_allpat", sep = ""))
pp.data.imputed = read.table(paste(data.path1,"/pp_data_imputed_allpat",sep = ""))

pp.data.imputed. = pp.data.imputed
data.pp.to.impute. = data.pp.to.impute.tetris.wide
pp.data.imputed.$origin = "imputed"
data.pp.to.impute.$origin = "original"

message("imputation removes ids?")
message("# all after imputation"); print(length(unique(pp.data.imputed$id)))
message("# all after imputation - HCC"); print(length(unique(intersect(pp.data.imputed$id, ids.ever.hcc))))

data.pp.all = rbind(pp.data.imputed., data.pp.to.impute.)

data.pp.all.long = data.pp.all <- data.pp.all %>%
    pivot_longer(cols = names(data.pp.all)[-grep("time|id|origin",names(data.pp.all))], names_to = "variable", values_to = "value")

id.check =  unique(data.pp.all$id)[sample(1:length(unique(data.pp.all$id)),1)]

pp = ggplot(subset(data.pp.all, id ==  id.check), aes(x = time, y = value, col = origin, pch = origin )) +
    facet_wrap(variable~id, scale = "free", ncol = 3) +
    geom_point(alpha = 0.25, size = 3) + theme_classic()
x11(title = "just and imputation case");print(pp)

png("imputation_outcome_example.png", width = 10, height = 10, res = 300, units = "in");print(pp);dev.off()

data.pp.tenofovir.timing = unique(data.pp[,c("id", "time", "before.TFV_TAF")])  
data.pp.tenofovir.timing = do.call(rbind, by(data.pp.tenofovir.timing, data.pp.tenofovir.timing$id, function(m) m[order(m$time),]))

tenofovir.timing = do.call(rbind,by(data.pp.tenofovir.timing, data.pp.tenofovir.timing$id, function(m) {
    m = subset(m, before.TFV_TAF == 0)
if(dim(m)[1]>0) {out = data.frame(id = m$id[1], time.tenofovir.start = min(m$time))}else{out = NULL}
    out}))

pp.data.imputed.tenofovir.start = merge(pp.data.imputed, tenofovir.timing, by = "id", all.x = T)

pp.data.imputed =  subset(pp.data.imputed.tenofovir.start, time.tenofovir.start<=time)[,-grep("time.tenofovir.start", names(pp.data.imputed.tenofovir.start))]

rm(pp.data.imputed.tenofovir.start);rm(tenofovir.timing);rm(data.pp.tenofovir.timing)
malloc.trim <- function() invisible(.C("malloc_trim", as.integer(0)));malloc.trim();gc()

lab.ts.after.tdf. = by(pp.data.imputed, pp.data.imputed$id, function(m) m)


## Reshape data for ggplot (long format)
long_data <- pp.data.imputed %>%
  pivot_longer(cols = names(pp.data.imputed)[-grep("time|id",names(pp.data.imputed))], names_to = "variable", values_to = "value")


# Plot the time series with one panel per variable
x11(title = "before interpolation")
ggplot(long_data, aes(x = time, y = value, color = id)) +
  geom_line() +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(
    title = "Multivariate Time Series for Individuals",
    x = "Time",
    y = "Value",
    color = "Individual"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right"
  ) + theme(legend.position = "none")

  ii.keep = which(unlist(lapply(lab.ts.after.tdf., nrow))>1) #more than one registry

print("before limiting the series to interpolate - tenofovir fup time left");print(length(lab.ts.after.tdf.))
print("before limiting the series to interpolate- tenofovir fup time left - HCC");print(unique(length(intersect(names(lab.ts.after.tdf.), unique(ids.ever.hcc)) )))

lab.ts.after.tdf. =  lab.ts.after.tdf.[ii.keep]

##the the difference between the first and las time is at least 1.5:
ii.keep = which(unlist(lapply(lab.ts.after.tdf., function(m) {
cond1 = ((max(m$time) - min(m$time))>1.5) 
cond2 = (min(m$time)<30)
(cond2 + cond1) == 2
})))

lab.ts.after.tdf. =  lab.ts.after.tdf.[ii.keep]
print("# with 1.5y fup & >=2 points"); print(length(ii.keep))
print("# with 1.5y fup & >=2 points + hcc"); print(length(intersect(names(ii.keep), unique(ids.ever.hcc))))

lab.ts.after.tdf.interpolated = interpolate_to_uniform_time(
    ts_list = lab.ts.after.tdf.,
    variable_types = c(rep("continous", length(names(lab.ts.after.tdf.[[1]])) -2 -4), rep("binary",4)),   step_size = 1/4)

##-check a random interpolation case
ii = sample(1:length(lab.ts.after.tdf.interpolated),1)
inter = as.data.frame(lab.ts.after.tdf.interpolated[[ii]])
datos = as.data.frame(lab.ts.after.tdf.[[ii]])

inter$origin = "inter"
datos$origin = "datos"

inter.datos = rbind(inter,datos)

long_data.inter.datos <- as.data.frame(inter.datos %>%
  pivot_longer(cols = names(inter.datos)[-grep("time|origin|id",names(inter.datos))], names_to = "variable", values_to = "value"))

x11(title = "check a random interpolation");
pp = ggplot(long_data.inter.datos, aes(x = as.numeric(time), y = as.numeric(value), color = as.factor( origin), lty = as.factor(origin) )) +
    geom_point(alpha = 0.5) +
    geom_line() +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(
    title = "Multivariate Time Series for Individuals",
    x = "Time",
    y = "Value",
    color = "Individual"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "right"
  ) + theme(legend.position = "right")

print(pp)

png("interpolation_random_example.png", width = 7, height = 12, res = 300, units = "in");print(pp);dev.off()

lab.ts.after.tdf.interpolated.ids = lab.ts.after.tdf.interpolated
write.table(do.call(rbind, lab.ts.after.tdf.interpolated.ids), "lab_ts_after_tdf_interpolated_ids")
lab.ts.after.tdf.interpolated = lapply(lab.ts.after.tdf.interpolated, function(x) x[,-grep("id",names(x))])
