############################## REDCap ##########################################

add_labels <- function(dataName, varName){
  factor(dataName[, varName],
         levels = attributes(dataName[, varName])$redcapLevels,
         labels = attributes(dataName[, varName])$redcapLabels)
}

process_a1c <- function(cohort, event){  
  
  if (cohort == "p.4T"){
    dat <- dat_p4T %>% select(-a1c_result) %>% rename(a1c_result = a1c_result_percent)
  } else if (cohort == "4T.1"){
    dat <- dat_4T1
  } else if (cohort == "hist"){
    dat <- dat_hist
    dat$redcap_event_name <- event
  } 
  
  dat_a1c <- dat %>%
    mutate(a1c_collection_date = as.Date(ifelse(!is.na(a1c_collection_date), 
                                                as.Date(a1c_collection_date), 
                                                as.Date(a1c_lab_result_date)), 
                                         origin="1970-01-01"),
           a1c_result = str_remove(a1c_result, "%"),   
           a1c_result = str_remove(a1c_result, "[+]"), 
           a1c_result = as.numeric(str_remove(a1c_result, ">"))) %>% 
    
    select(record_id, redcap_event_name, redcap_repeat_instrument, redcap_repeat_instance,
           a1c_collection_date, a1c_result, a1c_collection_method) %>%
    
    arrange(record_id, 
            a1c_collection_date,  
            desc(a1c_result))  
  
  if (cohort != "hist"){
    
    dat_a1c %<>% filter(redcap_event_name %in% event & 
                          redcap_repeat_instrument == "hba1c_tracking" &
                          !is.na(a1c_result)) 
  }
  
  dat_a1c.1 <- dat_a1c %>%
    arrange(record_id, a1c_collection_date, desc(a1c_result)) %>%
    group_by(record_id, a1c_collection_date) %>%
    slice_head() %>%
    ungroup() 
  
  return(list(dat_a1c, dat_a1c.1))
  
} 

############################## LOESS ###########################################

get_optimal_span <- function(x, y){
  
  df <- data.frame(x, y)
  span.seq <- seq(from = 0.25, to = 0.7, by = 0.05) 
  k <- 10        
  set.seed(1)
  folds <- sample(x = 1:k, size = length(x), replace = T)
  cv.error.mtrx <- matrix(rep(x = NA, times = k * length(span.seq)),
                          nrow = length(span.seq), ncol = k)
  
  suppressWarnings(
    for(i in 1:length(span.seq)) {
      for(j in 1:k) {
        loess.fit <- loess(formula = y ~ x, data = df[folds != j, ], span = span.seq[i])
        preds <- predict(object = loess.fit, newdata = df[folds == j, ])
        cv.error.mtrx[i, j] <- mean((df$y[folds == j] - preds)^2, na.rm = T)
      }
    }
  )
  
  cv.errors <- rowMeans(cv.error.mtrx)
  best.span.i <- which.min(cv.errors)
  span.seq[best.span.i]
  
  spannum = span.seq[best.span.i]  
  
  return(spannum)
  
}

plot_loess <- function(cohorts, titleText, col.hist, col.p4T, col.4T1, a1c.target, spanVal){
  
  cex <- 1.3
  
  dat <- dat_long.all %>% filter(cohort %in% cohorts & studyday <= 12*30) %>% 
    arrange(cohort, record_id) %>%
    mutate(cohort = factor(cohort))
  
  dat.b <- dat_base.all %>% filter(cohort %in% cohorts) %>% 
    mutate(cohort = factor(cohort))
  
  with(dat, 
       plot(studymon, cgm_hba1c, cex.lab = cex,
            xlim = c(0, 12), ylim = c(4, 20),  
            xlab = "Months since diabetes onset", ylab = "HbA1c (%)", 
            axes = F, main = titleText, 
            pch = ifelse(cohort == "4T.1", 20, 1), 
            col = rep(c(alpha(col.hist, 0.6),
                        alpha(col.p4T, 0.6),
                        alpha(col.4T1, 0.6)), table(dat$cohort))))
  
  cols <- c(col.hist, col.p4T, col.4T1)
  
  for (i in 1:length(cohorts)){
    dat.i <- subset(dat, cohort == cohorts[i])
    with(dat.i, lines(loess.smooth(studymon, cgm_hba1c, span = spanVal), 
                      lwd = 3, col = cols[i]))
  }
  
  axis(1, at = seq(0, 15, 3), cex.axis = cex)  
  rug(x = 1:12, ticksize = -0.015, side = 1)
  axis(2, at = seq(0, 20, 2), las = 2, cex.axis = cex)
  box()
  
  abline(h = a1c.target, col = "black", lty = 2, lwd = 1)
  text(-0.15, a1c.target + 0.3, paste0(a1c.target, "%"), 
       col = "black", cex = cex, font = 2)
  
  legend("topright", paste0(c("Historical", "Pilot 4T", "4T Study 1"), 
                            " (N=", table(dat.b$cohort), ")"), 
         lwd = 2, cex = cex, col = c(col.hist, col.p4T, col.4T1), bty = 'n')
  
} 

get_loess_means <- function(dat, timepts, spanVal, means = NULL){  

  loess.hist <- with(subset(dat, cohort == "hist"),
                     loess.smooth(I(studyday/30), cgm_hba1c, span = spanVal, 
                                  evaluation = 10*max(timepts)))
  loess.p4T <- with(subset(dat, cohort == "p.4T"),
                    loess.smooth(I(studyday/30), cgm_hba1c, span = spanVal, 
                                 evaluation = 10*max(timepts)))
  loess.4T1 <- with(subset(dat, cohort == "4T.1"),
                    loess.smooth(I(studyday/30), cgm_hba1c, span = spanVal, 
                                 evaluation = 10*max(timepts)))
  
  diff.4T1.minus.hist <- loess.4T1$y - loess.hist$y
  diff.4T1.minus.p4T <- loess.4T1$y - loess.p4T$y
  
  if (is.null(means)){ 
    out <- cbind(diff.4T1.minus.hist, diff.4T1.minus.p4T)[timepts*10,]
  } else { 
    out <- cbind(loess.hist$y, loess.p4T$y, loess.4T1$y)[timepts*10,]
  }
  
  return(out)
  
}

############################## bootstrap #######################################

organize_tab <- function(tab){
  
  diff.tab <- data.frame(do.call(rbind, tab))  
  
  out <- cbind(t(round(sapply((diff.tab %>% filter(substr(row.names(diff.tab), 1, 2) %in% "X6")), 
                              quantile, c(0.025, 0.975)), 2)), 
               t(round(sapply((diff.tab %>% filter(substr(row.names(diff.tab), 1, 2) %in% "X9")), 
                              quantile, c(0.025, 0.975)), 2)), 
               t(round(sapply((diff.tab %>% filter(substr(row.names(diff.tab), 1, 3) %in% "X12")), 
                              quantile, c(0.025, 0.975)), 2)))
  
  return(cbind(paste0("(", sprintf("%.2f", out[,1]), ", ", sprintf("%.2f", out[,2]), ")"),
               paste0("(", sprintf("%.2f", out[,3]), ", ", sprintf("%.2f", out[,4]), ")"),
               paste0("(", sprintf("%.2f", out[,5]), ", ", sprintf("%.2f", out[,6]), ")"))) 
  
}

combine_output <- function(tab.mean, tab.ci){
  
  out <- cbind(paste(tab.mean[,1], tab.ci[,1]),
               paste(tab.mean[,2], tab.ci[,2]))
  
  rownames(out) <- rownames(tab.ci)
  colnames(out) <- colnames(tab.ci)
  
  return(out)
  
}

bootstrap_ci <- function(dat.in, timepts = seq(6, 12, 3), n.times = 1000, seed = 2023){
  
  dat.nest <- dat.in %>% select(record_id, studyday, cgm_hba1c, cohort) %>% 
    nest(data = -record_id)
  
  set.seed(seed)
  dat.boot <- rsample::bootstraps(dat.nest, times = n.times) 
  
  diff.tab.cohort <- map(dat.boot$splits, ~as_tibble(.) %>% unnest() %>%
                         get_loess_means(., timepts = timepts, spanVal = span.mse))
  
  for (i in 1:n.times){
    rownames(diff.tab.cohort[[i]]) <- timepts
  }
  
  out <- organize_tab(tab = diff.tab.cohort)

  rownames(out) <- c("4T Study 1 minus Historical", "4T Study 1 minus Pilot 4T") 
  colnames(out) <- paste("Month", timepts)
  
  return(out)
  
}

############################## CGM metrics #####################################

plot_loess_ci <- function(xvar, yvar, spanval, addRef = NULL, maintxt, ylabtxt, yrange, yseq, 
                          addPilot = NULL, xvar.p4T = NULL, yvar.p4T = NULL, leg.pos = "bottomright"){
  
  col.p4T <- "#3C5488"
  col.4T1 <- "#D44D44" 
  
  plot(xvar, yvar, xlim = c(0, 64), ylim = yrange, 
       main = maintxt, col = alpha(col.4T1, 0.3),
       xlab = "Study Week", ylab = ylabtxt, axes = F)
  axis(1, seq(0, 65, 5))
  axis(2, yseq, las = 2)
  
  if (!is.null(addPilot)){
    
    points(xvar.p4T, yvar.p4T, col = alpha(col.p4T, 0.3))
    
    fit.loess <- predict(loess(yvar.p4T ~ xvar.p4T, span = spanval), se = T) 
    x <- xvar.p4T[!is.na(xvar.p4T) & !is.na(yvar.p4T)]
    fit <- fit.loess$fit
    CI.l <- fit.loess$fit - qt(0.975, fit.loess$df) * fit.loess$se
    CI.u <- fit.loess$fit + qt(0.975, fit.loess$df) * fit.loess$se
    
    dat.loess <- data.frame(cbind(x, fit, CI.l, CI.u))
    dat.loess <- dat.loess[order(dat.loess$x),]
    
    polygon(c(dat.loess$x, rev(dat.loess$x)),
            c(dat.loess$CI.l, rev(dat.loess$CI.u)), 
            col = alpha(col.p4T, 0.5), border = NA)
    
    lines(dat.loess$x, dat.loess$fit, col = "gray30", lwd = 2)
    
    legend(leg.pos, fill = c(col.4T1, col.p4T),
           c("4T Study 1", "Pilot 4T"), bty = 'n')
    
  }
  
  fit.loess <- predict(loess(yvar ~ xvar, span = spanval), se = T) 
  
  x <- xvar[!is.na(xvar) & !is.na(yvar)]
  fit <- fit.loess$fit
  CI.l <- fit.loess$fit - qt(0.975, fit.loess$df) * fit.loess$se
  CI.u <- fit.loess$fit + qt(0.975, fit.loess$df) * fit.loess$se
  
  dat.loess <- data.frame(cbind(x, fit, CI.l, CI.u))
  dat.loess <- dat.loess[order(dat.loess$x),]
  
  polygon(c(dat.loess$x, rev(dat.loess$x)),
          c(dat.loess$CI.l, rev(dat.loess$CI.u)), col = alpha(col.4T1, 0.5), border = NA)
  
  lines(dat.loess$x, dat.loess$fit, col = "gray30", lwd = 2)
  
  if (!is.null(addRef)){
    
    if (addRef %in% c(140, 180)){  
      
      abline(h = 70, lty = 2)
      abline(h = addRef, lty = 2)
      
      text(69, 70, "70", xpd = T, cex = 0.9)   
      text(69, addRef, as.character(addRef), xpd = T, cex = 0.9)
      
    } else {  
      
      abline(h = addRef, lty = 2)
      text(69, addRef, paste0(addRef, "%"), xpd = T, cex = 0.9) 
  
    }
  }
}

plot_agp <- function(ID.list, mainText, SG.n, SG.perc){
  
  thermColors <- c("#a11c10", "#ce231a", "#73b456", "#ffed00", "#fbb400")
  
  xnames <- paste("Month", c(1, 3, 6, 9, 12))
  for (i in 1:length(xnames)){
    npat.i <- SG.n[i]
    xnames[i] <- paste0(xnames[i], "\n(N=", npat.i, ")")
  }
  
  barplot(SG.perc, col = thermColors, main = mainText,
          xlim = c(0, 9), ylim = c(0, 100), 
          ylab = "% glucose TIR",beside = F, border = "white",  
          names.arg = xnames, cex.names = 1.05,
          legend.text = T, 
          args.legend = list(x = "topright", 
                             c("Very Low", "Low", "Target", "High", "Very High"),
                             horiz = F, cex = 0.9),
          space = rep(0.5, 5),
          axes = F)
  axis(2, seq(0, 100, 10), las = 2)
  
  percs <- round(SG.perc[3,])
  for (i in 1:length(percs)){
    text(i*1.5-0.5, percs[i]-4, paste0(percs[i], "%"), cex = 0.8)
  }
  
} 

############################## ADA targets ########################################

get_counts <- function(dat, target){
  
  tab <- dat %>%
    filter(!is.na(cgm_hba1c)) %>%
    
    filter(.data[[target]] %in% 1) %>%
    distinct(record_id, month, .keep_all = T) %>%
    group_by(cohort, month) %>%
    summarise(n = n()) %>% cbind(
      
      (dat %>%
         filter(!is.na(cgm_hba1c)) %>%
         distinct(record_id, month, .keep_all = T) %>%
         count(cohort, month) %>% rename(N = n) %>% select(N))
    ) %>%
    mutate(perc = round(n/N, 3))
  
  out <- cbind(tab[1:(nrow(tab)/3), c("month", "n", "N", "perc")],
               tab[(nrow(tab)/3+1):(nrow(tab)*2/3), c("n", "N", "perc")],
               tab[(nrow(tab)*2/3+1):nrow(tab), c("n", "N", "perc")])
  colnames(out) <- c("month", "n.hist", "N.hist", "p.hist",
                     "n.p4T", "N.p4T", "p.p4T",
                     "n.4T1", "N.4T1", "p.4T1")
  
  return(out)
  
}

create_barplot <- function(dat, thres, titletxt){
  
  tab <- get_counts(dat, thres)
  rownames(tab) <- tab$month 
  
  colors <- c("gray70", "cornflowerblue", "cornflowerblue")
  legText <- c("Historical", "Pilot 4T", "4T Study 1")
  
  barplot(t(100*tab[, c("p.hist", "p.p4T", "p.4T1")]), beside = T, main = titletxt, 
          ylim = c(-10, 100), ylab = "Percentage",
          col = colors, xlab = "Months since T1D onset",
          legend = legText, 
          args.legend = list(y = 105, bty = "n"),
          density = rep(c(NA, 30, NA), 5), 
          angle = rep(c(NA, 45, NA), 5), 
          las = 1)
  
  text(0.5, -5, "N", cex = 0.9, font = 3)
  
  text(seq(1.5, 21.5, 4), 100*tab[,"p.hist"]+4, paste0(round(100*tab[,"p.hist"]), "%"), cex = 0.78)  
  text(seq(2.5, 22.5, 4), 100*tab[,"p.p4T"]+4, paste0(round(100*tab[,"p.p4T"]), "%"), cex = 0.78)  
  text(seq(3.5, 23.5, 4), 100*tab[,"p.4T1"]+4, paste0(round(100*tab[,"p.4T1"]), "%"), cex = 0.78) 
  
  text(seq(1.5, 21.5, 4), -5, tab[,"N.hist"], cex = 0.8)  
  text(seq(2.5, 22.5, 4), -5, tab[,"N.p4T"], cex = 0.8)  
  text(seq(3.5, 23.5, 4), -5, tab[,"N.4T1"], cex = 0.8)  
  
} 

create_barplot.gmi <- function(dat, thres, titletxt, legtext = c("Pilot 4T", "4T Study 1")){

  tab <- dat %>% pull(all_of(thres)) %>% matrix(ncol = 2) %>%
    set_colnames(c("p.p4T", "p.4T1")) %>%
    set_rownames(unique(dat$month)) 
  
  colors <- rep("cornflowerblue", 2)
  
  barplot(t(tab[, c("p.p4T", "p.4T1")]), beside = T, main = titletxt, 
          ylim = c(-10, 100), ylab = "Percentage",
          col = colors, xlab = "Months since T1D onset",
          legend = legtext, 
          args.legend = list(y = 105, bty = "n"),
          density = rep(c(30, NA), 5), 
          angle = rep(c(45, NA), 5), 
          las = 1)
  
  text(0.6, -5, "N", cex = 0.9, font = 3)
  
  text(seq(1.5, 13.5, 3), tab[,"p.p4T"]+4, paste0(round(tab[,"p.p4T"]), "%"), cex = 0.78)  
  text(seq(2.5, 14.5, 3), tab[,"p.4T1"]+4, paste0(round(tab[,"p.4T1"]), "%"), cex = 0.78)  
  
  text(seq(1.5, 13.5, 3), -5, (dat %>% filter(cohort == "p.4T") %>% pull(totn)), cex = 0.8)
  text(seq(2.5, 14.5, 3), -5, (dat %>% filter(cohort == "4T.1") %>% pull(totn)), cex = 0.8)
  
} 

############################## regression ######################################

get_output <- function(fit){  

  tab <- summary(fit)$tTable 
  beta <- tab[,"Value"]; bse <- tab[,"Std.Error"]
  pval <- tab[,"p-value"]
  CI <- cbind(beta - qnorm(0.975)*bse, beta + qnorm(0.975)*bse)
  out <- round(cbind(beta, CI, pval), 3)
  colnames(out) <- c("Estimate", "95%L", "95%U", "P")
  return(out)
  
}   

betasum <- function(index, coef, var){
  
  beta <- sum(coef[index])
  var  <- var[index,index]   
  bse  <- as.double(sqrt(rep(1,length(index))%*%var%*%rep(1,length(index))))
  ci   <- beta+c(-1,1)*1.96*bse
  p    <- 2*(1-pnorm(abs(beta/bse)))
  
  out  <- round(c(beta,bse,ci,p), 3)
  names(out) <- c("Beta", "SE", "CI.L", "CI.U", "P")
  return(out)
  
}   

get_slopes <- function(fit, groupBy){    

  ind_0 <- which(names(fixef(fit)) == "I((studymon - time) * studymon.gt4)")
  ind_1 <- which(names(fixef(fit)) == "grouphist:I((studymon - time) * studymon.gt4)")
  ind_2 <- which(names(fixef(fit)) == "groupp.4T:I((studymon - time) * studymon.gt4)")
  
  tab.slope <- rbind( 
    betasum(c(ind_0, ind_1), fixef(fit), vcov(fit)), 
    betasum(c(ind_0, ind_2), fixef(fit), vcov(fit)),  
    betasum(c(ind_0), fixef(fit), vcov(fit)))        
  
  tab.slope <- cbind(tab.slope, c(get_output(fit)[ind_1, "P"], get_output(fit)[ind_2, "P"], NA))
  
  rownames(tab.slope) <- c("Historical", "Pilot 4T", "4T Study 1")
  colnames(tab.slope)[6] <- "P_interaction"
  out <- tab.slope[,-2]
  
  return(out)
  
}

get_diffs <- function(fit){
  
  ind.grp.hist = which(names(fixef(fit)) == "grouphist") 
  ind.grp.p4T = which(names(fixef(fit)) == "groupp.4T") 
  
  ind.x.hist = which(names(fixef(fit)) == "grouphist:I((studymon - time) * studymon.gt4)")
  ind.x.p4T = which(names(fixef(fit)) == "groupp.4T:I((studymon - time) * studymon.gt4)")
  
  tab.hist <- rbind(
    
    betasum(ind.grp.hist, fixef(fit), vcov(fit)),               
    betasum(c(ind.grp.hist, ind.x.hist), fixef(fit), vcov(fit)),        
    betasum(c(ind.grp.hist, rep(ind.x.hist, 2)), fixef(fit), vcov(fit)),
    betasum(c(ind.grp.hist, rep(ind.x.hist, 3)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.hist, rep(ind.x.hist, 4)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.hist, rep(ind.x.hist, 5)), fixef(fit), vcov(fit)),  
    betasum(c(ind.grp.hist, rep(ind.x.hist, 6)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.hist, rep(ind.x.hist, 7)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.hist, rep(ind.x.hist, 8)), fixef(fit), vcov(fit))) 
  
  tab.p4T <- rbind(
    
    betasum(ind.grp.p4T, fixef(fit), vcov(fit)),               
    betasum(c(ind.grp.p4T, ind.x.p4T), fixef(fit), vcov(fit)),        
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 2)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 3)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 4)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 5)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 6)), fixef(fit), vcov(fit)),  
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 7)), fixef(fit), vcov(fit)), 
    betasum(c(ind.grp.p4T, rep(ind.x.p4T, 8)), fixef(fit), vcov(fit))) 
  
  out <- rbind(rep(NA, 5), tab.hist, rep(NA, 5), tab.p4T)
  rownames(out) <- c("Historical minus 4T1", paste("Month", 4:12), 
                     "Pilot 4T minus 4T1", paste("Month", 4:12))
  
  return(out)
  
}

fit_models <- function(dat.mod = dat_long.all %>% mutate(cohort = factor(cohort)), 
                       n.months = 12, t0 = 4, groupBy = "cohort"){  

  dat.mod <- subset(dat.mod, studyday <= n.months*30 & studyday >= 0 & !is.na(cgm_hba1c)) 
  
  dat.mod$time <- t0
  dat.mod$studymon.gt4 <- with(dat.mod, as.integer(studymon > time))
  dat.mod$group <- factor(dat.mod[,groupBy],  
                          levels = c("4T.1", "hist", "p.4T"))
  
  dat.mod.censor <- subset(dat.mod, is.na(AdvHybridCL.date1) |
                             !is.na(AdvHybridCL.date1) & studyday < AdvHybridCL.date1)
  
  mod.0 <- as.formula("cgm_hba1c ~ group + studymon + 
                       I((studymon-time)*studymon.gt4) +
                       group * I((studymon-time)*studymon.gt4)")
  
  mod.1 <- as.formula("cgm_hba1c ~ group + studymon +                       
                       I((studymon-time)*studymon.gt4) +
                       group * I((studymon-time)*studymon.gt4) +
                       age.at.onset + cgm_sex + Hispanic + Public")
  
  fit0.lmm <- nlme::lme(fixed = mod.0,  
                        random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                        correlation = corCAR1(form = ~ 1 | record_id), data = dat.mod,
                        method = "REML")
  
  fit1.lmm <- nlme::lme(fixed = mod.1,
                        random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                        correlation = corCAR1(form = ~ 1 | record_id), data = dat.mod,
                        method = "REML", na.action = na.exclude)
  
  fit0.lmm.sens1 <- nlme::lme(fixed = mod.0,  
                              random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                              correlation = corCAR1(form = ~ 1 | record_id), data = dat.mod.censor,
                              method = "REML")
  
  fit1.lmm.sens1 <- nlme::lme(fixed = mod.1,
                              random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                              correlation = corCAR1(form = ~ 1 | record_id),
                              na.action = na.exclude, data = dat.mod.censor,
                              method = "REML")
  
  fit0.lmm.sens2 <- nlme::lme(fixed = mod.0,  
                              random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                              correlation = corCAR1(form = ~ 1 | record_id), 
                              data = subset(dat.mod, OpenLoopEver == 1),
                              method = "REML")
  
  fit1.lmm.sens2 <- nlme::lme(fixed = mod.1,
                              random = reStruct( ~ 1 + studymon | record_id, pdClass = "pdSymm", REML = T),
                              correlation = corCAR1(form = ~ 1 | record_id),
                              data = subset(dat.mod, OpenLoopEver == 1),
                              method = "REML", na.action = na.exclude)
  
  est.fit0 <- get_slopes(fit = fit0.lmm, groupBy)
  est.fit0.sens1 <- get_slopes(fit = fit0.lmm.sens1, groupBy)
  est.fit0.sens2 <- get_slopes(fit = fit0.lmm.sens2, groupBy)
  
  est.fit1 <- get_slopes(fit = fit1.lmm, groupBy)
  est.fit1.sens1 <- get_slopes(fit = fit1.lmm.sens1, groupBy)
  est.fit1.sens2 <- get_slopes(fit = fit1.lmm.sens2, groupBy)
  
  tab.slope <- cbind(rbind(rep(NA, 5), est.fit0, 
                           rep(NA, 5), est.fit0.sens1, 
                           rep(NA, 5), est.fit0.sens2), 
                     
                     rbind(rep(NA, 5), est.fit1, 
                           rep(NA, 5), est.fit1.sens1, 
                           rep(NA, 5), est.fit1.sens2))
  
  rownames(tab.slope)[c(2:4, 6:8, 10:12)] <- paste0("  ", rownames(tab.slope)[c(2:4, 6:8, 10:12)])
  rownames(tab.slope)[seq(1, nrow(tab.slope), 4)] <- c("Main", "Sensitivity #1", "Sensitivity #2")
  
  tab.diff <- cbind(get_diffs(fit = fit0.lmm), get_diffs(fit = fit1.lmm))
  out <- list(tab.slope, tab.diff)
  
  return(out)
  
}

create_forestplot <- function(tab.fit, maxCI, titletxt, ref){   

  k <- ifelse(titletxt == "Adjusted", 6, 1)
  beta <- sprintf("%.2f", round(tab.fit[,k],2))
  
  CI.l <- sprintf("%.2f", round(tab.fit[,k+1],2))
  CI.u <- sprintf("%.2f", round(tab.fit[,k+2],2))
  
  tab <- cbind(rownames(tab.fit), paste0(beta, " (", CI.l, ", ", CI.u, ")"))
  tab[tab %in% c("NA (NA, NA)")] <- ""
  
  b.name <- ifelse(ref == 0, "Slope", "OR")
  forestplot(labeltext = rbind(c("Cohort", paste0(b.name, " (95% CI)")), 
                               tab), 
             is.summary = c(T, rep(c(T, F, F, F), 3)),  
             mean = c(NA, tab.fit[, k]),
             lower = c(NA, tab.fit[, k+1]),
             upper = c(NA, tab.fit[, k+2]),
             xlab = "Change in HbA1c since Month 4", boxsize = 0.3, lwd.ci = 1.5, 
             xticks = seq(0, maxCI, 0.2),
             txt_gp = fpTxtGp(label = gpar(cex = 1),   
                              ticks = gpar(cex = 1),   
                              xlab = gpar(cex = 1)),    
             col = fpColors(box = "black", line = "black", zero = "blue"), new_page = T,
             grid = structure(seq(0, maxCI, 0.2), gp = gpar(lty = 3, col = "gray")),   
             clip = c(0.6, maxCI), 
             graphwidth = unit(100, 'mm'),
             zero = ref,
             title = titletxt)
  
}

############################## other ###########################################

get_1yr_metrics <- function(cohort, week){
  
  if (cohort == "p.4T"){
    dat.1yr <- dat_cgm.p4T %>% filter(Week.of.Visit %in% week) 
  } else if (cohort == "4T.1"){
    dat.1yr <- dat_cgm.4T1 %>% filter(Week.of.Visit %in% week) 
  }
  
  out <- cbind(sapply((dat.1yr %>% select(GMI:hypo.severe)), mean, na.rm = T),
               sapply((dat.1yr %>% select(GMI:hypo.severe)), sd, na.rm = T))
  colnames(out) <- c("mean", "sd")
  
  return(round(out, 3))
  
}
