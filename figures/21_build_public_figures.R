# Rebuild the ROC, tertile, balance and robustness figures from local analysis
# inputs. Never export patient-level values or observed ROC thresholds.
args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ra <- Sys.getenv("RA_ROOT", unset = if (length(args)) file.path(dirname(args[1]), "..") else ".")
ra <- normalizePath(ra, winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(ra, "Rlib"), .libPaths()))
tab <- file.path(ra, "outputs", "tables")
int <- file.path(ra, "outputs", "interim")
out <- Sys.getenv("FRONTMED_FIGURE_DIR", unset = file.path(ra, "outputs", "figures"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
which_figures <- strsplit(Sys.getenv("FRONTMED_PLOT_WHICH", unset = "2,3,4,5"), ",", fixed = TRUE)[[1]]
read_tab <- function(name) read.csv(file.path(tab, name), check.names = FALSE)
save_pair <- function(stem, draw, width, height) {
  png(file.path(out, paste0(stem, ".png")), width = width, height = height, res = 300)
  draw(); dev.off()
  tiff(file.path(out, paste0(stem, ".tif")), width = width, height = height,
       res = 300, compression = "lzw")
  draw(); dev.off()
}

if ("2" %in% which_figures) {
  suppressPackageStartupMessages(library(pROC))
  m <- read.csv(file.path(int, "cbc_export_v2.csv"), check.names = TRUE)
  e <- read_tab("eicu_sicu_external_validation_dataset.csv")
  e <- e[e$values_valid & e$outcome_available & e$alc > 0, ]
  # Recompute from source components to avoid CSV-rounded index ties.
  e$neutrophil <- e$wbc_k * e$polys_pct / 100
  e$lymphocyte <- e$wbc_k * e$lymphs_pct / 100
  amc <- e$wbc_k * e$monos_pct / 100
  e$NLR <- e$neutrophil / e$lymphocyte
  e$SII <- e$platelet_k * e$neutrophil / e$lymphocyte
  e$SIRI <- e$neutrophil * amc / e$lymphocyte
  e$PIV <- e$platelet_k * e$neutrophil * amc / e$lymphocyte
  e$PLR <- e$platelet_k / e$lymphocyte
  stopifnot(nrow(m) == 879, sum(m$hospital_expire_flag) == 87,
            nrow(e) == 2855, sum(e$death) == 276)
  vars_m <- c("nlr", "alc_single", "anc_single", "siri", "sii", "piv", "plr")
  # The canonical export names its selected analysis lymphocyte column 'alc'.
  if (!"alc_single" %in% names(m)) vars_m[2] <- "alc"
  vars_e <- c("NLR", "lymphocyte", "neutrophil", "SIRI", "SII", "PIV", "PLR")
  labels <- c("NLR", "Lymphocyte", "Neutrophil", "SIRI", "SII", "PIV", "PLR")
  cols <- c("#17365D", "#C83E4D", "#A7ADB4", "#686E75", "#E28A27", "#7651A6", "#C7CBD0")
  curves <- list(); aucs <- matrix(NA_real_, 7, 2)
  for (panel in 1:2) {
    d <- if (panel == 1) m else e
    yy <- if (panel == 1) d$hospital_expire_flag else d$death
    vv <- if (panel == 1) vars_m else vars_e
    curves[[panel]] <- lapply(seq_along(vv), function(j) {
      if (panel == 2) {
        score <- if (j == 2) log1p(d[[vv[j]]]) else -log1p(d[[vv[j]]])
        r <- roc(yy, score, direction = ">", quiet = TRUE)
      } else r <- roc(yy, d[[vv[j]]], direction = if (j == 2) ">" else "<", quiet = TRUE)
      aucs[j, panel] <<- as.numeric(auc(r)); r
    })
  }
  mt <- read_tab("table2_single_auc.csv")
  lt <- read_tab("lymphocyte_auc_positive_control.csv")
  et <- read_tab("eicu_sicu_single_auc.csv")
  expected_m <- c(mt$AUC[match("NLR", mt$metric)],
    lt$AUC[match("lymphocyte_absolute", lt$metric)],
    lt$AUC[match("neutrophil_absolute", lt$metric)],
    mt$AUC[match(c("SIRI", "SII", "PIV", "PLR"), mt$metric)])
  expected_e <- et$AUC[match(c("NLR", "lymphocyte", "neutrophil", "SIRI", "SII", "PIV", "PLR"), et$metric)]
  cat("Figure 2 MIMIC AUCs:", paste(round(aucs[,1],6),collapse=", "), "\n")
  cat("Figure 2 eICU AUCs:", paste(round(aucs[,2],6),collapse=", "), "\n")
  cat("Archived eICU AUCs:", paste(round(expected_e,6),collapse=", "), "\n")
  stopifnot(all(is.finite(aucs)), max(abs(aucs[,1]-expected_m)) < 1e-10,
            max(abs(aucs[,2]-expected_e)) < 1e-10)
  cat("Figure 2: all 14 AUCs agree with archived aggregate tables.\n")
  draw <- function() {
    layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(4, 1))
    for (panel in 1:2) {
      par(mar = c(4, 4, 3, 1))
      plot(0:1, 0:1, type = "n", xlim = c(0,1), ylim = c(0,1),
           xlab = "1 - Specificity", ylab = "Sensitivity", main = c("A  MIMIC-IV", "B  eICU-CRD")[panel])
      abline(0,1,lty=2,col="grey60")
      for (j in seq_along(labels)) {
        r <- curves[[panel]][[j]]
        lines(1-r$specificities, r$sensitivities, col=cols[j], lwd=if(j==1)2.5 else 1.5)
      }
    }
    par(mar=c(0,0,0,0));plot.new()
    legend("center", sprintf("%s (%.3f / %.3f)",labels,aucs[,1],aucs[,2]),
           col=cols,lwd=2,bty="n",ncol=4,cex=.8)
  }
  save_pair("Figure2_roc_two_panels",draw,3600,1620)
}

if ("3" %in% which_figures) {
  z <- read_tab("platelet_tertile_mortality.csv")
  draw <- function() {
    par(mar=c(5,5,4,1))
    x <- barplot(z$mortality_pct,names.arg=z$platelet_tertile,
      col=c("#17365D","#9AA3AD","#6E7882"),ylim=c(0,max(z$mortality_pct)*1.18),
      ylab="In-hospital mortality (%)",xlab="Preoperative platelet count tertile",
      main=sprintf("In-hospital mortality by platelet tertile (n = %d)",sum(z$n)))
    text(x,z$mortality_pct+.5,sprintf("%.2f%%",z$mortality_pct),font=2)
    text(x,rep(.7,nrow(z)),sprintf("%d/%d",z$events,z$n),col="white")
  }
  save_pair("FigureS1_platelet_tertile",draw,2160,1560)
}

if ("4" %in% which_figures) {
  z <- read_tab("ipw_balance.csv");z <- z[z$comparison=="NLR_vs_SII",]
  z <- z[order(abs(z$smd_unweighted)),]
  draw <- function() {
    par(mar=c(5,12,4,1));y <- seq_len(nrow(z))
    plot(abs(z$smd_unweighted),y,type="n",yaxt="n",xlim=c(0,max(abs(z$smd_unweighted))*1.15),
         xlab="Absolute standardized mean difference",ylab="",main="Covariate balance before and after IPW (NLR vs SII)")
    abline(v=.1,lty=2,col="#C83E4D")
    lines(abs(z$smd_unweighted),y,col="#9AA3AD")
    points(abs(z$smd_unweighted),y,pch=1,col="#6E7882")
    points(abs(z$smd_weighted),y,pch=16,col="#17365D")
    axis(2,at=y,labels=gsub("admission_type: |organ: ","",z$covariate),las=1,cex.axis=.8)
    legend("bottomright",c("Unweighted","After IPW"),pch=c(1,16),col=c("#6E7882","#17365D"),bty="n")
  }
  save_pair("FigureS2_ipw_loveplot",draw,2520,2280)
}

if ("5" %in% which_figures) {
  a <- read_tab("available_case_primary.csv"); b <- read_tab("common_sample_sensitivity.csv")
  cl <- read_tab("cluster_bootstrap_ci.csv"); fs <- read_tab("first_stay_sensitivity.csv")
  em <- read_tab("emergency_stratified.csv"); aa <- read_tab("reviewer_all_absolute_subgroup_ci.csv")
  ip <- read_tab("ipw_weighted_primary.csv"); ib <- read_tab("ipw_bootstrap_ci.csv")
  iy <- read_tab("ipw_outcome_conditioned_primary.csv"); iyb <- read_tab("ipw_outcome_conditioned_bootstrap_ci.csv")
  ei <- read_tab("eicu_sicu_paired_delong_bootstrap.csv")
  as_free <- function(z,j) c(-z$delta_AUC[j],-z$CI95_upper[j],-z$CI95_lower[j])
  labels <- c("Primary available-case","Common five-index sample","Subject-clustered bootstrap", "First stay per subject", "Emergency/urgent stratum", "All-absolute subgroup", "IPW, baseline covariates", "IPW, outcome-conditioned", "eICU directional consistency")
  result <- list()
  for (j in 1:2) {
    # Cluster resampling changes uncertainty, not the common-sample point estimate.
    v <- rbind(as_free(a,j),as_free(b,j),
      c(-b$delta_AUC[j],-cl$cluster_CI95_upper[j],-cl$cluster_CI95_lower[j]),
      as_free(fs,j),as_free(em,j),as_free(aa,j),
      c(ip$weighted_delta_AUC[j],ib$CI95_lower[j],ib$CI95_upper[j]),
      c(iy$weighted_delta_AUC[j],iyb$CI95_lower[j],iyb$CI95_upper[j]),as_free(ei,j))
    colnames(v) <- c("delta","lower","upper")
    result[[j]] <- v
  }
  draw <- function() {
    par(mfrow=c(1,2),mar=c(5,ifelse(1==1,10,2),4,3),oma=c(1,0,2,0))
    for(j in 1:2) {
      par(mar=c(5,if(j==1)10 else 1.5,3,4))
      v <- result[[j]];y <- rev(seq_len(nrow(v)))
      plot(v[,1],y,type="n",xlim=c(-.02,.245),ylim=c(.5,9.5),yaxt="n",
           xlab=expression(Delta*AUC~"(95% CI)"),ylab="",main=c("NLR vs SII","SIRI vs PIV")[j],bty="l")
      abline(v=0,lty=2,col="#C83E4D");abline(v=seq(.05,.2,.05),col="grey93")
      if(j==1)axis(2,at=y,labels=labels,las=1,tick=FALSE,cex.axis=.76)
      segments(v[,2],y,v[,3],y,col="#17365D",lwd=1.5)
      segments(v[,2],y-.1,v[,2],y+.1,col="#17365D")
      segments(v[,3],y-.1,v[,3],y+.1,col="#17365D")
      points(v[,1],y,pch=16,col="#17365D",cex=1)
      text(v[,3]+.005,y,sprintf("%.3f (%.3f, %.3f)",v[,1],v[,2],v[,3]),pos=4,offset=0,cex=.63,xpd=NA)
    }
    mtext(expression("Robustness of the paired "*Delta*AUC*" across analyses"),outer=TRUE,side=3,line=.4)
    mtext("Positive values favour the platelet-free index; eICU is a general surgical ICU direction-consistency cohort.",outer=TRUE,side=1,line=-.2,cex=.7)
  }
  save_pair("Figure3_robustness_forest",draw,3450,1920)
  # Aggregate plotting source, useful for independently checking every label.
  fwrite <- function(x,p) write.csv(x,p,row.names=FALSE)
  fwrite(rbind(data.frame(comparison="NLR vs SII",analysis=labels,result[[1]]),
               data.frame(comparison="SIRI vs PIV",analysis=labels,result[[2]])),
         file.path(out,"Figure3_plotting_values.csv"))
}
cat("Requested figure outputs rebuilt; no row-level values were exported.\n")
