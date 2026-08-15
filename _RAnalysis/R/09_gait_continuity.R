# 09_gait_continuity.R — is the CoM-dynamics feature space a continuum or a set of
# discrete clusters? Two independent, unsupervised summaries of the SAME standardized
# feature matrix Z (pendular recovery, KE-PE congruity, mechanical cost of transport):
#   (1) PAM clustering scored by average silhouette width over k = 2..6 (Table S3), which
#       tests for separable groups, plus the agreement between the k = 2 partition and
#       hodograph rotation sense;
#   (2) PCA on the same Z (Table S1), which tests for low-dimensional structure.
# PAM operates on Z directly; the PCA is NOT fed into PAM. A low, flat silhouette profile
# plus a leading gait-continuum component (PC1) indicates a continuum rather than discrete
# clusters. PC2 is a second, largely gait-independent axis (mechanical cost / unsteadiness;
# see the PC-space supplement, Fig S3).
#
# The hodograph signed area is NOT one of the clustered features, so the agreement below is a
# genuine test of an unsupervised split against an independent measure.
#
# The agreement with the CLASSIC criterion (congruity at 50%, the classifier defined in
# 02_clean.R) is also reported, but as internal consistency rather than validation, because
# congruity IS one of the clustered features.
#
# This script emits tables only. The clustering evidence is Table S3 and the agreement statistics
# below. The PCA here is also what 10_fig_gaitspace_planes.R draws as Fig 2; that script recomputes
# it from the same inputs and stops if the variance shares disagree with pca_summary.csv.
#
# Outputs
#   output/pca_summary.csv            loadings on all six PCs + variance and running total (Table S1)
#   output/pca_loadings_summary.csv   leading-two-PC form (retained for earlier references)
#   output/pam_clustering_summary.csv silhouette and cluster sizes by k (Table S3)
#   output/pam_rotation_agreement.csv k = 2 partition vs rotation sense: 2x2 table + statistics
#   output_internal/pam_classic_agreement.csv  k = 2 partition vs the classic criterion (internal
#                                     consistency only; congruity is a clustered feature)
#   output/gaitspace_stats.csv        one-line headline summary
suppressPackageStartupMessages({library(dplyr)})
if (!requireNamespace("cluster", quietly = TRUE))
  stop("09_gait_continuity: the 'cluster' package is required for the PAM analysis.")
if (!dir.exists("output")) dir.create("output")

# The six established continuous descriptors of CoM dynamics. Every one is defined for every step
# in the analysis sample, because descriptor completeness is part of the gate (04_analysis_sample.R),
# so this analysis runs on exactly the sample every other analysis uses and needs no filter of its
# own. `feats` is deliberately the same object 04_analysis_sample.R gates on.
#
# A dimensionality claim over six descriptors is only as strong as their independence, so the
# correlation matrix is emitted alongside (output/descriptor_correlations.csv). Three pairs are
# strongly related, and one of them is related by definition: the force- and velocity-weighted
# collision angle and the dimensionless mechanical cost of transport are two measures of a single
# theoretical quantity (Lee et al. 2011, 2013), so their loading together is expected and is not
# independent corroboration of an axis.
feats <- SAMPLE_DESCRIPTORS
# One source for the descriptor names: SAMPLE_DESCRIPTOR_LABELS in 04_analysis_sample.R, so the
# PCA table and the correlation matrix cannot drift apart in what they call the same descriptor.
featlab <- SAMPLE_DESCRIPTOR_LABELS
# no descriptor filter here: 04_analysis_sample.R already guarantees all of `feats` are finite
d <- step %>%
  filter(is.finite(hodoArea), is.finite(meanSpeed_n)) %>%
  mutate(rot = factor(ifelse(hodoArea > 0, "pendular (walk)", "bouncing (run)"),
                      levels = c("pendular (walk)", "bouncing (run)")))

Z <- scale(as.matrix(d[, feats]))

# ---- (1) PAM: cluster tendency over k = 2..6 -------------------------------------------
KS <- 2:6
pam_fits <- lapply(KS, function(k) cluster::pam(Z, k))
sil <- vapply(pam_fits, function(f) f$silinfo$avg.width, numeric(1))
pam_tbl <- data.frame(
  k = KS,
  # fixed decimal places, so a column does not mix 0.280 with 0.271 as "0.28" and "0.271"
  avg_silhouette = sprintf("%.3f", sil),
  cluster_sizes = vapply(pam_fits, function(f)
    paste(sort(as.integer(table(f$clustering)), decreasing = TRUE), collapse = " / "), character(1)))
readr::write_csv(pam_tbl, "output/pam_clustering_summary.csv")

# ---- Agreement between the k = 2 partition and hodograph rotation sense ----------------
# Both labels are binary, so agreement is summarized three ways: raw percent agreement,
# Cohen's kappa (chance-corrected, with a large-sample 95% CI half-width), and the phi
# coefficient, which for two binary variables IS the Pearson correlation between them.
# Cluster labels are arbitrary, so the cluster dominated by counterclockwise steps is
# aligned with the pendular sense before the table is formed.
cl2  <- pam_fits[[which(KS == 2)]]$clustering
pend <- d$hodoArea > 0                       # hodograph rotation sense: the independent measure
pend_classic <- d$congruity < 50             # the classifier of 02_clean.R (a clustered feature)
pend_cluster <- which.max(tapply(pend, cl2, mean))
clus_pend <- cl2 == as.integer(names(pend_cluster))

# Both labels are binary. agreement() returns the 2x2 counts and the three summaries, so the
# partition can be scored against the hodograph (a genuine test: hodoArea is not in Z) and
# against the classic criterion (internal consistency only: congruity IS in Z).
agreement <- function(x, y) {
  a  <- sum(x & y); b <- sum(x & !y); cc <- sum(!x & y); dd <- sum(!x & !y)
  n  <- a + b + cc + dd
  po <- (a + dd) / n
  p_row <- c(a + b, cc + dd) / n
  p_col <- c(a + cc, b + dd) / n
  pe <- sum(p_row * p_col)
  kappa <- (po - pe) / (1 - pe)
  # Fleiss, Cohen & Everitt (1969) large-sample variance of kappa. The off-diagonal term is
  # sum_{i != j} p_ij (p_.i + p_j.)^2, so the b cell (row 1, col 2) pairs p_.1 with p_2. and
  # the c cell (row 2, col 1) pairs p_.2 with p_1.. Checked against a 4000-draw nonparametric
  # bootstrap over the paired labels, which agreed; the emitted half-width is the analytic one.
  var_k <- ((a/n) * (1 - (p_row[1] + p_col[1]) * (1 - kappa))^2 +
            (dd/n) * (1 - (p_row[2] + p_col[2]) * (1 - kappa))^2 +
            (1 - kappa)^2 * ((b/n) * (p_col[1] + p_row[2])^2 + (cc/n) * (p_col[2] + p_row[1])^2) -
            (kappa - pe * (1 - kappa))^2) / (n * (1 - pe)^2)
  phi <- (a * dd - b * cc) / sqrt(as.numeric(a + b) * (cc + dd) * (a + cc) * (b + dd))
  # Wilson 95% CI on the percent agreement
  zz <- 1.959964; den <- 1 + zz^2 / n
  po_lo <- ((po + zz^2/(2*n)) - zz * sqrt(po*(1-po)/n + zz^2/(4*n^2))) / den
  po_hi <- ((po + zz^2/(2*n)) + zz * sqrt(po*(1-po)/n + zz^2/(4*n^2))) / den
  list(n = n, a = a, b = b, cc = cc, dd = dd, po = po, po_lo = po_lo, po_hi = po_hi,
       kappa = kappa, kappa_hw = 1.96 * sqrt(var_k), phi = phi, pe = pe)
}

agree_table <- function(g, lab2) data.frame(
  quantity = c("steps compared (n)",
               sprintf("both pendular (pendular cluster, pendular %s)", lab2),
               sprintf("cluster pendular, %s bouncing", lab2),
               sprintf("cluster bouncing, %s pendular", lab2),
               sprintf("both bouncing (bouncing cluster, bouncing %s)", lab2),
               # the total as well as its two parts: the agreeing count is what the prose quotes,
               # and two disjoint counts beside a percentage invite the reader to add them wrongly
               "steps agreeing (n)",
               "percent agreement",
               "percent agreement, Wilson 95% CI lower",
               "percent agreement, Wilson 95% CI upper",
               "Cohen's kappa",
               "Cohen's kappa, 95% CI half-width",
               "phi coefficient (Pearson r of the two binary labels)",
               "chance agreement expected under independence"),
  # counts as integers, everything below them at a fixed number of places, so a column does not
  # print the Wilson upper bound as "90" beside 88.8 and 87.5
  value = c(sprintf("%d", c(g$n, g$a, g$b, g$cc, g$dd, g$a + g$dd)),
            sprintf("%.1f", 100 * c(g$po, g$po_lo, g$po_hi)),
            sprintf("%.3f", c(g$kappa, g$kappa_hw, g$phi, g$pe))))

# The headline test: an unsupervised split that never saw the hodograph, scored against it.
gr <- agreement(clus_pend, pend)

# Interval on the percent agreement. The Wilson interval below assumes independent observations,
# which these are not, since each bird contributes many steps. It is kept for comparability with the
# literature, but the interval the paper reports is a cluster bootstrap over individuals, the
# same resampling unit 12_hodograph_validation.R uses for the crossover speeds. Reporting a
# Wilson interval here while bootstrapping there was an inconsistency: the two analyses made
# opposite assumptions about the same non-independence. Both intervals are emitted to
# pam_rotation_agreement.csv; the bootstrap one is the wider and is what the manuscript quotes.
NBOOT_AGREE <- 2000
set.seed(19)
agr_vec  <- clus_pend == pend
ids_a    <- unique(as.character(d$subjectID))
idx_a    <- split(seq_along(agr_vec), as.character(d$subjectID))
boot_po  <- replicate(NBOOT_AGREE,
  mean(agr_vec[unlist(idx_a[sample(ids_a, length(ids_a), replace = TRUE)], use.names = FALSE)]))
po_ci <- stats::quantile(boot_po, c(0.025, 0.975), na.rm = TRUE)

readr::write_csv(rbind(
    agree_table(gr, "rotation sense"),
    data.frame(quantity = c("percent agreement, cluster bootstrap 95% CI lower",
                            "percent agreement, cluster bootstrap 95% CI upper",
                            "cluster bootstrap replicates (over individuals)"),
               value = c(sprintf("%.1f", 100 * po_ci[1]), sprintf("%.1f", 100 * po_ci[2]),
                         sprintf("%d", NBOOT_AGREE)))),
  "output/pam_rotation_agreement.csv")
cat(sprintf("  agreement %.1f%%: Wilson %.1f to %.1f, cluster bootstrap %.1f to %.1f\n",
            100*gr$po, 100*gr$po_lo, 100*gr$po_hi, 100*po_ci[1], 100*po_ci[2]))

# Internal consistency with the classifier. Reported separately and NOT as validation,
# because congruity is one of the four clustered features.
gc_ <- agreement(clus_pend, pend_classic)
readr::write_csv(agree_table(gc_, "classic criterion"),
                 "output_internal/pam_classic_agreement.csv")

# Agreement within each gait of the classification, which is where disagreement concentrates.
by_gait <- data.frame(gait = as.character(d$gaitObjective), agree = clus_pend == pend) %>%
  group_by(gait) %>%
  summarise(n_steps = dplyr::n(), agree_pct = sprintf("%.1f", 100 * mean(agree)), .groups = "drop")
readr::write_csv(by_gait, "output/pam_rotation_agreement_by_gait.csv")

# names kept for the summary line and the numbers reference
n <- gr$n; po <- gr$po; kappa <- gr$kappa; kappa_hw <- gr$kappa_hw; phi <- gr$phi

# ---- (2) PCA on the same standardized matrix (independent of PAM) ----------------------
# PC sign is arbitrary, so PC1 is oriented with the pendular (walk) end positive by fixing
# pendular recovery > 0, and PC2 with the unsteady end positive.
pv <- prcomp(Z); ve <- 100 * pv$sdev^2 / sum(pv$sdev^2)
L <- pv$rotation; if (L["recovery", 1] < 0) L[, 1] <- -L[, 1]
readr::write_csv(data.frame(feature = unname(featlab[rownames(L)]),
                            PC1_loading = round(L[, 1], 3),
                            PC2_loading = round(L[, 2], 3), row.names = NULL),
                 "output/pca_loadings_summary.csv")

L2 <- L
if ("accClass" %in% names(d)) {
  iS <- which(d$accClass == "steady"); iU <- which(d$accClass != "steady")
  if (length(iS) && length(iU) && mean(pv$x[iS, 2]) > mean(pv$x[iU, 2])) L2[, 2] <- -L2[, 2]
}
# fixed decimal places throughout, so a column does not mix "-0.3" with "-0.393"
pca_full <- data.frame(feature = unname(featlab[rownames(L2)]),
                       apply(L2, 2, function(z) sprintf("%.3f", z)),
                       row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
names(pca_full)[-1] <- paste0("PC", seq_len(ncol(L2)))
pca_full <- rbind(pca_full,
  setNames(c(list("Variance explained (%)"), as.list(sprintf("%.1f", ve))), names(pca_full)),
  setNames(c(list("Cumulative variance (%)"), as.list(sprintf("%.1f", cumsum(ve)))), names(pca_full)))
readr::write_csv(pca_full, "output/pca_summary.csv")


readr::write_csv(data.frame(pca_n = nrow(d),
                            silhouette_max = round(max(sil), 3),
                            silhouette_k = KS[which.max(sil)],
                            k2_rotation_match_pct = round(100 * po, 1),
                            k2_rotation_kappa = round(kappa, 3),
                            k2_rotation_phi = round(phi, 3),
                            k2_classic_match_pct = round(100 * gc_$po, 1),
                            k2_classic_kappa = round(gc_$kappa, 3)),
                 "output/gaitspace_stats.csv")

cat(sprintf(paste0("09_gait_continuity: %d steps, %d hodograph-free features; PCA PC1=%.1f%% ",
                   "PC2=%.1f%% of variance; PAM silhouette max %.3f at k=%d.\n",
                   "  k=2 vs rotation sense (independent of Z): %.1f%% agreement, ",
                   "kappa %.3f +/- %.3f, phi %.3f.\n",
                   "  k=2 vs classic criterion (congruity is IN Z, consistency only): ",
                   "%.1f%% agreement, kappa %.3f.\n"),
            nrow(d), length(feats), ve[1], ve[2], max(sil), KS[which.max(sil)],
            100 * po, kappa, kappa_hw, phi, 100 * gc_$po, gc_$kappa))
print(pam_tbl); print(as.data.frame(by_gait))
