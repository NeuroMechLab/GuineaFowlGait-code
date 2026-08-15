# 18_sfig_pc_steadiness.R — Supplement: the CoM-dynamics PC space coloured by
# steadiness. Principal components of the six standardized CoM-dynamics descriptors
# (recovery, congruity, mechanical cost of transport, collision angle, touchdown energy,
# duty factor) across all classified
# steps. The feature set is the one Table S1 reports (09_gait_continuity.R), the six gated
# descriptors: the hodograph
# signed area is not among the inputs, because it is the measure tested against this space
# in 12_hodograph_validation.R. This panel and Table S1 must quote the same variance shares,
# so the feature list is kept identical to that script's `feats`.
# PC1 is the gait continuum (pendular <-> bouncing);
# PC2 is a mechanical-cost / unsteadiness axis, largely orthogonal to gait. Colouring
# the steps by their stride steadiness (steady / accelerating / decelerating) shows
# that steadiness separates along PC2: steady steps sit low, unsteady steps high;
# within each gait, accelerating steps sit high on both PC1 and PC2 while decelerating
# steps sit higher on PC2 and lower on PC1.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
if (!exists("theme_daley")) source("R/lib_theme.R")
feats <- SAMPLE_DESCRIPTORS   # the gated descriptor set (04_analysis_sample.R)
d <- step %>% filter(as.character(gaitObjective) %in% GAIT_LEVELS,
                     if_all(all_of(feats), is.finite), !is.na(accClass))
Z <- scale(as.matrix(d[, feats])); pc <- prcomp(Z)
v <- pc$sdev^2 / sum(pc$sdev^2)
# orient PC1 with the pendular (walk) end positive; PC2 with unsteady positive
s1 <- if (mean(pc$x[d$gaitObjective=="walk", 1]) < 0) -1 else 1
s2 <- if (mean(pc$x[d$accClass=="steady", 2]) > mean(pc$x[d$accClass!="steady", 2])) -1 else 1
d$PC1 <- pc$x[,1]*s1; d$PC2 <- pc$x[,2]*s2
L1 <- pc$rotation[,1]*s1; L2 <- pc$rotation[,2]*s2

cat(sprintf("18_sfig_pc_steadiness: PC1=%.1f%% PC2=%.1f%%\n", 100*v[1], 100*v[2]))
cat("PC1 loadings:", paste(sprintf("%s=%+.2f", feats, L1), collapse="  "), "\n")
cat("PC2 loadings:", paste(sprintf("%s=%+.2f", feats, L2), collapse="  "), "\n")
cen <- d %>% group_by(gaitObjective, accClass) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = dplyr::n(), .groups = "drop")
cat("\nWithin-gait PC means by steadiness (accel high on both; decel higher PC2, lower PC1):\n")
print(as.data.frame(cen %>% mutate(PC1 = round(PC1,2), PC2 = round(PC2,2))))

pal <- c(decelerating = "#0072B2", steady = "#6c757d", accelerating = "#E69F00")
glab <- c(walk = "walk", groundedRun = "grounded run", aerialRun = "aerial run")
d$gplot <- factor(glab[as.character(d$gaitObjective)], levels = glab)
cen$gplot <- factor(glab[as.character(cen$gaitObjective)], levels = glab)
p <- ggplot(d, aes(PC1, PC2, colour = accClass)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_point(data = cen, aes(PC1, PC2, fill = accClass), shape = 21, size = 3.2,
             colour = "black", stroke = 0.6, inherit.aes = FALSE) +
  scale_colour_manual(values = pal, name = "stride steadiness") +
  scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~gplot) +
  labs(x = sprintf("PC1 (%.0f%% of variance)", 100*v[1]),
       y = sprintf("PC2 (%.0f%% of variance)", 100*v[2])) +
  theme_daley()
ggsave("output/SFig_PCSteadiness.pdf", p, width = 10, height = 4, device = cairo_pdf)
ggsave("output/SFig_PCSteadiness.png", p, width = 10, height = 4, dpi = 200, bg = "white")
cat("18_sfig_pc_steadiness: wrote SFig_PCSteadiness.\n")
