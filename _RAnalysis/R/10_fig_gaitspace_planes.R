# 10_fig_gaitspace_planes.R — Fig 2: the gait space of the six established CoM-dynamics
# descriptors, drawn on two planes and under two colourings.
#
#   A  the two highest-|PC1-loading| descriptors against each other, coloured by dimensionless speed
#   B  PC1 against PC2, coloured by dimensionless speed
#   C  the same two descriptors, coloured by hodograph rotation sense
#   D  PC1 against PC2, coloured by rotation sense
#
# The left column is the descriptor plane a reader already knows how to read; the right column is
# the same steps after the PCA of 09_gait_continuity.R has compressed all six descriptors onto two
# axes. The top row shows the speed ordering along the continuum, the bottom row shows where the
# hodograph rotation sense divides it. The signed loop area is not among the PCA inputs, which is
# what makes the rotation-sense colouring an independent measure against this space rather than a
# variable read back against itself.
#
# The two descriptors on the left are CHOSEN FROM THE DATA, by ranking the PC1 loadings, and the
# ranking goes to output/gaitspace_figure_stats.csv so the choice can be read rather than assumed.
# The higher loading takes the x axis.
#
# The PCA repeats 09_gait_continuity.R: the same six gated descriptors (SAMPLE_DESCRIPTORS), the
# same analysis sample, the same standardization, and the same sign conventions (PC1 positive at
# the pendular end, PC2 positive at the unsteady end). The variance shares the axes print are
# therefore the Table S1 values, and the script stops if they disagree with pca_summary.csv.
#
# Reference lines are drawn only for descriptors carrying a defined division: duty factor 0.5 (the
# aerial-phase boundary) and congruity 50% (the divide the classifier of 02_clean.R uses). A
# different pair selected by the loading rank simply gets no line.
#
# Two display choices, both recorded in the stats file because both change what the reader sees.
# The PC panels are zoomed to PC_XLIM by PC_YLIM, without which a few steps high on PC2 flatten the
# cloud into a band; the number of steps left off the view is emitted. And in C and D the rotation
# classes are drawn minority last, so the pendular steps are not buried under the more numerous
# bouncing ones: draw order shows where each class sits, not how numerous it is.
#
# Outputs
#   output/Fig2_GaitSpace.{pdf,png}
#   output/gaitspace_figure_stats.csv   loading rank, plotted pair, separations, off-view count
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
if (!exists("theme_daley")) source("R/lib_theme.R")

feats <- SAMPLE_DESCRIPTORS
# Axis labels as the other figures write them (11_fig_energy_exchange.R, 16_fig_speed_relations.R).
featlab <- c(recovery = "Pendular energy recovery (%)", congruity = "KE-PE congruity (%)",
             CoTmech = "Mechanical cost of transport", collisionAngle = "Collision angle (deg)",
             E_TD_n = "Touchdown mechanical energy", dutyFactor = "Duty factor (per limb)")
# The PCA takes the collision angle in radians, the form the cost-of-transport equivalence needs;
# degrees are the readable unit, so a panel plots collisionAngle_deg (03_step_descriptors.R carries
# both). Every other descriptor plots the column it is clustered on.
PLOTCOL <- c(collisionAngle = "collisionAngle_deg")
plotcol <- function(v) if (v %in% names(PLOTCOL)) PLOTCOL[[v]] else v
REFLINE <- c(dutyFactor = 0.5, congruity = 50)
PC_XLIM <- c(-6, 6); PC_YLIM <- c(-3, 6)

d <- step %>%
  filter(is.finite(hodoArea), is.finite(meanSpeed_n)) %>%
  mutate(rot = factor(ifelse(hodoArea > 0, "pendular (walk)", "bouncing (run)"),
                      levels = c("pendular (walk)", "bouncing (run)")))

# ---- PCA, identical to 09_gait_continuity.R ---------------------------------------------
Z <- scale(as.matrix(d[, feats]))
pv <- prcomp(Z)
ve <- 100 * pv$sdev^2 / sum(pv$sdev^2)
L <- pv$rotation
if (L["recovery", 1] < 0) { L[, 1] <- -L[, 1]; pv$x[, 1] <- -pv$x[, 1] }
iS <- which(d$accClass == "steady"); iU <- which(d$accClass != "steady")
if (length(iS) && length(iU) && mean(pv$x[iS, 2]) > mean(pv$x[iU, 2])) {
  L[, 2] <- -L[, 2]; pv$x[, 2] <- -pv$x[, 2]
}
d$PC1 <- pv$x[, 1]; d$PC2 <- pv$x[, 2]

rank_tbl <- data.frame(descriptor = rownames(L), label = unname(featlab[rownames(L)]),
                       PC1_loading = round(L[, 1], 3), abs_PC1_loading = round(abs(L[, 1]), 3),
                       PC2_loading = round(L[, 2], 3), row.names = NULL) %>%
  arrange(desc(abs_PC1_loading)) %>%
  mutate(rank_on_PC1 = seq_len(dplyr::n()), plotted = rank_on_PC1 <= 2)
vx <- rank_tbl$descriptor[1]
vy <- rank_tbl$descriptor[2]

# The published variance shares (Table S1) must be the ones the axes print.
.pcaf <- "output/pca_summary.csv"
stopifnot(file.exists(.pcaf))
.s1 <- readr::read_csv(.pcaf, show_col_types = FALSE)
stopifnot(all(abs(round(ve[1:2], 1) -
                  as.numeric(.s1[.s1$feature == "Variance explained (%)", c("PC1", "PC2")])) < 0.05))

# How far each plane separates the rotation sense, as a point-biserial correlation (for a binary
# label this IS the Pearson correlation with its 0/1 coding), and how far each tracks speed.
pend <- as.numeric(d$hodoArea > 0)
off <- sum(d$PC1 < PC_XLIM[1] | d$PC1 > PC_XLIM[2] | d$PC2 < PC_YLIM[1] | d$PC2 > PC_YLIM[2])
readr::write_csv(
  rank_tbl %>% mutate(
    r_with_rotation_sense = round(abs(sapply(descriptor, function(v) cor(d[[v]], pend))), 3),
    r_with_speed = round(sapply(descriptor, function(v) cor(d[[v]], d$meanSpeed_n)), 3)) %>%
    bind_rows(data.frame(
      descriptor = c("PC1", "PC2"), label = sprintf("PC%d", 1:2),
      r_with_rotation_sense = round(c(abs(cor(d$PC1, pend)), abs(cor(d$PC2, pend))), 3),
      r_with_speed = round(c(cor(d$PC1, d$meanSpeed_n), cor(d$PC2, d$meanSpeed_n)), 3))) %>%
    mutate(steps = nrow(d), individuals = n_distinct(d$subjectID),
           var_pct_PC1 = round(ve[1], 1), var_pct_PC2 = round(ve[2], 1),
           plotted_x = vx, plotted_y = vy,
           n_pendular = sum(pend == 1), n_bouncing = sum(pend == 0),
           pc_view_x = paste(PC_XLIM, collapse = " to "),
           pc_view_y = paste(PC_YLIM, collapse = " to "), steps_off_pc_view = off),
  "output/gaitspace_figure_stats.csv")

# ---- panels -----------------------------------------------------------------------------
SPEED_RAMP <- c("#3D4A5C", "#1C9DA8", "#C0398B", "#D9A23B")
PT <- list(geom_point(alpha = 0.55, size = 0.8), theme_bio())
refline <- function(v, fn) if (v %in% names(REFLINE))
  fn(REFLINE[[v]], linetype = 3, colour = "grey55", linewidth = 0.3) else NULL

cx <- plotcol(vx); cy <- plotcol(vy)
desc_ann <- list(refline(vx, function(x, ...) geom_vline(xintercept = x, ...)),
                 refline(vy, function(y, ...) geom_hline(yintercept = y, ...)),
                 labs(x = featlab[[vx]], y = featlab[[vy]]))
pc_ann <- list(coord_cartesian(xlim = PC_XLIM, ylim = PC_YLIM),
               labs(x = sprintf("PC1 (%.0f%% of variance)", ve[1]),
                    y = sprintf("PC2 (%.0f%% of variance)", ve[2])))
by_speed <- list(scale_colour_gradientn(colours = SPEED_RAMP, name = "dimensionless\nspeed"))
by_rot <- list(scale_colour_manual(values = ROT_COLORS, name = "hodograph\nrotation sense"))

set.seed(11)
d_speed <- d[sample.int(nrow(d)), ]
d_rot <- d %>% arrange(factor(rot, levels = names(sort(table(rot), decreasing = TRUE))))

pA <- ggplot(d_speed, aes(.data[[cx]], .data[[cy]], colour = meanSpeed_n)) + PT + desc_ann + by_speed
pB <- ggplot(d_speed, aes(PC1, PC2, colour = meanSpeed_n)) + PT + pc_ann + by_speed
pC <- ggplot(d_rot, aes(.data[[cx]], .data[[cy]], colour = rot)) + PT + desc_ann + by_rot
pD <- ggplot(d_rot, aes(PC1, PC2, colour = rot)) + PT + pc_ann + by_rot

# Each top panel is drawn on the same x range as the panel beneath it, so the axis is carried
# once by the lower row and the two rows sit together.
fig <- (pA + bio_drop_x() | pB + bio_drop_x()) / (pC | pD) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right", legend.box = "vertical",
        legend.key.height = unit(22, "pt"), legend.key.width = unit(7, "pt"))
bio_save("Fig2_GaitSpace", fig, width = BIO_W2, height = 4.9)

cat(sprintf("10_fig_gaitspace_planes: %d steps, %d individuals; PC1 %.1f%%, PC2 %.1f%%.\n",
            nrow(d), n_distinct(d$subjectID), ve[1], ve[2]))
cat(sprintf("  PC1 loading rank: %s\n",
            paste(sprintf("%s %+.3f", rank_tbl$descriptor, rank_tbl$PC1_loading), collapse = ", ")))
cat(sprintf("  descriptor plane: x = %s, y = %s; PC view leaves %d of %d steps (%.2f%%) off.\n",
            vx, vy, off, nrow(d), 100 * off / nrow(d)))
cat(sprintf("  rotation-sense separation |r|: %s %.3f, %s %.3f, PC1 %.3f, PC2 %.3f\n",
            vx, abs(cor(d[[vx]], pend)), vy, abs(cor(d[[vy]], pend)),
            abs(cor(d$PC1, pend)), abs(cor(d$PC2, pend))))
