# =============================================================================
# 12_figure1_flow.R
# Draw the two-panel MIMIC-IV / eICU-CRD participant flow diagram (3900x2400,
# 300 dpi, cairo). Counts are hard-coded from the frozen analyses; this script
# does not read patient data. Writes figures/Figure1_flow.png.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

out <- file.path(P_FIGURES, "Figure1_flow.png")
png(out, width = 3900, height = 2400, res = 300, bg = "white", type = "cairo")
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(0.25, 0.25, 1.6, 0.25), family = "sans",
    xaxs = "i", yaxs = "i", xpd = NA, lheight = 1.0)

draw_box <- function(x, y, w, h, label, fill = "#F2F4F7", border = "#344054", cex = 0.88) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2,
       col = fill, border = border, lwd = 2.2)
  text(x, y, label, cex = cex, col = "#17212B", font = 1)
}
arrow_down <- function(x, y1, y2, col = "#344054") {
  arrows(x, y1, x, y2, length = 0.12, angle = 22, lwd = 2, col = col)
}
arrow_side <- function(x1, x2, y, col = "#7B8794") {
  arrows(x1, y, x2, y, length = 0.10, angle = 22, lwd = 1.8, col = col)
}
panel_setup <- function(title) {
  plot.new(); plot.window(xlim = c(0, 11), ylim = c(0.4, 10))
  title(main = title, cex.main = 1.45, font.main = 2, col.main = "#17212B",
        line = 0.2)
}

panel_setup("MIMIC-IV")
draw_box(4.0, 9.1, 6.4, 1.0,
         "Screened major abdominal surgery stays\nwith a preoperative platelet count\nn = 3,609 | deaths = 186", fill = "#E8EEF5", cex = 0.78)
arrow_side(7.2, 7.75, 9.1)
draw_box(9.35, 9.1, 3.1, 1.05,
         "Excluded after clinical\ncriteria\nn = 1,233 | deaths = 5", fill = "#F7F7F7", border = "#98A2B3", cex = 0.70)
arrow_down(4.0, 8.58, 8.05)
draw_box(4.0, 7.55, 6.4, 1.0,
         "Eligible cohort after exclusions\nn = 2,376 | deaths = 181", fill = "#DCE9F3", cex = 0.86)
arrow_side(7.2, 7.75, 7.55)
draw_box(9.35, 7.55, 3.1, 1.05,
         "No qualifying panel\nn = 1,490 | deaths = 94", fill = "#F7F7F7", border = "#98A2B3", cex = 0.72)
arrow_down(4.0, 7.02, 6.22)
draw_box(4.0, 5.7, 6.4, 1.0,
         "Qualifying same-panel differential CBC\nn = 886 | deaths = 87", fill = "#E5F1EC", cex = 0.84)
arrow_side(7.2, 7.75, 5.70)
draw_box(9.35, 5.70, 3.1, 1.05,
         "No usable monocyte\nn = 7 | deaths = 0", fill = "#F7F7F7", border = "#98A2B3", cex = 0.72)
arrow_down(4.0, 5.17, 4.28)
draw_box(4.0, 3.75, 6.4, 1.0,
         "Complete five-index analysis sample\nn = 879 | deaths = 87", fill = "#D6E9DE", cex = 0.86)
text(4.0, 2.45, "NLR/SII comparison: n = 886; SIRI/PIV comparison: n = 879",
     cex = 0.68, col = "#475467")
box()

panel_setup("eICU-CRD")
draw_box(4.0, 9.1, 6.4, 1.0,
         "First adult surgical ICU stays\nn = 9,824", fill = "#E8EEF5", cex = 0.86)
arrow_down(4.0, 8.58, 8.05)
draw_box(4.0, 7.55, 6.4, 1.0,
         "Complete five-item, same-offset panels\nn = 3,024 stays | deaths = 287", fill = "#DCE9F3", cex = 0.80)
arrow_side(7.2, 7.75, 7.55)
draw_box(9.35, 7.55, 3.1, 0.9,
         "Invalid numeric panel\nn = 20 | deaths = 4", fill = "#F7F7F7", border = "#98A2B3", cex = 0.72)
arrow_down(4.0, 7.02, 6.55)
draw_box(4.0, 6.0, 6.4, 1.0,
         "Valid numeric panels\nn = 3,004 stays | deaths = 283", fill = "#E5F1EC", cex = 0.84)
arrow_down(4.0, 5.47, 4.98)
draw_box(4.0, 4.45, 6.4, 1.0,
         "One index stay per patient\nn = 2,866 | deaths = 277", fill = "#E5F1EC", cex = 0.84)
arrow_side(7.2, 7.75, 4.45)
draw_box(9.35, 4.45, 3.1, 0.9,
         "Zero absolute lymphocytes\nn = 11 | deaths = 1", fill = "#F7F7F7", border = "#98A2B3", cex = 0.70)
arrow_down(4.0, 3.92, 3.45)
draw_box(4.0, 2.9, 6.4, 1.0,
         "Final patient-level analysis cohort\nn = 2,855 | deaths = 276", fill = "#D6E9DE", cex = 0.82)
box()

dev.off()
cat("Figure 1 written to", out, "\n")
