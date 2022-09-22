library(forcats)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)
library(openxlsx)

CD4_DEdata <- read.xlsx("../Dot Plot Data.xlsx",
                        sheet = "CD4+ T Cells", cols = c(2, 4:6),
                        rowNames = FALSE, colNames = TRUE)
colnames(CD4_DEdata) <- c("Name", "LFC", "P.val", "log10p.val")

head(CD4_DEdata)

# CD4 Cells ------------------------------------------------------------------
  CD4_plot <- ggplot(CD4_DEdata,
                     aes(x = P.val, y = fct_rev(Name),
                         size = 5, color = LFC)) +
    geom_point() +
    # scale_size_area(limits = c(20, 200)) +
    scale_colour_gradient2(low = "blue", mid = "white",
                           high = "red", limits = c(-3.5, 3.5)) +
    xlab("FDR") +
    scale_x_continuous(limits = c(0.00, 0.05)) +
    ylab("") +
    ggtitle("CD4+ T Cells\nCD8 depleted vs Control") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  CD4_plot


# # CD8 Cells ------------------------------------------------------------------
# p_hm_CD8 <- ggplot(hallmark_CD8,
#                    aes(x = FDR, y = fct_rev(Name),
#                        size = Size, color = NES)) +
#   geom_point() +
#   scale_size_area(max_size = 5, limits = c(20, 200)) +
#   scale_colour_gradient2(low = "blue", mid = "white",
#                          high = "red", limits = c(-3, 3)) +
#   xlab("FDR") +
#   scale_x_continuous(limits = c(0.00, 0.05)) +
#   ylab("") +
#   ggtitle("CD8+ T Cells") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")
# p_hm_CD8
#
#
#
# # NK Cells ------------------------------------------------------------------
# p_hm_NK <- ggplot(hallmark_NK,
#                   aes(x = FDR, y = fct_rev(Name),
#                       size = Size, color = NES)) +
#   geom_point() +
#   scale_size_area(max_size = 5, limits = c(20, 200)) +
#   scale_colour_gradient2(low = "blue", mid = "white",
#                          high = "red", limits = c(-3, 3)) +
#   xlab("FDR") +
#   scale_x_continuous(limits = c(0.00, 0.05)) +
#   ylab("") +
#   ggtitle("NK Cells") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "right")
# p_hm_NK
#
# # Combine all plots  ---------------------------------------------------------
# hallmark_plot <- grid.arrange(grobs = list(CD4_plot, p_hm_CD8, p_hm_NK),
#                               widths = c(2.5, 1, 1.5),
#                               layout_matrix = rbind(c(1, 2, 3)))
#
# hallmark_plot
# ggsave("Hallmark GSEA Dot Plots.pdf", plot = hallmark_plot,
#        units="in", width=10, height=7.5, dpi=300)
