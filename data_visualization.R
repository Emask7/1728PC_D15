library(forcats)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)
library(openxlsx)


# CD4 Cells ------------------------------------------------------------------
  CD4_DEdata <- read.xlsx("../Dot Plot Data.xlsx",
                          sheet = "CD4+ T Cells", cols = c(2, 4:6),
                          rowNames = FALSE, colNames = TRUE)
  colnames(CD4_DEdata) <- c("Name", "LFC", "P.val", "log10p.val")

  head(CD4_DEdata)

  CD4_plot <- ggplot(CD4_DEdata,
                     aes(x = P.val, y = fct_rev(Name), color = LFC)) +
    geom_point(size = 5) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    xlab("p-value") +
    scale_x_continuous(limits = c(0.00, 0.05)) +
    ylab("") +
    ggtitle("CD4+ T Cells\nCD8 depleted vs Control") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right")
  CD4_plot

  ggsave("CD4 T Cells DEGs.png", plot = CD4_plot,
         units="in", width=6, height=10, dpi=300)


# CD8 Cells ------------------------------------------------------------------
  CD8_DEdata <- read.xlsx("../Dot Plot Data.xlsx",
                          sheet = "CD8+ T Cells", cols = c(2, 4:6),
                          rowNames = FALSE, colNames = TRUE)
  colnames(CD8_DEdata) <- c("Name", "LFC", "P.val", "log10p.val")

  head(CD8_DEdata)

  CD8_plot <- ggplot(CD8_DEdata,
                     aes(x = P.val, y = fct_rev(Name), color = LFC)) +
    geom_point(size = 5) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    xlab("p-value") +
    scale_x_continuous(limits = c(0.00, 0.05)) +
    ylab("") +
    ggtitle("CD8+ T Cells\nCD8 depleted vs Control") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right")
  CD8_plot

  ggsave("CD8 T Cells DEGs.png", plot = CD8_plot,
         units="in", width=6, height=16, dpi=300)



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
# Combine all plots  ---------------------------------------------------------
comb_plot <- grid.arrange(grobs = list(CD4_plot, CD8_plot),
                              # widths = c(2.5, 1, 1.5),
                              layout_matrix = rbind(c(1, 2))
                          )

comb_plot
ggsave("Combined Dot Plots.pdf", plot = comb_plot,
       units="in", width=10, height=7.5, dpi=300)
