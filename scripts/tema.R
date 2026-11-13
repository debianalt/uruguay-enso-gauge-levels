# tema y función única para escribir figuras. La cargan scripts/figuras.R y cualquier figura nueva.
# Criterio: una sola tinta para los datos, un solo acento (el agua), grilla y ejes en segundo plano, rótulos
# directos en vez de leyendas, legible en escala de grises. Tamaños iguales en todas las figuras.

suppressPackageStartupMessages(library(ggplot2))

TINTA <- "#1b1b1b"; SECUNDARIA <- "#5f5f5f"; TENUE <- "#9a9a9a"; GRILLA <- "#e6e6e6"
AGUA <- "#2c6e9b"; AGUA_CLARA <- "#c9dcea"; BANDA <- "#d9d9d9"; TIERRA <- "#f4f2ee"; TIERRA_FOCO <- "#e7e2d8"
ANCHO_CM <- 17; BASE_PT <- 9
coma <- function(x, d = 1) formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = ".")
# figuras en inglés para el envío (11-nov-2026): punto decimal
punto <- function(x, d = 1) formatC(x, format = "f", digits = d, decimal.mark = ".", big.mark = ",")

tema <- function() {
  theme_minimal(base_size = BASE_PT) +
    theme(text = element_text(colour = TINTA),
          axis.text = element_text(colour = SECUNDARIA, size = BASE_PT - 1),
          axis.title = element_text(colour = TINTA, size = BASE_PT),
          axis.ticks = element_line(colour = TENUE, linewidth = 0.25), axis.ticks.length = unit(2, "pt"),
          axis.line = element_line(colour = TENUE, linewidth = 0.3),
          panel.grid.major = element_line(colour = GRILLA, linewidth = 0.25), panel.grid.minor = element_blank(),
          strip.text = element_text(colour = TINTA, size = BASE_PT, hjust = 0, face = "bold",
                                    margin = margin(b = 3)),
          panel.spacing = unit(10, "pt"), legend.position = "none",
          plot.margin = margin(4, 6, 4, 4))
}

# única función que escribe una figura: PNG a 600 ppp y PDF vectorial, con el nombre de la numeración del artículo
guardar <- function(g, n, alto_cm, ancho_cm = ANCHO_CM) {
  dir.create("figures", showWarnings = FALSE)
  nombre <- if (is.numeric(n)) sprintf("Fig%d", n) else n   # p. ej. "graphical_abstract"
  ggsave(sprintf("figures/%s.png", nombre), g, width = ancho_cm, height = alto_cm, units = "cm", dpi = 600,
         device = ragg::agg_png, bg = "white")
  ggsave(sprintf("figures/%s.pdf", nombre), g, width = ancho_cm, height = alto_cm, units = "cm", device = cairo_pdf,
         bg = "white")
  cat(sprintf("figures/%s.png y .pdf (%g × %g cm)\n", nombre, ancho_cm, alto_cm))
}
