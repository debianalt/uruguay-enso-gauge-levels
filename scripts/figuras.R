# las cuatro figuras del artículo, numeradas en orden de aparición.
#   Fig1  mapa del tramo y de los siete hidrómetros
#   Fig2  cómo se define una crecida: un año de El Niño y un año neutro en San Javier
#   Fig3  probabilidad anual de alcanzar la cota de evacuación en función del Niño 3.4
#   Fig4  la relación con El Niño (b) bajo cada prueba de robustez
#   graphical_abstract: el mapa con la probabilidad de evacuación en 0 y +1,5 °C (13-nov-2026); no es figura del
#   cuerpo. El fondo del mapa es mapa_base(), común a la Fig1 y al resumen gráfico.
#
# Usa el modelo y el remuestreo de scripts/modelo.R, los mismos de scripts/frecuencia.R. Lo que se dibuja sale de
# tables/ o se escribe ahí (fig2_*.csv, fig3_curva.csv, fig4_b.csv). Controles que cortan la corrida:
#   - la curva de la Fig3 en los valores de la tabla coincide con tables/prob_frecuencia_sola.csv (la tabla
#     principal desde el 11-nov-2026: frecuencia sola, d = 0);
#   - las crecidas marcadas en la Fig2 son las de data/crecidas_evento.csv.
# El bootstrap de la banda de la Fig3 es un sorteo aparte del de la tabla: sus bordes pueden diferir en centésimas.
#
# Rótulos en inglés y punto decimal desde el 11-nov-2026 (envío a Atmósfera).
#
# Uso (desde la raíz del repositorio):  Rscript scripts/figuras.R

source("scripts/modelo.R", encoding = "UTF-8")
source("scripts/tema.R", encoding = "UTF-8")
suppressPackageStartupMessages({library(sf); library(ggrepel); library(patchwork)})
sf_use_s2(FALSE)
M0 <- 3; N_BOOT <- 500

# ------------------------------------------------------------------------------------------------- Fig1: mapa
ne <- function(capa) suppressWarnings(st_read(file.path("data/geo", capa, paste0(capa, ".shp")), quiet = TRUE))
paises <- ne("ne_10m_admin_0_countries"); prov <- ne("ne_10m_admin_1_states_provinces")
rios <- ne("ne_10m_rivers_lake_centerlines")
rio_ign <- st_read("data/geo/rio_uruguay_ign.geojson", quiet = TRUE)
hid <- st_read("data/geo/hidrometros.geojson", quiet = TRUE)
caja <- c(xmin = -57.3, xmax = -52.4, ymin = -29.5, ymax = -26.2)
recorte <- st_as_sfc(st_bbox(caja, crs = 4326))
sa <- paises[paises$CONTINENT == "South America", ]
foco <- prov[prov$iso_a2 == "AR" & prov$name %in% c("Misiones", "Corrientes"), ]
prov_reg <- suppressWarnings(st_intersection(prov[prov$iso_a2 %in% c("AR", "BR", "PY", "UY"), "name"], recorte))
rio_ne <- suppressWarnings(st_intersection(rios[!is.na(rios$name) & rios$name %in% c("Uruguay", "Paraná", "Iguaçu"),
                                                c("name")], recorte))
rio_uy <- rio_ign[rio_ign$nombre == "Río Uruguay", ]
xy <- st_coordinates(hid); hid$lon <- xy[, 1]; hid$lat <- xy[, 2]
hid$rotulo <- ifelse(is.na(hid$evacuacion_m), hid$nombre, sprintf("%s · %s m", hid$nombre, punto(hid$evacuacion_m)))
km50 <- 50 / (111.32 * cos(28.9 * pi / 180))
etiquetas <- data.frame(
  x = c(-53.3, -56.55, -56.75, -54.75, -56.75), y = c(-27.85, -26.6, -28.85, -26.75, -28.35),
  t = c("BRAZIL", "PARAGUAY", "ARGENTINA", "Misiones", "Corrientes"), tipo = c("pais", "pais", "pais", "prov", "prov"))

# fondo del mapa (tierra, provincias foco, ríos, rótulos de país y provincia), común a la Fig1 y al resumen gráfico
mapa_base <- function(caja, etiquetas, rio_uy_xy = c(-53.05, -27.0), parana_xy = c(-55.62, -27.2)) {
  ggplot() +
  geom_sf(data = paises, fill = TIERRA, colour = NA) +
  geom_sf(data = foco, fill = TIERRA_FOCO, colour = NA) +
  geom_sf(data = foco, fill = NA, colour = "white", linewidth = 0.4) +
  geom_sf(data = paises, fill = NA, colour = TENUE, linewidth = 0.3) +
  geom_sf(data = rio_ne[rio_ne$name != "Uruguay", ], colour = AGUA_CLARA, linewidth = 0.5) +
  geom_sf(data = rio_ne[rio_ne$name == "Uruguay", ], colour = AGUA, linewidth = 0.9, lineend = "round") +
  geom_sf(data = rio_uy, colour = AGUA, linewidth = 0.9, lineend = "round") +
  geom_text(data = etiquetas[etiquetas$tipo == "pais", ], aes(x, y, label = t), colour = TENUE, size = 2.9,
            family = "sans", fontface = "plain") +
  geom_text(data = etiquetas[etiquetas$tipo == "prov", ], aes(x, y, label = t), colour = SECUNDARIA, size = 2.7,
            fontface = "italic") +
  annotate("text", x = rio_uy_xy[1], y = rio_uy_xy[2], label = "Uruguay River", colour = AGUA, size = 2.7,
           fontface = "italic") +
  annotate("text", x = parana_xy[1], y = parana_xy[2], label = "Paraná River", colour = "#7fa8c6", size = 2.5,
           fontface = "italic", angle = 33) +
  coord_sf(xlim = caja[c("xmin", "xmax")], ylim = caja[c("ymin", "ymax")], expand = FALSE) +
  labs(x = NULL, y = NULL) + tema() +
  theme(panel.background = element_rect(fill = "white", colour = NA), panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA, colour = TENUE, linewidth = 0.3), axis.line = element_blank())
}

mapa <- mapa_base(caja, etiquetas) +
  geom_point(data = hid, aes(lon, lat), shape = 21, fill = "white", colour = TINTA, size = 2.1, stroke = 0.7) +
  geom_text(data = hid, aes(lon + 0.11, lat - 0.085, label = rotulo), size = 2.75, colour = TINTA, hjust = 0,
            vjust = 0.5) +
  annotate("segment", x = -55.2, xend = -55.2 + km50, y = caja["ymin"] + 0.22, yend = caja["ymin"] + 0.22,
           colour = TINTA, linewidth = 0.5) +
  annotate("text", x = -55.2 + km50 / 2, y = caja["ymin"] + 0.33, label = "50 km", size = 2.5,
           colour = TINTA) +
  scale_x_continuous(breaks = seq(-57, -53, 1)) + scale_y_continuous(breaks = seq(-29, -27, 1))

situacion <- ggplot() +
  geom_sf(data = sa, fill = "#e9e6e0", colour = "white", linewidth = 0.15) +
  geom_sf(data = recorte, fill = NA, colour = TINTA, linewidth = 0.45) +
  coord_sf(xlim = c(-82, -34), ylim = c(-56, 13), expand = FALSE) + theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = TENUE, linewidth = 0.3))
fig1 <- mapa + inset_element(situacion, left = 0.80, bottom = 0.03, right = 0.985, top = 0.42, align_to = "panel")
guardar(fig1, 1, alto_cm = 12.2)

# --------------------------------------------------------------------- Fig2: qué es una crecida (San Javier)
SID <- 65
dia <- leer("diario"); dia$fecha <- as.Date(dia$fecha); dia <- dia[dia$serie == SID, ]
dia$ah <- anio_h(dia$fecha, M0); enso$ah <- anio_h(enso$mes, M0)
xa <- tapply(enso$x6, enso$ah, mean); na <- table(dia$ah)
cand <- data.frame(ah = as.integer(names(na)), n = as.integer(na))
cand$x <- xa[as.character(cand$ah)]; cand <- cand[cand$n >= 347 & cand$ah >= 1965 & !is.na(cand$x), ]
# regla, fijada antes de mirar las series: desde 1965, años hidrológicos con el 95 % de los días; el de El Niño es
# el de mayor media del índice y el neutro el de media más cercana a cero.
a_nino <- cand$ah[which.max(cand$x)]; a_neutro <- cand$ah[which.min(abs(cand$x))]
u <- serie[as.character(SID), "u_m"]; z_al <- serie[as.character(SID), "alerta_m"]
z_ev <- serie[as.character(SID), "evacuacion_m"]
evs <- ev[ev$serie == SID & ev$mes_valido, ]; evs$ah <- anio_h(evs$fecha_pico, M0)
panel <- function(a, clase) {
  d <- dia[dia$ah == a, ]; e <- evs[evs$ah == a, ]
  tit <- sprintf("%d–%02d · %s (mean index %s °C) · %d flood%s", a, (a + 1) %% 100, clase,
                 sub("^([1-9])", "+\\1", punto(cand$x[cand$ah == a], 1)), nrow(e), ifelse(nrow(e) == 1, "", "s"))
  list(d = transform(d, panel = tit, dia_h = as.integer(fecha - as.Date(sprintf("%d-03-01", a)))),
       e = if (nrow(e)) transform(e, panel = tit, dia_h = as.integer(fecha_pico - as.Date(sprintf("%d-03-01", a)))) else NULL)
}
p1 <- panel(a_nino, "El Niño"); p2 <- panel(a_neutro, "neutral")
dd <- rbind(p1$d, p2$d); ee <- rbind(p1$e, p2$e)
dd$panel <- factor(dd$panel, levels = unique(dd$panel)); if (!is.null(ee)) ee$panel <- factor(ee$panel, levels = levels(dd$panel))
stopifnot(nrow(p1$e) == sum(evs$ah == a_nino), all(ee$pico_m >= u))
write.csv(dd[, c("serie", "fecha", "pref", "ah")], "tables/fig2_serie.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ee[, c("serie", "fecha_pico", "pico_m", "ah")], "tables/fig2_crecidas.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
meses <- data.frame(dia_h = cumsum(c(0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31)),
                    m = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb"))
niveles <- data.frame(y = c(u, z_al, z_ev), t = c(sprintf("threshold, %s m", punto(u, 2)),
  sprintf("alert, %s m", punto(z_al)), sprintf("evacuation, %s m", punto(z_ev))), lt = c("dotted", "dashed", "solid"),
  panel = factor(levels(dd$panel)[1], levels = levels(dd$panel)))
ymax <- max(c(dd$pref, z_ev), na.rm = TRUE) + 1.1
fig2 <- ggplot(dd, aes(dia_h, pref)) +
  geom_ribbon(aes(ymin = pmin(pref, u), ymax = pref), fill = AGUA_CLARA) +
  geom_hline(data = transform(niveles, panel = NULL), aes(yintercept = y, linetype = lt), colour = SECUNDARIA,
             linewidth = 0.3) +
  geom_line(colour = AGUA, linewidth = 0.45) +
  geom_point(data = ee, aes(dia_h, pico_m), shape = 21, fill = "white", colour = TINTA, size = 1.9, stroke = 0.7) +
  geom_text(data = niveles, aes(x = 2, y = y, label = t), hjust = 0, vjust = -0.45, size = 2.6, colour = SECUNDARIA) +
  facet_wrap(~panel, ncol = 1) + scale_linetype_identity() +
  scale_x_continuous(breaks = meses$dia_h + 15, labels = meses$m, expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 14, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Height at the San Javier gauge (m)") + tema() +
  theme(panel.grid.major.x = element_blank(), axis.ticks.x = element_blank())
guardar(fig2, 2, alto_cm = 11.5)
cat(sprintf("Fig2: año de El Niño %d (x = %.2f, %d crecidas), neutro %d (x = %.2f, %d crecidas)\n", a_nino,
            cand$x[cand$ah == a_nino], nrow(p1$e), a_neutro, cand$x[cand$ah == a_neutro], NROW(p2$e)))

# ------------------------------------------------- Fig3: probabilidad de evacuación en función del Niño 3.4
set.seed(2026)
SIDS <- c(61, 65, 67, 68)                      # los cuatro con cota de evacuación que sostienen la relación solos
grilla <- seq(-1.5, 2.5, by = 0.1)
tabla <- read.csv("tables/prob_frecuencia_sola.csv", fileEncoding = "UTF-8")
x_ult <- max(tabla$x)
curvas <- list(); control <- list()
for (sid in SIDS) {
  D <- datos(sid, "x6"); f <- ajustar(D)
  z <- serie[as.character(sid), "evacuacion_m"]; nombre <- serie[as.character(sid), "nombre"]
  # banda con la frecuencia sola: mismo remuestreo por años hidrológicos que scripts/frecuencia_sola.R
  ye <- anio_h(D$ye, M0); ym <- anio_h(D$mm, M0); anios <- sort(unique(ym))
  idx_e <- split(seq_along(ye), ye); idx_m <- split(seq_along(ym), ym)
  pm <- matrix(NA_real_, length(grilla), N_BOOT)
  for (r in seq_len(N_BOOT)) {
    s <- remuestrear(idx_e, idx_m, sample(anios, length(anios), replace = TRUE))
    if (length(s$ie) < 10) next
    fr <- ajustar(list(xe = D$xe[s$ie], y = D$y[s$ie], xm = D$xm[s$im]))
    pm[, r] <- prob(z, D$u, fr$a, fr$b, fr$c0, 0, fr$xi0, grilla)
  }
  curvas[[length(curvas) + 1]] <- data.frame(serie = sid, nombre = nombre, z_m = z, x = grilla,
    p = prob(z, D$u, f$a, f$b, f$c0, 0, f$xi0, grilla),
    p_lo = apply(pm, 1, quantile, 0.05, na.rm = TRUE), p_hi = apply(pm, 1, quantile, 0.95, na.rm = TRUE), en_rango = grilla <= max(D$xm))
  t <- tabla[tabla$serie == sid & tabla$nivel == "evacuación", ]
  control[[length(control) + 1]] <- data.frame(serie = sid, nombre = nombre, x = t$x, tabla = t$p_frec,
    figura = prob(z, D$u, f$a, f$b, f$c0, 0, f$xi0, t$x))
}
cv <- do.call(rbind, curvas); ct <- do.call(rbind, control)
dif <- max(abs(ct$tabla - ct$figura)); stopifnot(dif < 1e-8)
cat(sprintf("Fig3, control: diferencia máxima con tables/prob_frecuencia_sola.csv = %.1e\n", dif))
write.csv(cv, "tables/fig3_curva.csv", row.names = FALSE, fileEncoding = "UTF-8")
# Garruchos lleva la marca de su reserva (serie no homogénea en aguas altas; nota de la Tabla II)
rot <- function(n, z) sprintf("%s · %s m%s", n, punto(z), ifelse(n == "Garruchos", " †", ""))
# la media más alta del índice en un año hidrológico completo desde 1907: más allá, la curva extrapola
xa_h <- tapply(enso$x6, enso$ah, mean); n_h <- tapply(!is.na(enso$x6), enso$ah, sum)
x_tope <- max(xa_h[n_h == 12 & as.integer(names(xa_h)) >= 1907])
cat(sprintf("Fig3: media anual más alta del índice desde 1907 = %.2f\n", x_tope)); stopifnot(abs(x_tope - 1.4) < 0.05)
cv$fuera <- cv$x >= x_tope - 1e-9
cv$panel <- factor(rot(cv$nombre, cv$z_m), levels = unique(rot(cv$nombre, cv$z_m)))
ct$z_m <- serie[as.character(ct$serie), "evacuacion_m"]; ct$panel <- factor(rot(ct$nombre, ct$z_m), levels = levels(cv$panel))
ct$lab <- sprintf("%d", round(ct$tabla * 100))
marca <- data.frame(panel = factor(levels(cv$panel)[1], levels = levels(cv$panel)), x = x_ult, y = 6,
                    t = "October\n2026")
marca_tope <- data.frame(panel = factor(levels(cv$panel)[2], levels = levels(cv$panel)), x = x_tope, y = 7,
                         t = "beyond any\nobserved year")
# la curva pasa por x_tope en los dos tramos: se interpola el punto de unión
union <- do.call(rbind, lapply(split(cv, cv$panel), function(d) {
  data.frame(panel = d$panel[1], x = x_tope, p = approx(d$x, d$p, x_tope)$y)
}))
dentro <- rbind(cv[!cv$fuera, c("panel", "x", "p")], union); dentro <- dentro[order(dentro$panel, dentro$x), ]
afuera <- rbind(union, cv[cv$fuera, c("panel", "x", "p")]); afuera <- afuera[order(afuera$panel, afuera$x), ]
fig3 <- ggplot(cv, aes(x, p * 100)) +
  annotate("rect", xmin = x_tope, xmax = Inf, ymin = 0, ymax = 100, fill = TINTA, alpha = 0.045) +
  geom_vline(xintercept = x_ult, colour = TENUE, linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = p_lo * 100, ymax = p_hi * 100), fill = BANDA) +
  geom_line(data = dentro, colour = TINTA, linewidth = 0.6) +
  geom_line(data = afuera, colour = TINTA, linewidth = 0.6, linetype = "22") +
  geom_text(data = marca_tope, aes(x, y, label = t), size = 2.5, colour = SECUNDARIA, hjust = -0.06, lineheight = 0.9) +
  geom_point(data = ct, aes(x, tabla * 100), shape = 21, fill = "white", colour = TINTA, size = 1.7, stroke = 0.6) +
  geom_text(data = ct, aes(x, tabla * 100, label = lab), size = 2.6, colour = TINTA, hjust = 1.35, vjust = -0.35) +
  geom_text(data = marca, aes(x, y, label = t), size = 2.5, colour = SECUNDARIA, hjust = 1.08, lineheight = 0.9) +
  facet_wrap(~panel, nrow = 2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-1, 2, 1), labels = function(v) sub("^([1-9])", "+\\1", sub("-", "–", format(v))),
                     expand = c(0.01, 0)) +
  labs(x = "Niño 3.4 index, six-month mean (°C)",
       y = "Annual probability of reaching the evacuation level (%)") + tema()
guardar(fig3, 3, alto_cm = 12.5)

# ------------------------------------------------------------------- graphical abstract: el mapa
# (13-nov-2026; elegido por el autor esa noche en lugar de la Fig3 reducida) probabilidad de alcanzar la cota de evacuación con el índice en 0 y en +1,5 °C, de la Tabla II.
# Llevan número sólo los cuatro hidrómetros de la Fig3 (SIDS), los que sostienen la relación solos; Alba Posse,
# Panambí y Pepirí Miní van como punto con nombre. Puntos sobre el río, nunca áreas: la unidad es el hidrómetro.
ga_t <- tabla[tabla$nivel == "evacuación" & tabla$serie %in% SIDS & tabla$x %in% c(0, 1.5), ]
stopifnot(nrow(ga_t) == 2 * length(SIDS))
pct <- function(sid, x) sprintf("%.0f%%", 100 * ga_t$p_frec[ga_t$serie == sid & ga_t$x == x])
ga <- hid
ga$con_num <- ga$serie %in% SIDS
ga$nom <- ifelse(ga$serie == 67, paste0(ga$nombre, "†"), ga$nombre)
ga$l1 <- ifelse(is.na(ga$evacuacion_m), ga$nom, sprintf("%s · %s m", ga$nom, punto(ga$evacuacion_m)))
ga$l2 <- ifelse(ga$con_num, vapply(ga$serie, function(s) paste(pct(s, 0), "→", pct(s, 1.5)), ""), NA)
caja_ga <- c(xmin = -56.75, xmax = -52.85, ymin = -29.05, ymax = -26.75)
etq_ga <- data.frame(x = c(-53.35, -56.35, -55.0, -55.75, -56.35), y = c(-28.75, -26.95, -27.2, -27.55, -28.25),
                     t = c("BRAZIL", "PARAGUAY", "ARGENTINA", "Misiones", "Corrientes"),
                     tipo = c("pais", "pais", "pais", "prov", "prov"))
# rótulos abajo a la derecha del punto (lado brasileño, perpendicular al cauce, que corre de SO a NE)
dx <- 0.07; dy <- 0.05; dl <- 0.10
gab <- mapa_base(caja_ga, etq_ga, rio_uy_xy = c(-53.2, -26.97), parana_xy = c(-56.0, -27.35)) +
  geom_point(data = ga, aes(lon, lat, fill = con_num), shape = 21, colour = TINTA, size = 2.3, stroke = 0.7) +
  scale_fill_manual(values = c(`TRUE` = AGUA, `FALSE` = "white")) +
  geom_text(data = ga, aes(lon + dx, lat - dy, label = l1), size = 2.75, colour = TINTA, hjust = 0, vjust = 1) +
  geom_text(data = ga[ga$con_num, ], aes(lon + dx, lat - dy - dl, label = l2), size = 3.1, colour = AGUA,
            fontface = "bold", hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  labs(title = "Annual probability of reaching the evacuation level",
       subtitle = "Niño 3.4 index (six-month mean) at 0 °C  →  at +1.5 °C",
       caption = paste("Argentine gauges on the Uruguay River, the Argentina–Brazil border. Conditional on the index over a",
                       "whole hydrological year; not a forecast.",
                       "Open circles: no official level, or the relationship does not hold at the gauge alone.",
                       "† Series not homogeneous at high water.", sep = "\n")) +
  theme(plot.title = element_text(size = BASE_PT + 1, face = "bold", colour = TINTA, margin = margin(b = 2)),
        plot.subtitle = element_text(size = BASE_PT - 0.5, colour = AGUA, margin = margin(b = 4)),
        plot.caption = element_text(size = BASE_PT - 1.5, colour = SECUNDARIA, hjust = 0, lineheight = 1.1),
        plot.title.position = "plot", plot.caption.position = "plot")
guardar(gab, "graphical_abstract", alto_cm = 10, ancho_cm = 13)
cat("Resumen gráfico (mapa):", paste(ga$nombre[ga$con_num], ga$l2[ga$con_num], collapse = "; "), "\n")

# ------------------------------------------------------------ Fig4: b bajo cada prueba de robustez (conjunto)
lee <- function(f) read.csv(file.path("tables", f), fileEncoding = "UTF-8")
aj <- lee("ajuste.csv"); pr <- aj[aj$serie == "uruguay_conjunto" & aj$variante == "principal", ]
um <- lee("robustez_umbral.csv"); um <- um[um$serie == "uruguay_conjunto", ]
re <- lee("robustez_represas.csv"); re <- re[re$serie == "uruguay_conjunto" & re$corte == "2000-03-01", ]
te <- lee("robustez_tendencia.csv"); t65 <- lee("robustez_tendencia_desde1965.csv")
k <- lee("robustez_1930_1966.csv"); av <- lee("robustez_1941_1982.csv"); es <- lee("robustez_estacion.csv")
ix <- lee("robustez_indice.csv"); ixf <- function(col) ix[ix$indice == col & ix$desde == "", ]
ix48 <- function(col) ix[ix$indice == col & ix$desde != "", ]
um <- um[!(um$percentil == 98 & um$dias == 10), ]     # es el ajuste principal: no se repite (12-nov-2026)
fila <- function(grupo, caso, b, lo, hi, n = NA) data.frame(grupo, caso, b, lo, hi, n)
fi <- function(caso, r) fila("Index", caso, r$b, r$b_boot_lo, r$b_boot_hi, r$n_eventos)
fb <- rbind(
  fila("Main", "98th percentile, 10 days, full series", pr$b, pr$b_boot_lo, pr$b_boot_hi, pr$n_eventos),
  do.call(rbind, lapply(seq_len(nrow(um)), function(i) fila("Flood\ndefinition",
    sprintf("%dth percentile, %d days", um$percentil[i], um$dias[i]), um$b[i], um$b_boot_lo[i], um$b_boot_hi[i],
    um$n_eventos[i]))),
  fi("lagged three months", ixf("x6_lag3")), fi("three-month mean", ixf("x3")), fi("twelve-month mean", ixf("x12")),
  fi("HadISST, from July 1948", ix48("x6")), fi("ERSST v6, from July 1948", ix48("x6_ersst")),
  fila("River\nregulation", "before March 2000", re$b_antes, re$b_antes_lo, re$b_antes_hi, re$n_antes),
  fila("River\nregulation", "from March 2000", re$b_despues, re$b_despues_lo, re$b_despues_hi, re$n_despues),
  fila("Trend", "with time term", te$b, te$b_boot_lo, te$b_boot_hi, pr$n_eventos),
  fila("Trend", "with time term, from 1965", t65$b, t65$b_boot_lo, t65$b_boot_hi, t65$n_eventos),
  fila("Season", "annual cycle, two harmonics", es$b, es$b_boot_lo, es$b_boot_hi, es$n_eventos),
  fila("Constancy of\nthe relationship", "1930–1966", k$b_1, k$b_1_lo, k$b_1_hi, k$n_1),
  fila("Constancy of\nthe relationship", "from 1967", k$b_2, k$b_2_lo, k$b_2_hi, k$n_2),
  fila("Constancy of\nthe relationship", "1941–1982", av$b_1, av$b_1_lo, av$b_1_hi, av$n_1),
  fila("Constancy of\nthe relationship", "from 1983", av$b_2, av$b_2_lo, av$b_2_hi, av$n_2))
write.csv(fb, "tables/fig4_b.csv", row.names = FALSE, fileEncoding = "UTF-8")
fb$grupo <- factor(fb$grupo, levels = unique(fb$grupo)); fb$orden <- rev(seq_len(nrow(fb)))
fb$rot <- sprintf("%s  (%s)", fb$caso, format(fb$n, big.mark = ",", decimal.mark = ".", trim = TRUE))
fb$es_pr <- fb$grupo == "Main"
fig4 <- ggplot(fb, aes(b, reorder(rot, orden))) +
  annotate("rect", xmin = pr$b_boot_lo, xmax = pr$b_boot_hi, ymin = -Inf, ymax = Inf, fill = "#efefef") +
  geom_vline(xintercept = 0, colour = TINTA, linewidth = 0.35) +
  geom_vline(xintercept = pr$b, colour = TENUE, linewidth = 0.3, linetype = "dashed") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = TINTA, linewidth = 0.45) +
  geom_point(aes(fill = es_pr), shape = 21, colour = TINTA, size = 2, stroke = 0.6) +
  scale_fill_manual(values = c(`TRUE` = TINTA, `FALSE` = "white")) +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = function(v) punto(v, 2), limits = c(-0.08, 1.0),
                     sec.axis = sec_axis(~exp(.), name = "Factor on flood frequency, per °C",
                                         breaks = c(1, 1.5, 2, 2.5), labels = function(v) punto(v, 1))) +
  labs(x = "b, per °C of the Niño 3.4 index", y = NULL) + tema() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1, face = "bold",
                                                                      margin = margin(r = 6)),
        panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = TINTA, size = BASE_PT - 1),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing = unit(7, "pt"))
guardar(fig4, 4, alto_cm = 12.2)
