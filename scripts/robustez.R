# robustez de la b del conjunto: ¿aguanta otro umbral, las represas y una tendencia?
#
# Las reglas (qué se prueba y cuándo pasa) están en el registro de decisiones del autor, entrada del 10-nov-2026,
# escritas antes de correr. Acá no se decide nada: se calcula y se compara contra esas reglas.
#
# Lo que se afirma en la comunicación es la frecuencia: con El Niño hay más crecidas sobre el umbral. Ese número
# es b, y b sale sólo de la parte de Poisson del modelo (cuántos eventos por año). Por eso el bootstrap de estas
# pruebas reajusta la tasa y no la GPD.
#
#   paso previo  rearmar los eventos desde la serie diaria y comprobar que salen los del principal
#   prueba 1     umbral P97, P98, P99 y separación de 5, 10 y 20 días
#   prueba 2     antes y después de las represas brasileñas (corte al empezar el año hidrológico 2000)
#   prueba 3     una tendencia en el tiempo al lado de El Niño
#
# Uso (desde la raíz del repositorio):  Rscript scripts/robustez.R [--boot 500]

args <- commandArgs(trailingOnly = TRUE)
N_BOOT <- if ("--boot" %in% args) as.integer(args[which(args == "--boot") + 1]) else 500L
source("scripts/modelo.R", encoding = "UTF-8")
set.seed(2026)

B_LO <- 0.30; B_HI <- 0.70          # IC90 de la b del principal (tables/ajuste.csv, uruguay_conjunto)
M0 <- 3                             # el año hidrológico del Uruguay empieza en marzo
q5 <- function(v) unname(quantile(v, 0.05)); q95 <- function(v) unname(quantile(v, 0.95))

diario <- read.csv("data/crecidas_diario.csv", stringsAsFactors = FALSE)
diario$fecha <- as.Date(diario$fecha)

# ---------- eventos desde la serie diaria (la misma regla que armó data/crecidas_evento.csv) -----------
# u = percentil pct de todas las lecturas; un evento nuevo empieza cuando pasan más de `gap` días sin una
# lectura sobre u; un mes es válido si tiene lectura en al menos el 80 % de sus días.
armar <- function(sid, pct, gap) {
  s <- diario[diario$serie == sid & !is.na(diario$pref), ]; s <- s[order(s$fecha), ]
  u <- unname(quantile(s$pref, pct / 100, type = 7))            # type 7 = np.percentile
  cal <- seq(min(s$fecha), max(s$fecha), by = "day")
  tot <- table(format(cal, "%Y-%m-01")); con <- table(format(s$fecha, "%Y-%m-01"))
  m <- data.frame(mes = as.Date(names(tot)), total = as.integer(tot))
  m$dias <- as.integer(con[names(tot)]); m$dias[is.na(m$dias)] <- 0L
  m$valido <- m$dias >= 0.80 * m$total
  sobre <- s[s$pref > u, ]
  id <- cumsum(c(FALSE, diff(as.integer(sobre$fecha)) > gap + 1))
  e <- do.call(rbind, lapply(split(sobre, id), function(x) {
    i <- which.max(x$pref)
    data.frame(inicio = min(x$fecha), fin = max(x$fecha), fecha_pico = x$fecha[i], pico_m = x$pref[i])
  }))
  e$mes_valido <- primer_dia(e$fecha_pico) %in% m$mes[m$valido]
  list(e = e, m = m, u = u)
}

# lo mismo que datos() de modelo.R, pero sobre eventos rearmados; desde/hasta recortan el período
datos_de <- function(A, desde = NULL, hasta = NULL) {
  e <- A$e[A$e$mes_valido, ]; m <- A$m[A$m$valido, ]
  if (!is.null(desde)) { e <- e[e$fecha_pico >= desde, ]; m <- m[m$mes >= desde, ] }
  if (!is.null(hasta)) { e <- e[e$fecha_pico < hasta, ]; m <- m[m$mes < hasta, ] }
  e$mes <- primer_dia(e$fecha_pico)
  e <- merge(e, enso[, c("mes", "x6")], by = "mes"); e <- e[!is.na(e$x6), ]
  m <- merge(m, enso[, c("mes", "x6")], by = "mes"); m <- m[!is.na(m$x6), ]
  list(xe = e$x6, y = e$pico_m - A$u, xm = m$x6, ye = e$fecha_pico, mm = m$mes, u = A$u)
}
recortar <- function(D, desde = NULL, hasta = NULL) {
  ie <- rep(TRUE, length(D$ye)); im <- rep(TRUE, length(D$mm))
  if (!is.null(desde)) { ie <- ie & D$ye >= desde; im <- im & D$mm >= desde }
  if (!is.null(hasta)) { ie <- ie & D$ye < hasta; im <- im & D$mm < hasta }
  list(xe = D$xe[ie], y = D$y[ie], xm = D$xm[im], ye = D$ye[ie], mm = D$mm[im], u = D$u)
}

# b del conjunto: una sola b para todos los hidrómetros, la a propia de cada uno (perfilada, fórmula cerrada)
b_conj <- function(Ds) {
  if (!length(Ds)) return(NA_real_)
  optimize(function(b) sum(vapply(Ds, function(D) pnll_pois(b, D$xe, D$xm), 0)), c(-6, 6), tol = 1e-10)$minimum
}
con_eventos <- function(Ds, k) Ds[vapply(Ds, function(D) length(D$xe) >= k, TRUE)]

# bootstrap por años hidrológicos enteros, todas las estaciones a la vez (la misma crecida pasa por las siete)
boot_anios <- function(Ds, estad, n = N_BOOT) {
  ye <- lapply(Ds, function(D) anio_h(D$ye, M0)); ym <- lapply(Ds, function(D) anio_h(D$mm, M0))
  anios <- sort(unique(unlist(ym)))
  idx <- lapply(seq_along(Ds), function(i) list(e = split(seq_along(ye[[i]]), ye[[i]]),
                                                 m = split(seq_along(ym[[i]]), ym[[i]])))
  out <- vector("list", n)
  for (r in seq_len(n)) {
    smp <- sample(anios, length(anios), replace = TRUE)
    Db <- lapply(seq_along(Ds), function(i) {
      s <- remuestrear(idx[[i]]$e, idx[[i]]$m, smp)
      list(xe = Ds[[i]]$xe[s$ie], xm = Ds[[i]]$xm[s$im], ye = Ds[[i]]$ye[s$ie], mm = Ds[[i]]$mm[s$im])
    })
    out[[r]] <- estad(Db)
  }
  do.call(rbind, out)
}

nombre_de <- function(sid) serie[as.character(sid), "nombre"]
niveles_de <- function(sid) c(alerta = serie[as.character(sid), "alerta_m"],
                              "evacuación" = serie[as.character(sid), "evacuacion_m"])
veredicto <- list()
anotar <- function(prueba, regla, valor, pasa) {
  veredicto[[length(veredicto) + 1]] <<- data.frame(prueba = prueba, regla = regla, valor = valor, pasa = pasa)
  cat(sprintf("  %s | %s | %s | %s\n", prueba, regla, valor, if (pasa) "PASA" else "NO PASA"))
}

# ---------- paso previo: con P98 y 10 días tienen que salir los eventos del principal ----------------------------
cat("== paso previo: eventos rearmados contra data/crecidas_evento.csv y crecidas_mes.csv\n")
mal <- 0
for (sid in URUGUAY) {
  A <- armar(sid, 98, 10)
  e0 <- ev[ev$serie == sid, ]; e0 <- e0[order(e0$fecha_pico), ]
  m0 <- mes[mes$serie == sid, ]; m0 <- m0[order(m0$mes), ]
  ok_e <- nrow(A$e) == nrow(e0) && all(A$e$fecha_pico == e0$fecha_pico) && all(abs(A$e$pico_m - e0$pico_m) < 1e-9) &&
    all(abs(A$u - e0$u_m) < 1e-9) && all(A$e$mes_valido == e0$mes_valido)
  ok_m <- nrow(A$m) == nrow(m0) && all(A$m$mes == m0$mes) && all(A$m$dias == m0$dias) && all(A$m$valido == m0$valido)
  cat(sprintf("  %d %-12s eventos %d / %d %s; meses %d / %d %s\n", sid, nombre_de(sid), nrow(A$e), nrow(e0),
              if (ok_e) "iguales" else "DISTINTOS", nrow(A$m), nrow(m0), if (ok_m) "iguales" else "DISTINTOS"))
  mal <- mal + !ok_e + !ok_m
}
if (mal > 0) { cat("el paso previo no pasa: se para acá\n"); quit(status = 1) }

# ---------- prueba 1: umbral y separación de eventos -----------------------------------------------------------
cat("== prueba 1: umbral y separación\n")
ref <- read.csv("tables/prob_condicional.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
f1 <- list(); p1 <- list()
for (pct in c(97, 98, 99)) for (gap in c(5, 10, 20)) {
  As <- lapply(URUGUAY, armar, pct = pct, gap = gap)
  Ds <- lapply(As, datos_de); names(Ds) <- URUGUAY
  usados <- con_eventos(Ds, 15)
  b <- b_conj(usados)
  bt <- boot_anios(usados, function(Db) b_conj(con_eventos(Db, 5)))
  pasa <- q5(bt) > 0 && b >= B_LO && b <= B_HI
  f1[[length(f1) + 1]] <- data.frame(percentil = pct, dias = gap, serie = "uruguay_conjunto",
    nombre = sprintf("Uruguay, %d hidrómetros", length(usados)), n_eventos = sum(vapply(usados, function(D) length(D$xe), 0)),
    u_m = NA, b = b, b_boot_lo = q5(bt), b_boot_hi = q95(bt), pasa = pasa)
  cat(sprintf("  P%d, %2d días: %4d eventos, b = %+.3f [%+.3f, %+.3f] %s\n", pct, gap,
              sum(vapply(usados, function(D) length(D$xe), 0)), b, q5(bt), q95(bt), if (pasa) "pasa" else "NO PASA"))
  for (sid in URUGUAY) {
    D <- Ds[[as.character(sid)]]; n <- length(D$xe)
    f1[[length(f1) + 1]] <- data.frame(percentil = pct, dias = gap, serie = as.character(sid), nombre = nombre_de(sid),
      n_eventos = n, u_m = D$u, b = if (n >= 15) b_conj(list(D)) else NA, b_boot_lo = NA, b_boot_hi = NA, pasa = NA)
    niv <- niveles_de(sid); niv <- niv[!is.na(niv)]
    if (!length(niv)) next
    f <- if (n >= 15) ajustar(D) else NULL
    for (i in seq_along(niv)) for (x in c(0, 1.5)) {
      r <- ref[ref$serie == sid & ref$nivel == names(niv)[i] & abs(ref$x - x) < 1e-9, ]
      calculable <- !is.null(f) && niv[[i]] > D$u
      p <- if (calculable) prob(niv[[i]], D$u, f$a, f$b, f$c, f$d, f$xi, x) else NA
      p_ref <- if (nrow(r) == 1) r$p else NA
      p1[[length(p1) + 1]] <- data.frame(percentil = pct, dias = gap, serie = sid, nombre = nombre_de(sid),
        nivel = names(niv)[i], z_m = niv[[i]], u_m = D$u, x = x, p = p, p_principal = p_ref,
        dif_puntos = 100 * (p - p_ref),
        nota = if (!calculable) "no calculable: la cota no supera el umbral o hay menos de 15 eventos"
               else if (!is.na(p_ref) && abs(p - p_ref) > 0.10) "corrimiento de más de 10 puntos" else "")
    }
  }
}
t1 <- do.call(rbind, f1); t1p <- do.call(rbind, p1)
write.csv(t1, "tables/robustez_umbral.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t1p, "tables/robustez_umbral_prob.csv", row.names = FALSE, fileEncoding = "UTF-8")
c1 <- t1[t1$serie == "uruguay_conjunto", ]
anotar("1 umbral y separación", "en los 9 casos el IC90 de b excluye el 0 y b está en 0,30–0,70",
       sprintf("b de %.3f a %.3f; IC90 más bajo %.3f; %d de 9 pasan", min(c1$b), max(c1$b), min(c1$b_boot_lo), sum(c1$pasa)),
       all(c1$pasa))
cat(sprintf("  probabilidades: corrimiento máximo %.1f puntos; %d de %d con más de 10 puntos; %d no calculables\n",
            max(abs(t1p$dif_puntos), na.rm = TRUE), sum(grepl("corrimiento", t1p$nota)), nrow(t1p),
            sum(grepl("no calculable", t1p$nota))))

# ---------- prueba 2: antes y después de las represas ----------------------------------------------------------
cat("== prueba 2: represas\n")
As <- lapply(URUGUAY, armar, pct = 98, gap = 10)
Ds <- lapply(As, datos_de); names(Ds) <- URUGUAY
x_ult <- tail(enso[!is.na(enso$x6), ], 1)$x6
f2 <- list(); p2 <- list()
for (corte in as.Date(c("2000-03-01", "2001-03-01"))) {
  corte <- as.Date(corte, origin = "1970-01-01")
  antes <- con_eventos(lapply(Ds, recortar, hasta = corte), 15)
  despues <- con_eventos(lapply(Ds, recortar, desde = corte), 15)
  ba <- b_conj(antes); bd <- b_conj(despues)
  bt <- boot_anios(Ds, function(Db) {
    a <- b_conj(con_eventos(lapply(Db, recortar, hasta = corte), 5))
    d <- b_conj(con_eventos(lapply(Db, recortar, desde = corte), 5))
    c(a, d, d - a)
  })
  f2[[length(f2) + 1]] <- data.frame(corte = as.character(corte), serie = "uruguay_conjunto",
    hidrometros_antes = length(antes), hidrometros_despues = length(despues),
    n_antes = sum(vapply(antes, function(D) length(D$xe), 0)), n_despues = sum(vapply(despues, function(D) length(D$xe), 0)),
    b_antes = ba, b_antes_lo = q5(bt[, 1]), b_antes_hi = q95(bt[, 1]),
    b_despues = bd, b_despues_lo = q5(bt[, 2]), b_despues_hi = q95(bt[, 2]),
    dif = bd - ba, dif_lo = q5(bt[, 3]), dif_hi = q95(bt[, 3]))
  cat(sprintf("  corte %s: antes b = %+.3f [%+.3f, %+.3f] (%d eventos); después b = %+.3f [%+.3f, %+.3f] (%d eventos); diferencia %+.3f [%+.3f, %+.3f]\n",
              corte, ba, q5(bt[, 1]), q95(bt[, 1]), sum(vapply(antes, function(D) length(D$xe), 0)),
              bd, q5(bt[, 2]), q95(bt[, 2]), sum(vapply(despues, function(D) length(D$xe), 0)),
              bd - ba, q5(bt[, 3]), q95(bt[, 3])))
  for (sid in URUGUAY) {
    Da <- recortar(Ds[[as.character(sid)]], hasta = corte); Dd <- recortar(Ds[[as.character(sid)]], desde = corte)
    f2[[length(f2) + 1]] <- data.frame(corte = as.character(corte), serie = as.character(sid),
      hidrometros_antes = NA, hidrometros_despues = NA, n_antes = length(Da$xe), n_despues = length(Dd$xe),
      b_antes = if (length(Da$xe) >= 15) b_conj(list(Da)) else NA, b_antes_lo = NA, b_antes_hi = NA,
      b_despues = if (length(Dd$xe) >= 15) b_conj(list(Dd)) else NA, b_despues_lo = NA, b_despues_hi = NA,
      dif = NA, dif_lo = NA, dif_hi = NA)
    z <- serie[as.character(sid), "evacuacion_m"]
    if (is.na(z) || length(Dd$xe) < 15 || z <= Dd$u) next
    f <- ajustar(Dd); f0 <- ajustar(Ds[[as.character(sid)]])
    for (x in c(0, 1.5, x_ult))
      p2[[length(p2) + 1]] <- data.frame(corte = as.character(corte), serie = sid, nombre = nombre_de(sid),
        nivel = "evacuación", z_m = z, x = x, n_despues = length(Dd$xe),
        p_despues = prob(z, Dd$u, f$a, f$b, f$c, f$d, f$xi, x),
        p_toda_la_serie = prob(z, Dd$u, f0$a, f0$b, f0$c, f0$d, f0$xi, x))
  }
}
t2 <- do.call(rbind, f2); t2p <- do.call(rbind, p2)
t2p$dif_puntos <- 100 * (t2p$p_despues - t2p$p_toda_la_serie)
write.csv(t2, "tables/robustez_represas.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t2p, "tables/robustez_represas_prob.csv", row.names = FALSE, fileEncoding = "UTF-8")
c2 <- t2[t2$serie == "uruguay_conjunto" & t2$corte == "2000-03-01", ]
anotar("2 represas (corte 1-mar-2000)", "el IC90 de b después del corte excluye el 0",
       sprintf("b después = %.3f [%.3f, %.3f]", c2$b_despues, c2$b_despues_lo, c2$b_despues_hi), c2$b_despues_lo > 0)
cat(sprintf("  ¿cambió la relación? diferencia %+.3f [%+.3f, %+.3f]: %s\n", c2$dif, c2$dif_lo, c2$dif_hi,
            if (c2$dif_lo > 0 || c2$dif_hi < 0) "el IC90 excluye el 0, cambió" else "el IC90 incluye el 0, no se distingue un cambio"))

# ---------- prueba 3: tendencia ----------------------------------------------------------------------------------
cat("== prueba 3: tendencia\n")
decadas <- function(f) (as.integer(format(f, "%Y")) + (as.integer(format(f, "%m")) - 0.5) / 12 - 1970) / 10
con_t <- function(D) { D$te <- decadas(D$ye); D$tm <- decadas(D$mm); D }
# λ = exp(a + b·x + g·t): para b y g dadas, la a de cada estación sale cerrada, igual que en el principal
nll_bg <- function(p, Ds) sum(vapply(Ds, function(D) {
  n <- length(D$xe); le <- p[1] * D$xm + p[2] * D$tm
  a <- log(max(n, 0.5)) - log(sum(exp(le)) / 12)
  -(sum(a + p[1] * D$xe + p[2] * D$te) - sum(exp(a + le)) / 12)
}, 0))
fit_bg <- function(Ds) optim(c(0.5, 0), nll_bg, Ds = Ds, method = "BFGS", control = list(reltol = 1e-12))$par
usados <- con_eventos(Ds, 15)
b0 <- b_conj(usados)                                             # control: sin tendencia tiene que dar la b del principal
bg <- fit_bg(lapply(usados, con_t))
bt <- boot_anios(usados, function(Db) fit_bg(lapply(con_eventos(Db, 5), con_t)))
ex <- enso[!is.na(enso$x6) & enso$mes >= min(do.call(c, lapply(usados, function(D) min(D$mm)))), ]
r_xt <- cor(ex$x6, decadas(ex$mes))
t3 <- data.frame(b_sin_tendencia = b0, b = bg[1], b_boot_lo = q5(bt[, 1]), b_boot_hi = q95(bt[, 1]),
                 g_por_decada = bg[2], g_boot_lo = q5(bt[, 2]), g_boot_hi = q95(bt[, 2]),
                 factor_por_decada = exp(bg[2]), cor_x6_tiempo = r_xt, meses = nrow(ex))
write.csv(t3, "tables/robustez_tendencia.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  sin tendencia b = %+.3f; con tendencia b = %+.3f [%+.3f, %+.3f], g = %+.4f por década [%+.4f, %+.4f]; cor(x6, tiempo) = %+.3f\n",
            b0, bg[1], q5(bt[, 1]), q95(bt[, 1]), bg[2], q5(bt[, 2]), q95(bt[, 2]), r_xt))
anotar("3 tendencia", "el IC90 de b con tendencia excluye el 0 y b está en 0,30–0,70",
       sprintf("b = %.3f [%.3f, %.3f]", bg[1], q5(bt[, 1]), q95(bt[, 1])),
       q5(bt[, 1]) > 0 && bg[1] >= B_LO && bg[1] <= B_HI)

# exploratorio, agregado después de ver la prueba 3 (log del 10-nov-2026): la misma tendencia sólo desde julio de
# 1965, después del salto de escala de San Javier, Garruchos y Santo Tomé. No entra en el veredicto.
desde65 <- as.Date("1965-07-01")
us65 <- con_eventos(lapply(Ds, recortar, desde = desde65), 15)
bg65 <- fit_bg(lapply(us65, con_t))
bt65 <- boot_anios(us65, function(Db) fit_bg(lapply(con_eventos(Db, 5), con_t)))
t3b <- data.frame(desde = as.character(desde65), n_eventos = sum(vapply(us65, function(D) length(D$xe), 0)),
                  b_sin_tendencia = b_conj(us65), b = bg65[1], b_boot_lo = q5(bt65[, 1]), b_boot_hi = q95(bt65[, 1]),
                  g_por_decada = bg65[2], g_boot_lo = q5(bt65[, 2]), g_boot_hi = q95(bt65[, 2]))
write.csv(t3b, "tables/robustez_tendencia_desde1965.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  exploratorio, desde %s (%d eventos): b = %+.3f [%+.3f, %+.3f], g = %+.4f por década [%+.4f, %+.4f]
",
            desde65, t3b$n_eventos, bg65[1], q5(bt65[, 1]), q95(bt65[, 1]), bg65[2], q5(bt65[, 2]), q95(bt65[, 2])))

# agregado por la lectura de Krepper et al. (2003, p. 114), con regla en el log del 10-nov-2026: ¿la b es la misma
# en 1930-1966, cuando ellos no ven respuesta del caudal a El Niño, que desde 1967? No entra en el veredicto.
k0 <- as.Date("1930-03-01"); k1 <- as.Date("1967-03-01")
tramo1 <- con_eventos(lapply(Ds, recortar, desde = k0, hasta = k1), 15)
tramo2 <- con_eventos(lapply(Ds, recortar, desde = k1), 15)
bk1 <- b_conj(tramo1); bk2 <- b_conj(tramo2)
btk <- boot_anios(Ds, function(Db) {
  a <- b_conj(con_eventos(lapply(Db, recortar, desde = k0, hasta = k1), 5))
  d <- b_conj(con_eventos(lapply(Db, recortar, desde = k1), 5))
  c(a, d, d - a)
})
tk <- data.frame(tramo1 = "1930-03 a 1967-02", tramo2 = "desde 1967-03",
  hidrometros_1 = length(tramo1), hidrometros_2 = length(tramo2),
  n_1 = sum(vapply(tramo1, function(D) length(D$xe), 0)), n_2 = sum(vapply(tramo2, function(D) length(D$xe), 0)),
  b_1 = bk1, b_1_lo = q5(btk[, 1]), b_1_hi = q95(btk[, 1]), b_2 = bk2, b_2_lo = q5(btk[, 2]), b_2_hi = q95(btk[, 2]),
  dif = bk2 - bk1, dif_lo = q5(btk[, 3]), dif_hi = q95(btk[, 3]))
write.csv(tk, "tables/robustez_1930_1966.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("== por Krepper (2003): 1930-1966 b = %+.3f [%+.3f, %+.3f] (%d hidrómetros, %d eventos); desde 1967 b = %+.3f [%+.3f, %+.3f] (%d, %d); diferencia %+.3f [%+.3f, %+.3f]
",
            bk1, q5(btk[, 1]), q95(btk[, 1]), length(tramo1), tk$n_1, bk2, q5(btk[, 2]), q95(btk[, 2]), length(tramo2), tk$n_2,
            bk2 - bk1, q5(btk[, 3]), q95(btk[, 3])))

# agregado por la lectura de Antico y Vuille (2022), con regla en el log del 10-nov-2026 (noche): ¿la b es la misma
# en 1941-1982, cuando ellos ven al Paraná desconectado de El Niño, que desde 1983? No entra en el veredicto.
a0 <- as.Date("1941-03-01"); a1 <- as.Date("1983-03-01")
tra1 <- con_eventos(lapply(Ds, recortar, desde = a0, hasta = a1), 15)
tra2 <- con_eventos(lapply(Ds, recortar, desde = a1), 15)
ba1 <- b_conj(tra1); ba2 <- b_conj(tra2)
bta <- boot_anios(Ds, function(Db) {
  a <- b_conj(con_eventos(lapply(Db, recortar, desde = a0, hasta = a1), 5))
  d <- b_conj(con_eventos(lapply(Db, recortar, desde = a1), 5))
  c(a, d, d - a)
})
ta <- data.frame(tramo1 = "1941-03 a 1983-02", tramo2 = "desde 1983-03",
  hidrometros_1 = length(tra1), hidrometros_2 = length(tra2),
  n_1 = sum(vapply(tra1, function(D) length(D$xe), 0)), n_2 = sum(vapply(tra2, function(D) length(D$xe), 0)),
  b_1 = ba1, b_1_lo = q5(bta[, 1]), b_1_hi = q95(bta[, 1]), b_2 = ba2, b_2_lo = q5(bta[, 2]), b_2_hi = q95(bta[, 2]),
  dif = ba2 - ba1, dif_lo = q5(bta[, 3]), dif_hi = q95(bta[, 3]))
write.csv(ta, "tables/robustez_1941_1982.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("== por Antico y Vuille (2022): 1941-1982 b = %+.3f [%+.3f, %+.3f] (%d hidrómetros, %d eventos); desde 1983 b = %+.3f [%+.3f, %+.3f] (%d, %d); diferencia %+.3f [%+.3f, %+.3f]\n",
            ba1, q5(bta[, 1]), q95(bta[, 1]), length(tra1), ta$n_1, ba2, q5(bta[, 2]), q95(bta[, 2]), length(tra2), ta$n_2,
            ba2 - ba1, q5(bta[, 3]), q95(bta[, 3])))

v <- do.call(rbind, veredicto)
write.csv(v, "tables/robustez_veredicto.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("== puerta de robustez: %s (%d de %d pruebas pasan)\n", if (all(v$pasa)) "PASA" else "NO PASA", sum(v$pasa), nrow(v)))
