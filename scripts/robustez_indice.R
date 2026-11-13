# la b del conjunto con otras formas del índice, y los intervalos perfilados al 90 %.
#
# Regla en el registro de decisiones del autor (12-nov-2026), escrita antes de correr.
#   prueba 6, rezago   el índice de cada mes es la media de seis meses que termina tres meses antes (x6 de t − 3).
#                      Pasa si el IC90 de la b del conjunto excluye el cero. No se prueban otros rezagos.
#   descriptivo        b del conjunto con ventana de 3 y de 12 meses, y con ERSST v6 desde julio de 1948 (contra
#                      HadISST en el mismo tramo).
#   perfiles al 90 %   b y d por hidrómetro (χ² de 1 gl = 2,71), y la p de la prueba conjunta de b y d.
#
# No toca robustez.R ni sus tablas. b_conj y boot_anios son los de robustez.R, copiados sin cambios.
# Uso (desde la raíz del repositorio):  Rscript scripts/robustez_indice.R [--boot 500]

args <- commandArgs(trailingOnly = TRUE)
N_BOOT <- if ("--boot" %in% args) as.integer(args[which(args == "--boot") + 1]) else 500L
source("scripts/modelo.R", encoding = "UTF-8")
set.seed(2026)
M0 <- 3
q5 <- function(v) unname(quantile(v, 0.05)); q95 <- function(v) unname(quantile(v, 0.95))

enso$x6_lag3 <- c(rep(NA, 3), head(enso$x6, -3))

b_conj <- function(Ds) {
  if (!length(Ds)) return(NA_real_)
  optimize(function(b) sum(vapply(Ds, function(D) pnll_pois(b, D$xe, D$xm), 0)), c(-6, 6), tol = 1e-10)$minimum
}
con_eventos <- function(Ds, k) Ds[vapply(Ds, function(D) length(D$xe) >= k, TRUE)]
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

conj <- function(caso, xcol, desde = NULL) {
  Ds <- con_eventos(lapply(URUGUAY, function(s) datos(s, xcol, desde)), 15)
  b <- b_conj(Ds)
  bt <- boot_anios(Ds, function(Db) b_conj(con_eventos(Db, 5)))
  r <- data.frame(caso = caso, indice = xcol, desde = if (is.null(desde)) "" else as.character(desde),
                  hidrometros = length(Ds), n_eventos = sum(vapply(Ds, function(D) length(D$xe), 0)),
                  b = b, b_boot_lo = q5(bt), b_boot_hi = q95(bt))
  cat(sprintf("  %-34s %d hidrómetros, %4d eventos, b = %+.3f [%+.3f, %+.3f]\n", caso, r$hidrometros, r$n_eventos,
              b, r$b_boot_lo, r$b_boot_hi))
  r
}

cat("== b del conjunto según el índice\n")
E48 <- as.Date("1948-07-01")
t <- rbind(
  conj("principal (control)", "x6"),
  conj("rezago de tres meses", "x6_lag3"),
  conj("ventana de tres meses", "x3"),
  conj("ventana de doce meses", "x12"),
  conj("HadISST desde 1948-07", "x6", E48),
  conj("ERSST v6 desde 1948-07", "x6_ersst", E48))
write.csv(t, "tables/robustez_indice.csv", row.names = FALSE, fileEncoding = "UTF-8")
if (abs(t$b[1] - 0.51702835346803) > 1e-6) { cat("el control no reproduce la b del principal: se para acá\n"); quit(status = 1) }

rz <- t[t$indice == "x6_lag3", ]
ex <- enso[!is.na(enso$x6) & !is.na(enso$x6_lag3), ]
v <- data.frame(prueba = "6 rezago de tres meses", regla = "el IC90 de la b del conjunto con el índice rezagado excluye el 0",
                valor = sprintf("b = %.3f [%.3f, %.3f]; cor(x6, x6 rezagado) = %.3f", rz$b, rz$b_boot_lo, rz$b_boot_hi,
                                cor(ex$x6, ex$x6_lag3)),
                pasa = rz$b_boot_lo > 0)
write.csv(v, "tables/robustez_rezago.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("== prueba 6: %s | %s\n", v$valor, if (v$pasa) "PASA" else "NO PASA"))

cat("== perfiles al 90 % por hidrómetro\n")
aj <- read.csv("tables/ajuste.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
aj <- aj[aj$variante == "principal" & aj$modelo == "pot", ]
filas <- list()
for (sid in URUGUAY) {
  D <- datos(sid, "x6"); f <- ajustar(D)
  CHI1 <- qchisq(0.95, 1)
  ib95 <- perfil(function(b) pnll_pois(b, D$xe, D$xm), f$b, f$nll_p)
  CHI1 <- qchisq(0.90, 1)
  ib <- perfil(function(b) pnll_pois(b, D$xe, D$xm), f$b, f$nll_p)
  id <- perfil(function(d) fit_gpd_d(d, D$y, D$xe)$nll, f$d, f$nll_g)
  r0 <- aj[aj$serie == as.character(sid), ]
  ok <- abs(f$b - r0$b) < 1e-6 && abs(ib95[1] - r0$b_lo) < 1e-4 && abs(ib95[2] - r0$b_hi) < 1e-4
  filas[[length(filas) + 1]] <- data.frame(serie = sid, nombre = serie[as.character(sid), "nombre"],
    n_eventos = length(D$xe), b = f$b, b_lo90 = ib[1], b_hi90 = ib[2], d = f$d, d_lo90 = id[1], d_hi90 = id[2],
    p_conjunta = pchisq(2 * (f$l1 - f$l0), 2, lower.tail = FALSE), lodo_pos = r0$lodo_pos, lodo_n = r0$lodo_n,
    veredicto = r0$veredicto, reproduce_ajuste = ok)
  cat(sprintf("  %d %-12s b = %.2f (%.2f a %.2f); d = %+.2f (%+.2f a %+.2f); p = %.4f; %s%s\n", sid,
              serie[as.character(sid), "nombre"], f$b, ib[1], ib[2], f$d, id[1], id[2],
              filas[[length(filas)]]$p_conjunta, r0$veredicto, if (ok) "" else "  NO REPRODUCE ajuste.csv"))
}
write.csv(do.call(rbind, filas), "tables/ajuste_perfil90.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("== última lectura por hidrómetro\n")
diario <- read.csv("data/crecidas_diario.csv", stringsAsFactors = FALSE)
print(aggregate(fecha ~ serie, diario[!is.na(diario$pref), ], max))
