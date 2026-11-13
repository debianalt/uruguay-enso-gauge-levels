# diagnóstico de la falla de calibración de Garruchos (cota de alerta).
# Regla en el registro de decisiones del autor (11-nov-2026), escrita antes de correr: G1 agrupamiento, G2 cola de la GPD,
# G3 homogeneidad en aguas altas, G4 ajuste posterior al quiebre de 1965.
# Uso (desde la raíz del repositorio):  Rscript scripts/calibracion_garruchos.R

source("scripts/modelo.R", encoding = "UTF-8")
set.seed(2026); N_BOOT <- 500L; N_MIN <- 10L; CORTE <- as.Date("1965-07-01")
tabla <- read.csv("tables/prob_frecuencia_sola.csv", fileEncoding = "UTF-8")
S_gpd <- function(z, u, c0, xi0) {
  s <- exp(c0); t <- 1 + xi0 * (z - u) / s
  if (abs(xi0) > 1e-6) ifelse(t > 0, pmax(t, 0)^(-1 / xi0), 0) else exp(-(z - u) / s)
}

# calibración de un ajuste f sobre los datos D: años completos; Poisson y binomial negativa; año y crecida
calibrar <- function(D, f, m0, cotas) {
  am <- anio_h(D$mm, m0); ae <- anio_h(D$ye, m0)
  comp <- as.integer(names(which(table(am) == 12)))
  L <- tapply(exp(f$a + f$b * D$xm) / 12, am, sum)[as.character(comp)]
  O_n <- sapply(comp, function(a) sum(ae == a)); n <- length(comp)
  phi <- sum((O_n - L)^2 / L) / (n - 2)
  ex <- sum((O_n - L)^2) - sum(L); k <- if (ex > 0) sum(L^2) / ex else Inf
  do.call(rbind, lapply(names(cotas), function(niv) {
    z <- cotas[[niv]]; S <- S_gpd(z, D$u, f$c0, f$xi0)
    obs <- sapply(comp, function(a) any(D$y[ae == a] + D$u >= z))
    Pp <- 1 - exp(-L * S); Pn <- if (is.finite(k)) 1 - (1 + L * S / k)^(-k) else Pp
    zeta <- function(P) (sum(obs) - sum(P)) / sqrt(sum(P * (1 - P)))
    Oev <- sum(sapply(comp, function(a) sum(D$y[ae == a] + D$u >= z))); Eev <- sum(L * S)
    data.frame(nivel = niv, z_m = z, anios = n, phi = phi, k = k, obs_anios = sum(obs),
      esp_poisson = sum(Pp), zeta_poisson = zeta(Pp), esp_nb = sum(Pn), zeta_nb = zeta(Pn),
      obs_crecidas = Oev, esp_crecidas = Eev, zeta_crecidas = (Oev - Eev) / sqrt(max(phi, 1) * Eev))
  }))
}

# ---- G1 y G2: hidrómetros con cota y al menos diez años completos
g12 <- list()
for (sid in unique(tabla$serie)) {
  m0 <- serie[as.character(sid), "mes_inicio_hidro"]; D <- datos(sid, "x6"); f <- ajustar(D)
  if (sum(table(anio_h(D$mm, m0)) == 12) < N_MIN) next
  t <- unique(tabla[tabla$serie == sid, c("nivel", "z_m")]); cotas <- setNames(as.list(t$z_m), t$nivel)
  g12[[length(g12) + 1]] <- cbind(serie = sid, nombre = serie[as.character(sid), "nombre"], tramo = "toda la serie",
                                  calibrar(D, f, m0, cotas))
}
g12 <- do.call(rbind, g12)

# ---- G4: Garruchos desde julio de 1965, mismo u
D4 <- datos(67, "x6", desde = CORTE); f4 <- ajustar(D4)
g4 <- cbind(serie = 67, nombre = "Garruchos", tramo = "desde julio de 1965",
            calibrar(D4, f4, 3, list(alerta = 12, "evacuación" = 14)))
xs <- sort(unique(tabla$x))
p4 <- expand.grid(nivel = c("alerta", "evacuación"), x = xs, stringsAsFactors = FALSE)
p4$z_m <- ifelse(p4$nivel == "alerta", 12, 14)
p4$p_post1965 <- prob(p4$z_m, D4$u, f4$a, f4$b, f4$c0, 0, f4$xi0, p4$x)
p4 <- merge(p4, tabla[tabla$serie == 67, c("nivel", "x", "p_frec")], by = c("nivel", "x"))
p4$dif_puntos <- 100 * (p4$p_post1965 - p4$p_frec)

# ---- G3: tasa después/antes del corte, remuestreo conjunto de años hidrológicos
G3 <- c(65, 67, 68); DD <- lapply(G3, function(s) datos(s, "x6")); names(DD) <- G3
razon <- function(D, w) {           # w: peso de cada año hidrológico (veces que salió en el sorteo)
  we <- w[as.character(anio_h(D$ye, 3))]; we[is.na(we)] <- 0
  wm <- w[as.character(anio_h(D$mm, 3))]; wm[is.na(wm)] <- 0
  (sum(we[D$ye >= CORTE]) / (sum(wm[D$mm >= CORTE]) / 12)) / (sum(we[D$ye < CORTE]) / (sum(wm[D$mm < CORTE]) / 12))
}
todos <- sort(unique(unlist(lapply(DD, function(D) anio_h(D$mm, 3)))))
est <- sapply(DD, razon, w = setNames(rep(1, length(todos)), todos))
bt <- t(replicate(N_BOOT, { w <- table(sample(todos, length(todos), replace = TRUE))
  sapply(DD, razon, w = setNames(as.numeric(w), names(w))) }))
q <- function(v) quantile(v, c(0.05, 0.95), na.rm = TRUE)
g3 <- data.frame(medida = c(paste("tasa después/antes,", serie[as.character(G3), "nombre"]),
                            "razón de razones Garruchos/San Javier", "razón de razones Garruchos/Santo Tomé"),
  estimacion = c(est, est[["67"]] / est[["65"]], est[["67"]] / est[["68"]]),
  ic_lo = c(apply(bt, 2, q)[1, ], q(bt[, "67"] / bt[, "65"])[1], q(bt[, "67"] / bt[, "68"])[1]),
  ic_hi = c(apply(bt, 2, q)[2, ], q(bt[, "67"] / bt[, "65"])[2], q(bt[, "67"] / bt[, "68"])[2]))

diag <- rbind(g12, g4)
write.csv(diag, "tables/calibracion_diagnostico.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(g3, "tables/garruchos_homogeneidad.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(p4, "tables/garruchos_post1965_prob.csv", row.names = FALSE, fileEncoding = "UTF-8")
options(width = 220)
cat("--- G1 (nb) y G2 (crecidas); G4 en la última fila de Garruchos\n"); print(format(diag[, -1], digits = 3), row.names = FALSE)
cat(sprintf("\nG4: b desde julio de 1965 = %.3f; xi0 = %.3f; crecidas %d\n", f4$b, f4$xi0, length(D4$xe)))
cat("\n--- G3\n"); print(format(g3, digits = 3), row.names = FALSE)
cat("\n--- G4, probabilidades (frecuencia sola)\n"); print(format(p4, digits = 3), row.names = FALSE)
