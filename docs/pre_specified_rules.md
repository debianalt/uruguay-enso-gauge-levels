# Pre-specified decision rules

Each block below is copied verbatim, in the original Spanish, from the author's dated working log. Every block was
written before the analysis it governs was run. Result paragraphs that the log added after each run are left out;
where a block itself records something added after seeing a result, it says so in its own words.

The rules cover: the R reproduction of the model (9 Nov 2026); the robustness tests for flood definition, dams and
trend (10 Nov); the constancy test after Krepper et al. (2003) and after Antico and Vuille (2022) (10 Nov); the table
with the height effect fixed at zero (10 Nov); calibration and the Poisson check (11 Nov); the Garruchos diagnosis
(11 Nov); the lagged index (12 Nov); and the seasonal cycle (13 Nov).

## 9 nov 2026: reglas de la reproducción en R, escritas ANTES de correr

- **Qué se reproduce**: el modelo principal de `neahub/pipeline/crecidas/frecuencia.py` (docstring, commit 1cdbca3),
  sin cambios:
  - tasa λ(x) = exp(a + b·x);
  - excedente GPD con σ(x) = exp(c + d·x) y ξ constante, acotada en [−0,5, 1];
  - x = x6 del Niño 3.4;
  - eventos y meses válidos tal como están en `crecidas.evento` y `crecidas.mes`;
  - el Uruguay en conjunto con a, c y ξ propios de cada estación y b y d compartidos;
  - probabilidad anual 1 − exp(−λ·S).

  No se cambia la definición de los eventos: se leen hechos.
- **Tolerancia** contra `crecidas.ajuste` y `crecidas.prob_condicional`:
  - b y d puntuales: ±0,02;
  - probabilidades puntuales: ±2 puntos;
  - intervalos (perfilados y bootstrap): ±0,05 en b y d, ±4 puntos en las probabilidades, porque el bootstrap
    depende de la semilla.

  Lo que se pase se investiga y se escribe aquí; no se ajusta el R para que coincida sin entender por qué difiere.
- **Prueba sintética**, la misma del Python: 200 réplicas con λ = exp(0 + 0,3·x), σ = exp(log 2 + 0,15·x), ξ = 0,1
  y la x6 real de 1950–2009. Pasa si el intervalo perfilado del 95 % cubre la b y la d verdaderas en ≥ 90 % de las
  réplicas.
- **Bootstrap**: 500 réplicas, años hidrológicos enteros en todas las estaciones a la vez, semilla 2026.
- **Más tolerancias**, agregadas antes de correr:
  - GEV de control: μ1 ±0,05 m/°C;
  - los valores *p* se informan sin juzgarlos, pero los veredictos tienen que coincidir uno a uno;
  - el número de eventos y de décadas evaluadas tiene que ser idéntico.
- **Cómo se implementa, sin cambiar el modelo.**
  - Para una b dada, la a de máxima verosimilitud tiene fórmula cerrada (exp(a) = n / Σ exp(b·x_m)/12). R
    perfila a y optimiza b en una dimensión; el Python optimizaba las dos con Nelder-Mead.
  - En el conjunto, d se perfila: para cada d, la c y la ξ de cada estación se ajustan por separado. El Python usaba
    Powell sobre 15 parámetros a la vez. El máximo es el mismo; el camino es más fácil de leer y más estable.
  - La GEV se escribe a mano: no se instala extRemes. Se valida contra la GEV de Python, que ya se había comparado
    con scipy.
- **Alcance de la reproducción**: los siete hidrómetros del Uruguay y el conjunto, con todas las variantes
  (principal, +1 m, ERSST v6, después del quiebre, ventanas de 3 y 12 meses). El Paraná y el Iguazú no se
  reproducen: no entran en la comunicación.
- **Exportación**: `scripts/exportar_crecidas.py` saca a `data/` las tablas del esquema `crecidas` tal como las dejó
  el Python, `crecidas.homogeneidad` incluida.

## 10 nov 2026: la puerta antes de escribir. Reglas de la robustez, escritas ANTES de correr
- **Decisión del autor**: «procede con tu recomendación». La recomendación era no eliminar la carpeta y no
  redactar todavía: correr primero la robustez pendiente y hacer las dos lecturas (Bischoff 2000, Krepper 2003).
  Si la puerta pasa, se escribe; si no, se estaciona con el motivo acá.
- **Qué se prueba**: la b del conjunto, que es lo único que se afirma. b sale sólo de la parte de Poisson del
  modelo, así que el bootstrap de estas pruebas reajusta la tasa y no la GPD. Semilla 2026, 500 réplicas, años
  hidrológicos enteros con todas las estaciones a la vez. Los intervalos son del 90 %, percentiles 5 y 95.
- **Paso previo, de control**: rearmar en R los eventos y los meses válidos desde la serie diaria
  (`crecidas.diario`, columna `pref`) con la regla de `neahub/pipeline/crecidas/build_crecidas.py` (P98 de todas
  las lecturas, evento nuevo cuando pasan más de 10 días sin lectura sobre u, mes válido con ≥ 80 % de los días).
  Con P98 y 10 días tienen que salir los de `data/crecidas_evento.csv` y `crecidas_mes.csv`, uno a uno. Si no
  salen, se para.
- **Prueba 1, umbral y separación**: P97, P98, P99 × 5, 10, 20 días (9 casos).
  - Pasa si en los 9 el IC90 de b excluye el cero y la b puntual queda dentro del IC90 del principal (0,30–0,70).
  - Probabilidad puntual de alerta y evacuación por pueblo en neutro y en Niño +1,5: se informa. Un corrimiento
    de más de 10 puntos contra el principal se marca.
  - Una cota por debajo del umbral se informa como no calculable.
- **Prueba 2, represas**: corte al empezar el año hidrológico 2000 (1-mar-2000), por el llenado del embalse de
  Itá. <!-- VERIFICAR: las fuentes encontradas no coinciden (llenado en 2000 según Wikipedia en portugués, sin
  fuente; «entró en operación en 2001» según un estudio de 2023 en Novos Cadernos NAEA). Falta ANEEL u ONS. -->
  Por eso se corre también con corte al 1-mar-2001, como control que se informa y no decide.
  - Se ajusta el conjunto antes y después del corte y se sortea b_después − b_antes.
  - Pasa si el IC90 de b_después excluye el cero.
  - Si el IC90 de la diferencia excluye el cero, la relación cambió con la regulación: se informa, y sacar la
    tabla del período regulado es un cambio de análisis que decide el autor.
  - Probabilidad de evacuación por pueblo con el ajuste posterior al corte, al lado de la de toda la serie: se
    informa.
- **Prueba 3, tendencia**: λ = exp(a + b·x + g·t), t en décadas desde 1970, b y g compartidas y a de cada
  estación. Se informa g con su IC90 y la correlación entre x6 y el tiempo.
  - Pasa si el IC90 de b excluye el cero y la b puntual queda dentro de 0,30–0,70.
- **Veredicto**: la puerta de robustez pasa si pasan las tres. Las reglas no se tocan después de ver los números.
- **Código**: las funciones de `scripts/frecuencia.R` pasan a `scripts/modelo.R` para que `scripts/robustez.R` las
  use. Control: `tables/ajuste.csv` y `tables/prob_condicional.csv` tienen que salir idénticos después del cambio.
- **Agregado después de ver los resultados (exploratorio, no entra en el veredicto)**: la prueba 3 dio una
  tendencia con IC90 que excluye el cero. San Javier, Garruchos y Santo Tomé tienen un salto de escala en junio
  de 1965, que por sí solo puede fabricar una tendencia en las superaciones de un umbral fijo. Se corre la misma
  prueba 3 sólo con datos desde el 1-jul-1965. Se informa b y g con su IC90; no hay regla de pasa o no pasa.

## 10 nov 2026: lectura de Bischoff (2000) y Krepper (2003), y una prueba que sale de ahí
- El autor bajó los dos PDF; están en la biblioteca (`Bischoff2000-16`, `Krepper2003-18`). Leídos completos.
- **Regla escrita antes de correr, motivada por la lectura.** Krepper et al. (2003, p. 114) dicen que el modo de
  unos 3 años del caudal, el que asocian a El Niño, se atenúa entre 1940 y 1960: entre 1930 y 1966 hubo siete
  eventos de El Niño «without an apparent response» en esa oscilación, y el 76 % de las amplitudes grandes cae
  después de 1967. El modelo de la comunicación supone una b constante de 1907 a 2026. Se prueba:
  - b del conjunto en 1-mar-1930 a 28-feb-1967 contra b desde el 1-mar-1967, con los hidrómetros que tengan 15
    eventos o más en cada tramo, y la diferencia sorteada por años (semilla y réplicas como en el resto).
  - Si el IC90 de la diferencia incluye el cero: no se distingue un cambio de b entre los dos tramos, y se dice
    así, citando a Krepper.
  - Si lo excluye: la b no es constante. La tabla tendría que salir del tramo posterior a 1967 y eso lo decide
    el autor.
  - No entra en el veredicto de la puerta del 10-nov, que ya estaba cerrado; se informa aparte.
- **Qué hace Bischoff et al. (2000)** (*Water International* 25(3), 446–456; cinco autores: Bischoff, García,
  Vargas, Jones y Conway, confirmado en el artículo).
  - Caudales mensuales (m³/s) de 1909–1988 en Santo Tomé, Paso de los Libres, Monte Caseros y Concordia (Tabla
    1, p. 448). Separa los años en fase cálida, fría y neutra y compara el caudal medio mensual de cada fase
    (Fig. 8, p. 452).
  - Resultado: con fase cálida el caudal tiende a estar por encima de la media desde la primavera hasta marzo
    (p. 452); las medias de noviembre y diciembre difieren entre fase cálida y fría al 1 % (p. 452).
  - El resumen habla de «maximum floods», pero el cuerpo trabaja el caudal medio mensual por fase. No cuenta
    crecidas, no usa alturas ni cotas, no da probabilidades y no llega aguas arriba de Santo Tomé.
- **Qué hace Krepper et al. (2003)** (*International Journal of Climatology* 23, 103–115).
  - Análisis espectral singular del caudal medio anual en Itá, Santo Tomé, Paso de los Libres y Concordia
    (Tabla I, p. 105) y de la lluvia por subcuenca. Dos modos, de unos 6 y unos 3,5 años, explican cerca de un
    tercio de la variabilidad (Tabla II, p. 109; p. 114). Asocian a El Niño el de 3,5 años.
  - Ese modo está modulado: atenuado entre 1940 y 1960 y amplificado después de 1960 (p. 114).
  - Tendencia positiva del caudal desde 1950, con cambio hacia 1970 (p. 106, Fig. 2).
  - No cuenta crecidas, no usa cotas, no da probabilidades por localidad.
- **Novedad: confirmada.** Ninguno de los dos, ni Vargas et al. (2010), da la probabilidad anual de superar las
  cotas oficiales de cada pueblo según El Niño. La afirmación 10 deja de estar en suspenso. La nota de novedad
  del 9-nov atribuía a Bischoff «las crecidas máximas»: se corrige a caudal medio mensual por fase.

## 10-nov-2026 (noche): Antico y Vuille (2022) leído; regla de una prueba más, escrita antes de correr

- **Antico y Vuille (2022)**, *International Journal of Climatology* 42, 7269–7279, leído completo.
  - Caudal medio anual del Paraná en Rosario, 1876–2016, año hidrológico de septiembre a agosto (p. 7271).
    Correlación móvil de 31 años con Niño 3.4 y SOI.
  - Conectados en 1876–1940 y 1983–2016; no en 1941–1982 (Fig. 3b, p. 7272; conclusión 1, p. 7276).
  - Lo atribuyen a un Niño de poca amplitud hacia 1960 y a la fase más fría de la IPO (pp. 7275–7276).
  - Las crecidas más graves del Paraná desde 1875 (1878, 1905, 1983, 1992, 1998) caen todas fuera de
    1941–1982 (p. 7274).
  - Es el Paraná y es caudal medio anual; no trata el Uruguay ni cuenta crecidas.
- **Regla de la prueba (exploratoria, no entra en el veredicto de la puerta).** Mismo método que la de Krepper:
  b del conjunto en el tramo 1941-03 a 1983-02 contra desde 1983-03, hidrómetros con 15 eventos o más en el
  tramo, bootstrap por años hidrológicos, semilla y réplicas del script. El año hidrológico de acá empieza en
  marzo y el de ellos en septiembre: los cortes son los más cercanos, no idénticos.
  - Se dice «no se distingue un cambio» si el IC90 de la diferencia incluye el cero.
  - Se dice «en 1941–1982 la relación no se distingue de cero» si el IC90 de b de ese tramo incluye el cero. En
    ese caso el texto lo declara como límite y la tabla de probabilidades se acompaña de la del tramo desde 1983.
  - Las tablas ya escritas (`robustez_*.csv`) tienen que salir idénticas; se compara por suma de control.

## 10-nov-2026 (noche): auditoría de consistencia del borrador; regla de una tabla más, escrita antes de correr

- **Regla de la tabla «frecuencia sola»** (no cambia el veredicto del 9-nov ni la puerta del 10-nov):
  - mismos hidrómetros, cotas y valores del índice que `tables/prob_condicional.csv`;
  - tasa con la a y la b del ajuste de cada hidrómetro (la parte de Poisson no depende de la GPD);
  - altura con la GPD estacionaria del mismo ajuste (c0 y ξ0, d = 0);
  - intervalo del 90 % con el mismo remuestreo por años hidrológicos, 500 réplicas, semilla 2026;
  - sale a `tables/prob_frecuencia_sola.csv`, con la diferencia en puntos contra la tabla principal.
  - **Qué se hace con el resultado:** si en evacuación a +1,5 °C la diferencia pasa de 5 puntos en algún
    hidrómetro, el texto da las dos versiones y dice cuál corresponde a lo que se afirma. Cuál va como Tabla 2
    lo decide el autor; hasta entonces el borrador muestra las dos.

## 11-nov-2026: decisiones del autor (8, 6 y 9) y regla de calibración y Poisson, escrita antes de correr

- **Decisiones del autor** (por diálogo, 11-nov): (8) la tabla principal es la de frecuencia sola y la del modelo
  completo va al suplemento; (6) el argumento de 4.1 se mantiene, anclado en esa tabla (con d = 0 la probabilidad
  de alcanzar una cota alta igual sube con el índice: sale de la fórmula de 2.5, no pide análisis nuevo); (9) se
  corre calibración y supuesto de Poisson. Los otros cuatro análisis de la auditoría (estacionalidad, asimetría,
  b común, trayectoria mensual en las tablas) no se hacen y quedan como límites declarados.
- **Regla de calibración y Poisson** (aprobada por el autor con el plan del 11-nov, antes de correr):
  - **Datos**: por hidrómetro, años hidrológicos (marzo–febrero) con los doce meses válidos. Tasa esperada del
    año con la trayectoria mensual observada del índice: Λ_t = Σ_m exp(a + b·x_m)/12, con a y b del ajuste
    principal.
  - **A. Calibración** (seis hidrómetros con cota; alerta y evacuación; sin la alerta de Alba Posse). Observado:
    el máximo diario del año alcanza la cota. Predicho con frecuencia sola: P_t = 1 − exp(−Λ_t·S(z − u)). Se
    compara O = Σ observados con E = ΣP_t, varianza V = ΣP_t(1 − P_t). **Calibrado si |O − E|/√V < 1,645.** Se
    informa también por tercil del índice medio del año, sin veredicto. Es dentro de la muestra y se dice así.
  - **B. Poisson** (siete hidrómetros). Dispersión de Pearson de los conteos anuales, φ = Σ(O_t − Λ_t)²/Λ_t
    /(n − 2), con valor p por simulación de 2.000 series Poisson(Λ_t), semilla 2026. **Se rechaza el supuesto en
    un hidrómetro si p < 0,05.**
  - **Qué se hace con el resultado**: si A pasa en todas las cotas de los cuatro hidrómetros que sostienen la
    relación y tienen cota de evacuación (El Soberbio, San Javier, Garruchos, Santo Tomé) y B no rechaza en
    ellos, 4.5 cambia «ese supuesto no se puso a prueba» por el resultado. Si B rechaza: probabilidad con conteo
    binomial negativo, P = 1 − (1 + λS/k)^(−k), con k de la dispersión, y diferencia en puntos contra la tabla
    principal; hasta 5 puntos en evacuación a +1,5 °C va como límite cuantificado en 4.5. **Más de 5 puntos, o A
    falla: se para y se reporta al autor antes de tocar el texto.**

## 11-nov-2026: Garruchos, regla para diagnosticar la falla de calibración (escrita antes de correr)

El autor pidió seguir con el plan y abordar la falla de Garruchos «de manera científica». Dato previo a la regla,
descriptivo: el salto de 1965 en aguas bajas deja a Garruchos 0,5–0,6 m más bajo respecto de sus dos vecinos, pero
su tasa de crecidas pasa de 19 en 34,4 años válidos (0,55 por año) a 114 en 60,1 (1,9 por año). Un corrimiento del
cero hacia abajo daría menos crecidas, no más, así que «cambio de cero» no alcanza como explicación. Tres fuentes
candidatas y un ajuste de control, con el mismo criterio |ζ| < 1,645 (`scripts/calibracion_garruchos.R`; corre en
todas las cotas, el veredicto es sobre Garruchos):
- **G1, agrupamiento**: la calibración A con P_t binomial negativa, P_t = 1 − (1 + Λ_t·S/k)^(−k), con la k de B.
- **G2, cola de la GPD**: a nivel de crecida y no de año. O = crecidas con pico ≥ z en los años completos;
  E = Σ_t Λ_t·S; ζ = (O − E)/√(φ·E), con la φ de B (1 si es menor).
- **G3, homogeneidad en aguas altas**: tasa de crecidas después/antes de julio de 1965 en San Javier, Garruchos y
  Santo Tomé, y razón de razones Garruchos/San Javier y Garruchos/Santo Tomé, con IC90 por remuestreo conjunto de
  años hidrológicos (500, semilla 2026). Si las dos razones de razones excluyen el 1, Garruchos tiene una
  inhomogeneidad propia en aguas altas; si no, el cambio es regional (tendencia; Barros et al. 2015).
- **G4, ajuste posterior al quiebre**: Garruchos ajustado desde julio de 1965 (mismo u), calibrado en los años
  completos de ese tramo con Poisson y con binomial negativa.

Qué se hace con el resultado:
- G1 pasa: la falla se atribuye al agrupamiento, que el texto ya informa como límite cuantificado; la fila no
  lleva nota aparte.
- G1 falla y G3 da inhomogeneidad propia: la fila de Garruchos lleva nota (serie no homogénea en aguas altas, con
  O y E). Si G4 pasa, se informan sus probabilidades; reemplazar la fila por el ajuste posterior a 1965 lo decide el
  autor.
- G2 falla: el texto dice que la GPD sobrestima la cola en esa cota.
- En todos los casos 4.5 da O y E de la cota de alerta de Garruchos.

### 12-nov-2026 — evaluación adversarial y regla del rezago (escrita antes de correr)

**Regla de la prueba 6, rezago (antes de correr).** El artículo presenta el índice como disparador, pero el modelo usa
el índice simultáneo al pico. Se repite el ajuste del conjunto (los 657 eventos del principal, mismos meses válidos)
con el índice rezagado: a cada mes se le asigna la media de seis meses que termina **tres meses antes** (`x6` del mes
t − 3), para los eventos y para la exposición. Mismo remuestreo de años hidrológicos, 500 réplicas, semilla 2026.
- **Pasa** si el IC90 de la *b* del conjunto con el índice rezagado excluye el cero.
- Si pasa: el texto dice que la relación se sostiene con el índice conocido tres meses antes, con esa *b* y su
  intervalo, y el reclamo del disparador queda. Si no pasa: se eliminan del resumen y de 4.4 las oraciones del
  disparador.
- No se prueban otros rezagos después de ver el resultado. Las probabilidades de la Tabla II no se recalculan con el
  índice rezagado.

## 13-nov-2026 (noche): prueba 7, estacionalidad — regla escrita antes de correr
Objeción (evaluación GAN, revisor independiente): la tasa λ(x) = exp(a + b·x) no tiene ciclo anual; las crecidas del
Uruguay y El Niño tienen estación, y la *b* podría estar tomando parte de la estación. El autor pidió correr la prueba.
- Modelo: λ(x, m) = exp(a_g + b·x + Σ_{k=1,2} [s_k·sin(2πkm/12) + c_k·cos(2πkm/12)]), m = mes calendario del pico (eventos)
  o del mes válido (exposición). *a* por hidrómetro (cerrada dada el resto, como `a_hat`), *b* y los cuatro coeficientes
  armónicos comunes a los hidrómetros con 15 eventos o más (como `b_conj`). Índice x6, serie completa.
- Intervalo: bootstrap por años hidrológicos, 500 réplicas, semilla 2026, hidrómetros con 5 eventos o más en cada
  réplica (como `boot_anios` de `scripts/robustez_indice.R`).
- Pasa si el IC90 de *b* excluye el cero. Se informa además si la *b* cae dentro del IC90 principal (0,30–0,70) y la
  razón de verosimilitudes de los armónicos (descriptiva).
- Puerta: si la *b* sale del IC90 principal o pierde el signo, se para y se lleva al autor antes de tocar el texto.
- No se recalculan la Tabla II ni las probabilidades con el modelo estacional (la probabilidad anual integra la
  estación); la prueba es sólo sobre *b*. Script `scripts/robustez_estacion.R`, salida `tables/robustez_estacion.csv`.
