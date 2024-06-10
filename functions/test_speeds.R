test_speeds <- function(id, trace = TRUE, cores = 1) {
  tel <- as.telemetry(d$tel[[which(d$animal == id)]], mark.rm = TRUE)
  guess <- ctmm.guess(tel, interactive = FALSE)
  m_0 <- ctmm.select(tel, CTMM = guess, trace = trace, cores = cores)
  summary(m_0, units = FALSE)$CI['speed (meters/second)', ]
}
