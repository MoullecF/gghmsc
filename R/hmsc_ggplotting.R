#' Plot model convergence indexes
#'
#' @param beta if TRUE, plots the beta (env. filters) parameters
#' @param V if TRUE, plots the V parameters
#' @param gamma if TRUE, plots the gamma (traits) paramters
#' @param omega if TRUE, plots the omega (spp associations) parameters
#' @param title character string to customize
#' @export
ggplot_convergence <- function(Hm, beta = TRUE, V=FALSE, gamma = FALSE,
                               omega=FALSE, title = "Model Convergence"){

  mpost <- convertToCodaObject(Hm)

  d <-
    bind_rows(
      effectiveSize(mpost$Beta) %>%
        as_tibble() %>% mutate(fit_statistic = "ess", variable = "beta"),
      gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
        as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
        mutate(variable = "beta", fit_statistic = "psrf")
    )

  if(V) {
    d <- d %>%
      bind_rows(
        effectiveSize(mpost$V) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "V")) %>%
      bind_rows(gelman.diag(mpost$V, multivariate=FALSE)$psrf%>%
                  as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
                  mutate(variable = "V", fit_statistic = "psrf")
      )
  }

  if(gamma) {
    d <- d %>%
      bind_rows(
        effectiveSize(mpost$Gamma) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "gamma")) %>%
      bind_rows(gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf%>%
                  as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
                  mutate(variable = "gamma", fit_statistic = "psrf")
      )
  }

  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }

    d <- d %>%
      bind_rows(
        effectiveSize(tmp) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "omega")) %>%
      bind_rows(gelman.diag(tmp, multivariate=FALSE)$psrf%>%
                  as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
                  mutate(variable = "omega", fit_statistic = "psrf")
      )
  }

  vline_df <- data.frame(fit_statistic = c("ess", "psrf"),
                         xintercept = c(length(mpost$Beta)*nrow(mpost$Beta[[1]]),
                                        1.01))


  ggplot(d, aes(x=value)) +
    geom_histogram(bins=70) +
    geom_vline(data = vline_df, aes(xintercept = xintercept), color="red", lty=2)+
    facet_grid(variable~fit_statistic, scales='free') +
    ggtitle(title)

}

