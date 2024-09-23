#' Plot model convergence indexes
#'
#' @param beta if TRUE, plots the beta (env. filters) parameters
#' @param V if TRUE, plots the V parameters
#' @param gamma if TRUE, plots the gamma (traits) paramters
#' @param omega if TRUE, plots the omega (spp associations) parameters
#' @param title character string to customize
#' @export
gghmsc_convergence <- function(Hm,
                               beta = TRUE,
                               V=FALSE,
                               gamma = FALSE,
                               omega=FALSE,
                               title = "Model Convergence"){
  requireNamespace("Hmsc")
  requireNamespace("coda")
  requireNamespace("dplyr")
  requireNamespace("tibble")

  mpost <- Hmsc::convertToCodaObject(Hm)

  d <-
    bind_rows(
      coda::effectiveSize(mpost$Beta) |>
        tibble::as_tibble() |>
        dplyr::mutate(fit_statistic = "ess", variable = "beta"),
      coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf|>
        as_tibble() |> dplyr::rename(value = `Point est.`) |>
        mutate(variable = "beta", fit_statistic = "psrf")
    )

  if(V) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$V) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "V")) |>
      dplyr::bind_rows(coda::gelman.diag(mpost$V, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                  dplyr::mutate(variable = "V", fit_statistic = "psrf")
      )
  }

  if(gamma) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$Gamma) |>
          tibble::as_tibble() |> dplyr::mutate(fit_statistic = "ess", variable = "gamma")) |>
      dplyr::bind_rows(coda::gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                 dplyr::mutate(variable = "gamma", fit_statistic = "psrf")
      )
  }

  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }

    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(tmp) |>
          tibble::as_tibble() |> dplyr::mutate(fit_statistic = "ess", variable = "omega")) |>
      dplyr::bind_rows(coda::gelman.diag(tmp, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                  dplyr::mutate(variable = "omega", fit_statistic = "psrf")
      )
  }

  vline_df <- data.frame(fit_statistic = c("ess", "psrf"),
                         xintercept = c(length(mpost$Beta)*nrow(mpost$Beta[[1]]),
                                        1.01))


  ggplot2::ggplot(d, ggplot2::aes(x=value)) +
    geom_histogram(bins=70) +
    geom_vline(data = vline_df, ggplot2::aes(xintercept = xintercept), color="red", lty=2)+
    facet_grid(variable~fit_statistic, scales='free') +
    ggtitle(title)

}
#' Trace plots
#'
#' @param Hm Hmsc object
#' @param which Can be "beta", "gamma", or "v"
#' @export
gghmsc_traceplot <- function(Hm,
                        which = "beta"){
  requireNamespace("Hmsc")
  requireNamespace("ggmcmc")
  co <- Hmsc::convertToCodaObject(Hm)
  if(which == "beta") return(plot(co$Beta))
  if(which == "gamma") return(plot(co$Gamma))
  if(which == "v") return(plot(co$V))
}

#' Plot effective sample size
#'
#' @export
gghmsc_ess <- function(Hm,
                       beta = TRUE,
                       V=FALSE,
                       gamma = FALSE,
                       omega=FALSE,
                       title = "Model Convergence"){

  mpost <- Hmsc::convertToCodaObject(Hm)

  d <- coda::effectiveSize(mpost$Beta) |>
    tibble::as_tibble() |>
    dplyr::mutate(fit_statistic = "ess", variable = "beta")

  if(V) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$V) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "V"))

  }

  if(gamma) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$Gamma) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "gamma"))

  }

  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }

    d <- d |>
      bind_rows(
        coda::effectiveSize(tmp) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "omega"))
  }

  vline_df <- data.frame(fit_statistic = "ess",
                         xintercept = length(mpost$Beta)*nrow(mpost$Beta[[1]]))


  ggplot2::ggplot(d,ggplot2::aes(x=value)) +
    ggplot2::geom_histogram(bins=70) +
    ggplot2::geom_vline(data = vline_df,
                        ggplot2::aes(xintercept = xintercept), color="red", lty=2)+
    ggplot2::facet_grid(variable~fit_statistic, scales='free') +
    ggplot2::ggtitle(title)

}

#' Plot variance partitioning
#'
#'
#' @export
gghmsc_vp <- function(Hm,
                      title = "Variance Explained",
                      cols = NULL,
                      lut_varnames = NULL,
                      lut_sppnames = NULL){
  requireNamespace("ggtext")
  requireNamespace("Hmsc")
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  VP <- Hmsc::computeVariancePartitioning(Hm)
  mpost <- Hmsc::convertToCodaObject(Hm)

  prevalence <- colSums(Hm$Y) |>
    tibble::as_tibble(rownames = "Species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(desc(prevalence))

  mf_df <- data.frame(Species = colnames(Hm$Y)) |>
    dplyr::left_join(prevalence)

  vp_df <- VP$vals|>
    tibble::as_tibble(rownames = "variable") |>
    tidyr::pivot_longer(cols = names(as.data.frame(VP$vals)),
                 names_to = "Species",
                 values_to = "value") |>
    dplyr::left_join(prevalence) |>
    na.omit()

  if(is.vector(lut_varnames)) vp_df <- vp_df |> dplyr::mutate(variable = lut_varnames[variable])
  if(is.vector(lut_sppnames)) vp_df <- vp_df |> dplyr::mutate(Speices = lut_varnames[Species])

  vp_summary <- vp_df |>
    dplyr::group_by(variable) |>
    dplyr::summarise(value_pct = mean(value) * 100) |>
    dplyr::ungroup() |>
    dplyr::mutate(variable_pct = paste0(variable, " (", round(value_pct,1), "%)")) |>
    dplyr::select(-value_pct)
  vspp <- vp_df |>
    dplyr::filter(variable == dplyr::first(vp_df$variable |> unique())) |>
    dplyr::arrange(prevalence) |>
    dplyr::pull(Species)
  vp_order <- vp_df |>
    dplyr::filter(variable == dplyr::first(vp_df$variable |> unique())) |>
    dplyr::arrange(prevalence) |>
    dplyr::mutate(Species_f = factor(Species, levels = vspp)) |>
    dplyr::select(Species, Species_f)

  p_vp <- dplyr::left_join(vp_df, vp_order) |>
    dplyr::left_join(vp_summary) |>
    dplyr::mutate(variable = factor(variable),
           value = value) |>
    ggplot2::ggplot(ggplot2::aes(x=value,y=Species_f, fill = variable_pct)) +
    ggplot2::geom_bar(stat="identity", color = "black")+
    ggplot2::theme_classic() +
    ggplot2::scale_fill_discrete(name = "Variable\n (Avg Variance Explained)") +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of Variance Explained") +
    ggplot2::theme(legend.position = "right",
          legend.text = ggtext::element_markdown(),
          legend.title = ggplot2::element_blank(),
          legend.justification = c(1,0),
          legend.background = ggplot2::element_rect(color="black"))

  if(is.vector(cols)) p_vp <- p_vp + ggplot2::scale_fill_manual(values = cols)
  if(is.vector(title)) p_vp <- p_vp + ggplot2::ggtitle(title)

  return(p_vp)

}

#' Plot beta posterior estimates as colored boxes
#'
#' @param order_x order the variables on the x axis
#' @param grouping_var_y group the species on the Y axis by traits
#' @examples
#' data("Hm")
#' library(ggplot2)
#' label_df <- data.frame(x = rep(12.75), y = c(1.5,8.5, 10.5, 16.5, 19.5),
#'                        label = c("IF", "IG", "NF", "NG", "NW"))
#' gghmsc_beta(Hm, grouping_var_y = Hm$TrData$fg,
#'             spp_exclude = c("Agropyron cristatum", "Bassia prostrata")) +
#'   geom_hline(yintercept = c(7.5, 9.5, 15.5, 18.5)) +
#'   geom_text(data = label_df, aes(x=x,y=y, label=label), angle =90)
#'
#'
#' @export
gghmsc_beta <- function(Hm,
                        order_x = NA,
                        grouping_var_y = NA,
                        grouping_var_y_2 = NA,
                        n_gvy = 1,
                        spp_exclude = NA,
                        support_level = 0.89,
                        lut_varnames = NULL,
                        lut_sppnames = NULL,
                        no_intercept = TRUE,
                        title = NA){
  requireNamespace("Hmsc")
  requireNamespace("tibble")
  requireNamespace('dplyr')
  requireNamespace('ggplot2')
  requireNamespace('ggtext')
  postBeta <- Hmsc::getPostEstimate(Hm, parName = "Beta")

  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }

  means <- postBeta$mean |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = c(covNames)) |>
    tidyr::pivot_longer(cols=colnames(postBeta$mean),
                        names_to = "Species", values_to = "Mean")

  supported <- postBeta$support |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = covNames) |>
    tidyr::pivot_longer(cols=colnames(postBeta$support),
                 names_to = "Species",
                 values_to = "Support") |>
    dplyr::filter(Support > support_level | Support < (1-support_level),
           env_var != "(Intercept)") |>
    dplyr::left_join(means, by = c("env_var", "Species"))|>
    dplyr::mutate(sign = ifelse(Mean>0, "+", "-"))

  sp_sorted <- colSums(Hm$Y) |>
    tibble::as_tibble(rownames = "Species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(prevalence) |>
    dplyr::pull(Species)

  if(!is.na(grouping_var_y[1])){
    if(n_gvy == 1){
    sp_sorted <- Hm$TrData |>
      tibble::as_tibble(rownames = "species") |>
      dplyr::arrange(grouping_var_y) |>
      dplyr::pull(species)}else{
        sp_sorted <- Hm$TrData |>
          tibble::as_tibble(rownames = "species") |>
          dplyr::arrange(grouping_var_y, grouping_var_y_2) |>
          dplyr::pull(species)
      }
  }

  vp_order <-   colSums(Hm$Y) |>
    tibble::as_tibble(rownames = "Species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(prevalence) |>
    dplyr::mutate(Species_f = factor(Species, levels = sp_sorted)) |>
    dplyr::filter(Species %in% supported$Species)


  supported <- supported |>
    dplyr::left_join(vp_order)#

  if(!is.na(order_x[1])){
    var_order <- tibble::tibble(var1 = order_x) |>
      dplyr::mutate(order = 1:length(order_x)) |>
      dplyr::mutate(order_f = factor(var1, levels = order_x)) |>
      dplyr::rename(env_var = var1)

    supported <- supported |>
      left_join(var_order) |>
      dplyr::select(-env_var) |>
      dplyr::rename(env_var = order_f)

  }

  if(is.vector(lut_varnames)) supported <- supported |> dplyr::mutate(env_var = lut_varnames[env_var])
  if(is.vector(lut_sppnames)) supported <- supported |> dplyr::mutate(Species = lut_varnames[Species])
  if(no_intercept) supported <- supported |> dplyr::filter(env_var != "(Intercept)")
  if(class(spp_exclude) == "character") supported <- dplyr::filter(supported, !Species %in% spp_exclude)
  p_beta <- supported |>
    ggplot2::ggplot(ggplot2::aes(x=env_var,y=reorder(Species_f,Species))) +
    ggplot2::geom_tile(lwd=.5,ggplot2::aes(fill = Mean, color = sign)) +
    ggplot2::theme_classic()+
    ggplot2::scale_fill_steps2() +
    ggplot2::scale_color_manual(values = c(("red"), ("blue"))) +
    ggplot2::guides(color = "none")+
    ggplot2::scale_x_discrete(expand = c(0,1)) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown(angle=45, vjust=1,hjust = 1),
          legend.position = "right",
          panel.grid.major.y = ggplot2::element_line(color = "grey", linetype=3),
          plot.background = ggplot2::element_rect(color="black"),
          plot.title = ggplot2::element_text(hjust = 1, face = "bold")) +
    ggplot2::xlab("Environmental Filters")+
    ggplot2::ylab("Species")

  if(!is.na(title)) p_beta <- p_beta + ggplot2::ggtitle(title)

  return(p_beta)
}

#' Plot beta estimates using PDFs
#'
#' @export
gghmsc_beta2 <- function(Hm,
                         lut_gensp=NA,
                         included_variables = NA,
                         lut_ivars = NA){
  require(dplyr)
  require(ggthemes)
  require(tidyr)
  require(ggplot2)
  require(ggnewscale)
  require(ggmcmc)
  c<-convertToCodaObject(Hm)
  mbc <- ggmcmc::ggs(c$Beta) |>
    separate(.,
             col = "Parameter",
             into = c("var", "x1", "gen", "sp", "x2"),
             sep = " ") |>
    dplyr::select(-x1, -x2,-sp) |>
    mutate(gensp = paste(gen),
           var = str_remove_all(var, "B\\["),
           gensp = str_remove_all(gensp, " \\(S\\d{2}\\)\\]"))|>
    filter(var != "(Intercept)") |>
    group_by(var, gensp, Chain) |>
    mutate(value = scale(value,center = F),
           sign = ifelse(value>0, "positive", "negative"),
           median_value = median(value)) |>
    filter(value<4 & value>-4) |>
    ungroup() |>
    left_join(Hm$TrData |> tibble::rownames_to_column("gensp"))

  if(!is.na(lut_gensp)){mbc <- mbc |>
    mutate(gensp = lut_gensp[gensp])}

  if(any(!is.na(included_variables))){
    mbc <- filter(mbc, var %in% included_variables)
  }

  prevalence <- Hm$Y |>
    as_tibble(rownames = "plot") |>
    pivot_longer(cols = names(.)[2:ncol(.)], names_to = "gen") |>
    group_by(gen) |>
    summarise(prevalence = sum(value),
              prev_pct = sum(value)/n()*100) |>
    mutate(prev_pct = ifelse(prev_pct<1, round(prev_pct,1), round(prev_pct))) |>
    ungroup()

  mbc <- mbc |>
    left_join(prevalence) |>
    mutate(gensp = paste0(gensp, " (", prev_pct, ")")) |>
    mutate(gensp = str_replace_all(gensp,"0.5", ".5"))

  vp_order <- mbc |>
    left_join(prevalence)  |>
    filter(var == first(mbc$var |> unique()),
           Iteration ==1, Chain==1)|>
    arrange(prevalence) |>
    mutate(gensp_f = factor(gensp, levels = .$gensp)) |>
    dplyr::select(gensp, gensp_f, gen)

  if(any(!is.na(included_variables))){
    mbc <- mutate(mbc, var = lut_ivars[var])
  }
  p <- ggplot(mbc |> left_join(vp_order),
             ggplot2::aes(x=value, y = gensp,#_f,
                  group=as.factor(Chain))) +
    scale_color_manual(values = (c("white", "grey90")))+
    ggdist::stat_slab(height=2,  lwd = .5, #alpha = 0.95,
                      color = "black",
                     ggplot2::aes(fill = after_stat(x>0),
                          alpha = exp(abs(median_value))))+
    facet_wrap(~var, scales = "free_x", nrow=1, ncol=length(unique(mbc$var))) +
    scale_alpha_continuous(range = c(0,1.25*(1/length(unique(mbc$Chain)))))+
    theme_classic() +
    guides(fill="none", alpha="none", color = "none")+
    geom_vline(xintercept=0, col="black", lty=2) +
    ggnewscale::new_scale_fill() +
    xlab("Scaled Effect on Occurrence Probability") +
    ylab("Species or Species Group (% Prevalence)") +
    theme(panel.spacing.x = unit(-1, "lines"),
          # panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12))#;p
  return(p)
}

#' plot trait covariate relationships
#'
#' @export
gghmsc_gamma <- function(Hm,
                         support_level = 0.89,
                         no_intercept = TRUE,
                         title = "Effects on Traits"){
  requireNamespace('magrittr')
  requireNamespace('stringr')
  requireNamespace('Hmsc')
  requireNamespace('dplyr')
  requireNamespace('tibble')
  requireNamespace('ggplot2')
  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }

  trNames = character(Hm$nt)
  trNamesNumbers = c(T,F)
  for (i in 1:Hm$nt) {
    sep = ""
    if (trNamesNumbers[1]) {
      trNames[i] = paste(trNames[i], Hm$trNames[i], sep = sep)
      sep = " "
    }
    if (trNamesNumbers[2]) {
      trNames[i] = paste(trNames[i], sprintf("(T%d)", i),
                         sep = sep)
    }
  }

  trNames <- stringr::str_remove_all(trNames, "yes") |>
    stringr::str_remove("pp")

  postGamma = Hmsc::getPostEstimate(Hm, parName="Gamma")

  means_gamma <- postGamma$mean |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = c(covNames))

  means_gammal <- means_gamma |>
    tidyr::pivot_longer(cols=names(means_gamma)[2:ncol(means_gamma)],
                        names_to = "Trait", values_to = "Mean")

  lut_gamma <- trNames
  names(lut_gamma) <- unique(means_gammal$Trait)

  supported_gamma <- postGamma$support |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = covNames) %>%
    tidyr::pivot_longer(cols=names(.)[2:ncol(.)],
                 names_to = "Trait",
                 values_to = "Support") |>
    dplyr::filter(Support > support_level |Support< (1-support_level),
           env_var != "(Iintercept)") |>
    dplyr::left_join(means_gammal, by = c("env_var", "Trait"))|>
    dplyr::mutate(sign = ifelse(Mean>0, "+", "-"),
           Trait = lut_gamma[Trait])|>
    dplyr::filter(Trait != "(Intercept)",
                  env_var != "(Intercept)")

  p_gamma <- supported_gamma |>
    ggplot2::ggplot(ggplot2::aes(x=env_var,y=(Trait), fill = Mean, color = sign)) +
    ggplot2::geom_tile(lwd=.5) +
    ggpubr::theme_pubclean()+
    ggplot2::scale_fill_steps2() +
    ggplot2::scale_color_manual(values = c(("red"), ("blue"))) +
    ggplot2::guides(color = "none")+
    ggplot2::theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
          # axis.title = element_blank(),
          legend.position = "right",
          plot.background = element_rect(color="black"),
          plot.title = element_text(hjust = 1, face = "bold")) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("Environmental Filters") +
    ggplot2::ylab("Traits")

  return(p_gamma)
}

#' plot trait covariate associations
#'
#' @export
gghmsc_gamma2 <- function(Hm,
                          lut_varnames = NA){
  c<-convertToCodaObject(Hm)
  mbc <- ggmcmc::ggs(c$Gamma) |>
    separate(.,
             col = "Parameter",
             into = c("var", "x1", "trait", "x2"),
             sep = " ") |>
    dplyr::select(-x1, -x2) |>
    mutate(var = str_remove_all(var, "G\\["),
           trait = str_remove_all(trait, "yes"),
           trait = str_remove_all(trait, "cots"),
           trait = str_replace_all(trait, "originN", "native"),
           trait = str_to_title(trait))|>
    filter(var != "(Intercept)", trait != "(Intercept)")

  if(!is.na(lut_varnames)){
    mbc <- mutate(mbc, var = lut_varnames[var])
  }

  p <- ggplot(mbc,
              ggplot2::aes(x=value, y = trait,
                  fill=as.factor(Chain))) +
    ggdist::stat_dist_interval(alpha=0.5) +
    facet_wrap(~var, scales = "free_x", nrow=2,
               ncol=ceiling(length(unique(mbc$var))/2)) +
    theme_classic() +
    guides(fill="none")+
    geom_vline(xintercept=0, col="black", lty=2) +
    xlab("Effect on Occurrence Probability") +
    ylab("Trait") +
    theme(strip.text = element_markdown())
  return(p)
}

#' create a correlation plot for species associations
#'
#' @export
gghmsc_omega <- function(Hm,
                         support_level = 0.89,
                         hc.method = "single",
                         hc.order = TRUE,
                         lut_gensp = NA,
                         axis_text_colors_x = "black",
                         axis_text_colors_y = "black",
                         title = "Residual Species Associations"){
  requireNamespace("ggcorrplot")
  OmegaCor = Hmsc::computeAssociations(Hm)

  hmdf_mean <- OmegaCor[[1]]$mean |>
    as.matrix()
  if(!is.na(lut_gensp[1])){
  rownames(hmdf_mean) <- lut_gensp[rownames(hmdf_mean)]
  colnames(hmdf_mean) <- lut_gensp[colnames(hmdf_mean)]
}

  pcor1<- ggcorrplot::ggcorrplot(hmdf_mean,
                                 type = "lower",
                                 hc.order = hc.order,
                                 hc.method = hc.method,
                                 colors = c("red", "grey90", "blue"),
                                 title = title,
                                 tl.srt = 90,
                                 legend.title = "Residual\nCorrelation") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme(plot.background = ggplot2::element_rect(color="black"),
          legend.position = c(0,1),
          axis.text.x = ggplot2::element_text(color = axis_text_colors_x,
                                     vjust = .05,
                                     face = "italic"),
          axis.text.y = ggplot2::element_text(color = axis_text_colors_y,
                                     face = "italic"),
          legend.justification = c(0,1),
          legend.background = ggplot2::element_rect(color="black"),
          plot.title = ggplot2::element_text(hjust = 1, face = "bold"))

  return(pcor1)
}

#' get Tjur R2 for each species
#'
#' @param which either all, r2 or named. all is a histogram with R2 and RSME. r2 is just R2. named is bar plot each species named on the  y axis.
#' @export
gghmsc_fit <- function(Hm, which = "r2", sp_names = "none",
                       title = "Variance Explained"){
  requireNamespace("Hmsc")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  requireNamespace("ggplot2")
  requireNamespace("ggtext")
  requireNamespace("dplyr")
  requireNamespace("magrittr")

  mpost <- Hmsc::convertToCodaObject(Hm)
  preds <- Hmsc::computePredictedValues(Hm)
  MF <- Hmsc::evaluateModelFit(hM=Hm, predY=preds)
  df <- MF |>
    tibble::as_tibble()
  df <- tidyr::pivot_longer(df, cols = names(df))
  spp <- colnames(Hm$Y)

  if(which == "named"){
    means <- df |>
      dplyr::filter(name %in% c("TjurR2", "R2")) |>
      dplyr::mutate(species = spp |> str_replace_all("_", ' ')) |>
      na.omit()
    if(sp_names[1] != "none") means <- dplyr::mutate(means, species = sp_names[species])
    return(ggplot2::ggplot(means, ggplot2::aes(x=value, y=reorder(species, value))) +
             ggplot2::geom_bar(stat = "identity") +
             ggplot2::ggtitle(title) +
             ggplot2::ylab("Species") +
             ggplot2::xlab("R<sup>2</sup>") +
             ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                            axis.text.y = ggplot2::element_text(face = "italic")))
  }

  if(which == "all"){
    means <- df |>
      dplyr::group_by(name) |>
      dplyr::summarise(mean = mean(value)) |>
      dplyr::ungroup()

    return(ggplot2::ggplot(df) +
             ggplot2::geom_histogram(aes(x=value)) +
             ggplot2::facet_wrap(~name) +
             ggplot2::geom_text(data = means, ggplot2::aes(label = paste("Avg =", round(mean, 2))), x=.75, y=4) +
             ggplot2::ggtitle(title))}

  if(which == "r2"){
    means <- df |>
      dplyr::filter(name == "R2") |>
      dplyr::group_by(name) |>
      dplyr::summarise(mean = mean(value),
                max = max(value),
                min = min(value)) |>
      dplyr::ungroup()

    return(ggplot2::ggplot(df%>% dplyr::filter(name == "R2") ) +
             ggplot2::geom_histogram(aes(x=value),bins = 15) +
             ggplot2::ggtitle(paste("Tjur R<sup>2</sup>, Avg:", round(means$mean, 2),
                                    ", Range: ", round(means$min, 2)," - ",round(means$max, 2))) +
             ggplot2::ggtitle(title) +
             ggplot2::theme(plot.title = ggtext::element_markdown()))
  }
}



lut_varnames <- c("Cheatgrass Cover" = "B_tectorum",
                     "Shrub Cover" = "shrub_pq",
                     "Native Perennial Cover" = "perennial_herbaceous",
                     "Grazing Intensity" = "grazing_intensity",
                     "Later Ignitions" = "StartMonth",
                     "Time Since Fire" = "tsf",
                     'Pre-Fire AET' = 'max_aet_z_preceeding_fire',
                     'Pre-Fire Tmin' = 'median_tmn_z_preceeding_fire',
                     'Post-Fire AET' = 'min_aet_z_after_fire',
                     'Pre-Fire AET' = 'min_aet_z_preceeding_fire',
                     'max_def_z_after_fire' = 'Post-Fire CWD',
                     "elevation_m" = 'Elevation',
                     'Warmer Aspects' = "folded_aspect")
