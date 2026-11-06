# Compare transient sensitivities, elasticities and random-design tLTREs
# across IPM posterior outputs, using side-by-side posterior violins.
#
# Packages: tidyverse, patchwork, paletteer
# Outputs: 6 PDFs (summed and age-specific for each of 3 metrics)
#
# IMPORTANT: Sensitivity derivatives below assume the demographic equation
# lambda = pYF*sYF + pSA*sSA + sum_a pAD[a] *
#          (sAD[a] + 0.5*BR[a]*sPY[a]*(sAD[a] if maternalSurvival else 1)).
# Verify this equation against EACH fitted IPM before interpreting comparisons.

compareModels_violin <- function(
    postPaths, modelNames, plotFolder = "figures/paramPop",
    nYear = 18, minYear = 2008, maxYear = minYear + nYear - 1,
    nAge = 19, t.period = NULL,
    splitAdults = TRUE, oneProp = TRUE,
    maternalSurvival = rep(TRUE, length(postPaths)),
    plotDraws = 3000, seed = 1,
    returnSumData = TRUE, returnPlotData = FALSE) {
  
  
  # for testing purposes
  nYear = 18
  minYear = 2008
  maxYear = minYear + nYear - 1
  nAge = 19
  t.period = NULL
  seed = 1
  postPaths = c(
    "results/IPM_CJSen_RSen_AB_DynDens_dCJS_Dave12_S33_R22_B22_areaHR.rds",
    "results/IPM_CJSen_RSen_AB_DynDens_dCJS_Dave12_S33_R22_B22_femaleNt.rds",
    "results/IPM_CJSen_RSen_AB_DynDens_dCJS_Dave12_S33_R22_B22_survMoms.rds",
    "results/IPM_CJSen_RSen_AB_DynDens_dCJS_Dave12_S33_R22_B22_femaleNt_survMoms.rds"
  )
  modelNames = c("IPM_areaHR", "IPM_femaleNt", "IPM_survMoms", "IPM_both")
  plotFolder = "figures/paramPop"
  splitAdults = TRUE
  oneProp = TRUE
  # Set these per model ONLY if the demographic equations warrant it:
  maternalSurvival = c(FALSE, FALSE, TRUE, TRUE)
  plotDraws = 3000
  returnSumData = TRUE
  returnPlotData = FALSE
  
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(paletteer)
  })

  # ---- Checks and setup --------------------------------------------------
  nModels <- length(postPaths)
  stopifnot(nModels > 0, length(modelNames) == nModels,
            length(maternalSurvival) == nModels,
            !anyNA(maternalSurvival), is.logical(maternalSurvival),
            length(unique(modelNames)) == nModels,
            nYear >= 3, nAge >= 2,
            maxYear - minYear + 1 == nYear)
  if (splitAdults && nAge < 10) {
    stop("splitAdults = TRUE requires nAge >= 10.")
  }
  missingFiles <- postPaths[!file.exists(postPaths)]
  if (length(missingFiles)) {
    stop("Posterior file(s) not found:\n", paste(missingFiles, collapse = "\n"))
  }
  if (is.null(t.period)) {
    transitions <- seq_len(nYear - 1)
  } else {
    transitions <- sort(unique(as.integer(t.period)))
    if (!length(transitions) || anyNA(transitions) ||
        any(transitions < 1 | transitions > (nYear - 1))) {
      stop("t.period must contain transition indices between 1 and nYear - 1.")
    }
  }
  if (length(transitions) < 2) {
    stop("At least two transitions are needed to estimate tLTRE covariance.")
  }
  if (!is.numeric(plotDraws) || length(plotDraws) != 1 ||
      is.na(plotDraws) || plotDraws < 1) {
    stop("plotDraws must be a positive number (or Inf).")
  }
  dir.create(plotFolder, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  modelCols <- stats::setNames(
    as.character(paletteer::paletteer_c("grDevices::Temps", nModels)),
    modelNames
  )

  # ---- Helpers -----------------------------------------------------------
  posterior_matrix <- function(path) {
    post <- readRDS(path)
    if (is.matrix(post)) {
      out <- post
    } else if (inherits(post, "mcmc.list") || is.list(post)) {
      out <- do.call(rbind, lapply(post, as.matrix))
    } else {
      out <- as.matrix(post)
    }
    if (!is.numeric(out) || is.null(colnames(out))) {
      stop("Posterior RDS must provide a numeric matrix with parameter names: ", path)
    }
    # Accept both 'BR[1, 2]' and 'BR[1,2]' styles.
    colnames(out) <- gsub("[[:space:]]", "", colnames(out))
    if (anyDuplicated(colnames(out))) stop("Duplicate posterior column names: ", path)
    out
  }

  get_columns <- function(out, keys) {
    ix <- match(keys, colnames(out))
    if (anyNA(ix)) {
      stop("Missing posterior parameter(s): ",
           paste(head(keys[is.na(ix)], 8), collapse = ", "))
    }
    out[, ix, drop = FALSE]
  }

  get_time <- function(out, parameter, nT) {
    get_columns(out, paste0(parameter, "[", seq_len(nT), "]"))
  }

  get_age_time <- function(out, parameter, nT) {
    keys <- unlist(lapply(seq_len(nT), function(t) {
      paste0(parameter, "[", seq_len(nAge), ",", t, "]")
    }), use.names = FALSE)
    array(get_columns(out, keys), dim = c(nrow(out), nAge, nT))
  }

  mean_age <- function(x, times) {
    apply(x[, , times, drop = FALSE], c(1, 2), mean)
  }
  mean_time <- function(x, times) {
    rowMeans(x[, times, drop = FALSE])
  }

  # Flatten one posterior vector per parameter/age and then make sums DRAW
  # BY DRAW (never sum medians or credible-interval endpoints).
  flatten_metric <- function(metric) {
    flat <- list()
    for (nm in names(metric)) {
      x <- metric[[nm]]
      if (is.matrix(x)) {
        for (a in seq_len(ncol(x))) flat[[paste0(nm, "_", a)]] <- x[, a]
      } else {
        flat[[nm]] <- as.numeric(x)
      }
    }
    for (nm in c("BR", "sPY", "sAD", "pAD")) {
      sum_keys <- function(ageIndex) {
        keys <- paste0(nm, "_", ageIndex)
        Reduce(`+`, flat[keys])
      }
      flat[[paste0(nm, "_all")]] <- sum_keys(seq_len(nAge))
      if (nAge >= 9) flat[[paste0(nm, "_2to9")]] <- sum_keys(2:9)
      if (nAge >= 10) flat[[paste0(nm, "_10up")]] <- sum_keys(10:nAge)
    }
    # Combined age-structure contribution is the sum for each draw.
    # For sensitivities/elasticities it is a descriptive sum only.
    flat$p_all <- flat$pYF + flat$pSA + flat$pAD_all
    flat
  }

  summary_flat <- function(flat, model) {
    dplyr::bind_rows(lapply(names(flat), function(v) {
      x <- flat[[v]]
      q <- stats::quantile(x, c(.025, .5, .975), na.rm = TRUE)
      tibble(Model = model, Variable = v,
             Parameter = sub("_.*$", "", v),
             Age = if (grepl("_", v, fixed = TRUE)) sub("^[^_]*_", "", v) else NA_character_,
             lCI = unname(q[1]), median = unname(q[2]), uCI = unname(q[3]))
    }))
  }

  # Only thin the draws used to DRAW violins. Summaries use all draws.
  long_for_plot <- function(flat, model, drawIndex, keep) {
    dplyr::bind_rows(lapply(intersect(keep, names(flat)), function(v) {
      tibble(Model = model, draw = drawIndex,
             Variable = v, Value = flat[[v]][drawIndex])
    }))
  }

  # ---- Categories, in the same order/style as plotLTRE_random ------------
  if (splitAdults) {
    adultS <- c("sAD_2to9", "sAD_10up")
    adultP <- c("pAD_2to9", "pAD_10up")
    adultTypesS <- c("Survival of adults (2–9)", "Survival of adults (10+)")
    adultTypesP <- c("Proportion of adults (2–9)", "Proportion of adults (10+)")
    adultLabelsS <- c("Survival of\nadults (2–9)", "Survival of\nadults (10+)")
    adultLabelsP <- c("Prop. of\nadults (2–9)", "Prop. of\nadults (10+)")
  } else {
    adultS <- "sAD_all"
    adultP <- "pAD_all"
    adultTypesS <- "Survival of adults"
    adultTypesP <- "Proportion of adults"
    adultLabelsS <- "Survival of\nadults"
    adultLabelsP <- "Prop. of\nadults"
  }

  sumKeys <- c("BR_all", "sPY_all", "sYF", "sSA", adultS)
  sumTypes <- c("Birth rate", "Survival of pouch young",
                "Survival of young-at-foot", "Survival of subadults", adultTypesS)
  sumLabels <- c("Birth\nrate", "Survival of\npouch young",
                 "Survival of\nyoung-at-foot", "Survival of\nsubadults", adultLabelsS)
  if (oneProp) {
    sumKeys <- c(sumKeys, "p_all")
    sumTypes <- c(sumTypes, "Age structure")
    sumLabels <- c(sumLabels, "Age\nstructure")
  } else {
    sumKeys <- c(sumKeys, "pYF", "pSA", adultP)
    sumTypes <- c(sumTypes, "Proportion of young-at-foot",
                  "Proportion of subadults", adultTypesP)
    sumLabels <- c(sumLabels, "Prop. of\nyoung-at-foot",
                   "Prop. of\nsubadults", adultLabelsP)
  }

  ageKeys <- list(
    survival = c("sYF", "sSA", paste0("sAD_", 2:nAge)),
    birth = paste0("BR_", 2:nAge),
    pouch = paste0("sPY_", 2:nAge),
    structure = c("pYF", "pSA", paste0("pAD_", 2:nAge))
  )
  ageLabels <- list(
    survival = parse(text = c("italic(S)[0]", "italic(S)[1]",
                              paste0("italic(S)[", 2:nAge, "]"))),
    birth = parse(text = paste0("italic(B)[", 2:nAge, "]")),
    pouch = parse(text = paste0("italic(SP)[", 2:nAge, "]")),
    structure = parse(text = c("italic(P)[0]", "italic(P)[1]",
                               paste0("italic(P)[", 2:nAge, "]")))
  )
  keepVars <- unique(c(sumKeys, unlist(ageKeys, use.names = FALSE)))

  # ---- Process each model ------------------------------------------------
  allSum <- list(sensitivities = list(), elasticities = list(), ltre = list())
  allPlot <- list(sensitivities = list(), elasticities = list(), ltre = list())

  for (i in seq_len(nModels)) {
    message("Processing model ", i, "/", nModels, ": ", modelNames[i])
    out <- posterior_matrix(postPaths[i])
    nSamples <- nrow(out)
    if (nSamples < 2) stop("At least two posterior draws needed.")

    # Posterior arrays, dimensions: [draw, age, year/transition].
    BR  <- get_age_time(out, "BR", nYear - 1)
    sPY <- get_age_time(out, "sPY", nYear - 1)
    sAD <- get_age_time(out, "sAD", nYear - 1)
    sYF <- get_time(out, "sYF", nYear - 1)
    sSA <- get_time(out, "sSA", nYear - 1)
    nAD <- get_age_time(out, "nAD", nYear)
    nYF <- get_time(out, "nYF", nYear)
    nSA <- get_time(out, "nSA", nYear)
    nTOT <- get_time(out, "nTOT", nYear)

    pYF <- nYF / nTOT
    pSA <- nSA / nTOT
    pAD <- sweep(nAD, c(1, 3), nTOT, "/")
    lambda <- nTOT[, 2:nYear, drop = FALSE] /
      nTOT[, 1:(nYear - 1), drop = FALSE]

    # Match original defaults: rate and lambda averages over nYear-1
    # transitions; population proportions over all nYear census years.
    # For explicit t.period, align p at t+1 with vital rates at t.
    propTimes <- if (is.null(t.period)) seq_len(nYear) else transitions + 1L
    B <- mean_age(BR, transitions)
    R <- mean_age(sPY, transitions)
    A <- mean_age(sAD, transitions)
    Y <- mean_time(sYF, transitions)
    SA <- mean_time(sSA, transitions)
    PY <- mean_time(pYF, propTimes)
    PSA <- mean_time(pSA, propTimes)
    PA <- mean_age(pAD, propTimes)
    L <- mean_time(lambda, transitions)

    # Transient sensitivities, evaluated per posterior draw.
    # Formula switches for models where recruitment requires maternal survival.
    recruitment <- 0.5 * B * R
    S <- list(
      BR  = PA * (0.5 * R) * (if (maternalSurvival[i]) A else 1),
      sPY = PA * (0.5 * B) * (if (maternalSurvival[i]) A else 1),
      sYF = PY,
      sSA = PSA,
      sAD = PA * (1 + if (maternalSurvival[i]) recruitment else 0),
      pYF = Y,
      pSA = SA,
      pAD = A + recruitment * (if (maternalSurvival[i]) A else 1)
    )
    # Age class 1 is not an adult class in the original model.
    for (nm in c("BR", "sPY", "sAD", "pAD")) S[[nm]][, 1] <- 0

    # Transient elasticities, divided by observed lambda per draw.
    E <- list(
      BR  = S$BR  * B / L,
      sPY = S$sPY * R / L,
      sYF = S$sYF * Y / L,
      sSA = S$sSA * SA / L,
      sAD = S$sAD * A / L,
      pYF = S$pYF * PY / L,
      pSA = S$pSA * PSA / L,
      pAD = S$pAD * PA / L
    )

    # Random-design tLTRE for each draw:
    # C_i = s_i * sum_j[ cov_t(theta_i, theta_j) * s_j ].
    # Survival and reproduction use t, proportions use t+1, as in the
    # original function. Dimensions are time x demographic parameter.
    contNames <- c(paste0("BR_", seq_len(nAge)),
                   paste0("sPY_", seq_len(nAge)), "sYF", "sSA",
                   paste0("sAD_", seq_len(nAge)), "pYF", "pSA",
                   paste0("pAD_", seq_len(nAge)))
    C <- matrix(NA_real_, nrow = nSamples, ncol = length(contNames),
                dimnames = list(NULL, contNames))
    age_time <- function(x, s, tt) {
      # x[s, age, time] -> matrix [time, age], even for one selected time.
      t(matrix(x[s, , tt, drop = FALSE], nrow = nAge, ncol = length(tt)))
    }

    for (s in seq_len(nSamples)) {
      X <- cbind(
        age_time(BR, s, transitions),
        age_time(sPY, s, transitions),
        sYF[s, transitions], sSA[s, transitions],
        age_time(sAD, s, transitions),
        pYF[s, transitions + 1L], pSA[s, transitions + 1L],
        age_time(pAD, s, transitions + 1L)
      )
      stopifnot(ncol(X) == length(contNames))
      V <- stats::cov(X)
      d <- c(S$BR[s, ], S$sPY[s, ], S$sYF[s], S$sSA[s],
             S$sAD[s, ], S$pYF[s], S$pSA[s], S$pAD[s, ])
      C[s, ] <- d * as.vector(V %*% d)
    }
    LTRE <- as.list(as.data.frame(C, check.names = FALSE))

    metrics <- list(sensitivities = S, elasticities = E, ltre = LTRE)
    drawIndex <- if (is.infinite(plotDraws) || plotDraws >= nSamples) {
      seq_len(nSamples)
    } else {
      sort(sample.int(nSamples, size = as.integer(plotDraws)))
    }
    for (metricName in names(metrics)) {
      flat <- flatten_metric(metrics[[metricName]])
      allSum[[metricName]][[i]] <- summary_flat(flat, modelNames[i])
      allPlot[[metricName]][[i]] <- long_for_plot(flat, modelNames[i],
                                                  drawIndex, keepVars)
    }
    rm(out, BR, sPY, sAD, sYF, sSA, nAD, nYF, nSA, nTOT,
       pYF, pSA, pAD, S, E, C, LTRE, metrics)
  }

  resultsSum <- lapply(allSum, dplyr::bind_rows)
  resultsPlot <- lapply(allPlot, function(z) {
    dplyr::bind_rows(z) %>%
      dplyr::mutate(Model = factor(Model, levels = modelNames))
  })

  # ---- Plotting -----------------------------------------------------------
  violin_plot <- function(dat, keys, xLabels, yLabel, title = NULL,
                          showLegend = TRUE) {
    pp <- dat %>%
      dplyr::filter(Variable %in% keys, is.finite(Value)) %>%
      dplyr::mutate(Variable = factor(Variable, levels = keys))
    ggplot(pp, aes(x = Variable, y = Value, fill = Model, colour = Model)) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
      geom_violin(alpha = 0.5, scale = "width", trim = TRUE,
                  width = 0.85, linewidth = 0.25,
                  position = position_dodge(width = 0.88)) +
      stat_summary(fun = median, geom = "point", colour = "grey15",
                   size = 0.65, position = position_dodge(width = 0.88),
                   show.legend = FALSE) +
      scale_fill_manual(values = modelCols, drop = FALSE) +
      scale_colour_manual(values = modelCols, drop = FALSE) +
      scale_x_discrete(drop = FALSE, labels = xLabels) +
      scale_y_continuous(expand = expansion(mult = c(0.04, 0.05))) +
      labs(x = NULL, y = yLabel, title = title, fill = "Model", colour = "Model") +
      theme_bw(base_size = 11) +
      theme(
        legend.position = if (showLegend) "top" else "none",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 11, face = "plain"),
        plot.margin = margin(3, 5, 3, 5)
      )
  }

  plotFiles <- list()
  metricLabels <- c(sensitivities = "Transient sensitivity",
                    elasticities = "Transient elasticity",
                    ltre = "tLTRE contribution")
  filePrefixes <- c(sensitivities = "Compare_Sensitivities",
                    elasticities = "Compare_Elasticities",
                    ltre = "Compare_LTREContributions")

  for (metricName in names(resultsPlot)) {
    dat <- resultsPlot[[metricName]]
    label <- unname(metricLabels[metricName])
    prefix <- unname(filePrefixes[metricName])

    pSum <- violin_plot(dat, sumKeys, sumLabels, label)
    sumPath <- file.path(plotFolder, paste0(prefix, "_sum.pdf"))
    ggsave(sumPath, pSum, width = if (oneProp) 10 else 13,
           height = 5.2, units = "in", device = "pdf")

    panels <- list(
      violin_plot(dat, ageKeys$survival, ageLabels$survival,
                  label, "a) Survival", showLegend = FALSE),
      violin_plot(dat, ageKeys$birth, ageLabels$birth,
                  label, "b) Birth rate", showLegend = FALSE),
      violin_plot(dat, ageKeys$pouch, ageLabels$pouch,
                  label, "c) Survival of pouch young", showLegend = FALSE),
      violin_plot(dat, ageKeys$structure, ageLabels$structure,
                  label, "d) Population proportions", showLegend = TRUE)
    )
    pAge <- patchwork::wrap_plots(panels, ncol = 1)
    agePath <- file.path(plotFolder, paste0(prefix, "_age.pdf"))
    ggsave(agePath, pAge, width = 15, height = 15,
           units = "in", device = "pdf")

    plotFiles[[metricName]] <- c(summary = sumPath, age_specific = agePath)
  }

  if (returnSumData) {
    result <- c(resultsSum, list(plots = plotFiles))
    if (returnPlotData) result$plotData <- resultsPlot
    return(result)
  }
  invisible(plotFiles)
}



