# Helper functions for analysis.Rmd

prepare <- function(colon, endpoint) {
  # endpoint: 1 = recurrence, 2 = death. (Index explicitly: inside subset(), `etype` would refer to the column.)
  d <- colon[colon$etype == endpoint, ]
  stopifnot(nrow(d) == length(unique(d$id)))
  d$years <- d$time / 365.25
  d$rx <- factor(d$rx, levels = c("Obs", "Lev", "Lev+5FU"))
  d$node4 <- factor(ifelse(d$node4 == 1, ">4 nodes", "1-4 nodes"), levels = c("1-4 nodes", ">4 nodes"))
  d$extent_grp <- factor(ifelse(d$extent >= 3, "Serosa or beyond", "Submucosa/muscle"),
                         levels = c("Submucosa/muscle", "Serosa or beyond"))
  d$differ_grp <- factor(ifelse(is.na(d$differ), "Unknown", ifelse(d$differ == 3, "Poor", "Well/moderate")),
                         levels = c("Well/moderate", "Poor", "Unknown"))
  d$sex <- factor(ifelse(d$sex == 1, "Male", "Female"), levels = c("Female", "Male"))
  d$age10 <- d$age / 10
  d
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

baseline_table <- function(d) {
  arms <- levels(d$rx)
  num_row <- function(label, x) {
    vals <- sapply(arms, function(a) sprintf("%.1f (%.1f)", mean(x[d$rx == a], na.rm = TRUE), sd(x[d$rx == a], na.rm = TRUE)))
    p <- summary(aov(x ~ d$rx))[[1]][["Pr(>F)"]][1]
    data.frame(characteristic = label, t(vals), p_value = fmt_p(p), check.names = FALSE)
  }
  pct_row <- function(label, flag) {
    vals <- sapply(arms, function(a) sprintf("%d (%.0f%%)", sum(flag[d$rx == a], na.rm = TRUE),
                                             100 * mean(flag[d$rx == a], na.rm = TRUE)))
    p <- chisq.test(table(d$rx, flag))$p.value
    data.frame(characteristic = label, t(vals), p_value = fmt_p(p), check.names = FALSE)
  }
  n_row <- data.frame(characteristic = "Patients, n", t(sapply(arms, function(a) as.character(sum(d$rx == a)))),
                      p_value = "", check.names = FALSE)
  rbind(n_row,
        num_row("Age, years (mean, SD)", d$age),
        pct_row("Male", d$sex == "Male"),
        num_row("Positive lymph nodes (mean, SD)", d$nodes),
        pct_row(">4 positive nodes", d$node4 == ">4 nodes"),
        pct_row("Bowel obstruction", d$obstruct == 1),
        pct_row("Perforation", d$perfor == 1),
        pct_row("Tumour adherent to nearby organs", d$adhere == 1),
        pct_row("Poorly differentiated", d$differ_grp == "Poor"),
        pct_row("Extends to serosa or beyond", d$extent_grp == "Serosa or beyond"),
        pct_row("Surgery >20 days before registration", d$surg == 1))
}

ARM_COLORS <- c("Obs" = "#b8b7b2", "Lev" = "#eb6834", "Lev+5FU" = "#2a78d6")

km_plot <- function(fit, title, ylab) {
  td <- broom::tidy(fit)
  td$arm <- sub("rx=", "", td$strata)
  start <- data.frame(time = 0, estimate = 1, conf.low = 1, conf.high = 1, arm = unique(td$arm))
  td <- dplyr::bind_rows(start, td)
  ggplot(td, aes(time, estimate, colour = arm, fill = arm)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.10, colour = NA,
                stat = "identity", position = "identity") +
    scale_colour_manual(values = ARM_COLORS) + scale_fill_manual(values = ARM_COLORS) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = title, x = "Years since randomization", y = ylab, colour = "Arm", fill = "Arm") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top", panel.grid.minor = element_blank())
}

km_summary <- function(fit, times) {
  s <- summary(fit, times = times)
  out <- data.frame(arm = sub("rx=", "", s$strata), year = s$time,
                    survival = sprintf("%.1f%%", 100 * s$surv),
                    ci_95 = sprintf("%.1f%% - %.1f%%", 100 * s$lower, 100 * s$upper),
                    at_risk = s$n.risk)
  med <- summary(fit)$table[, "median"]
  attr(out, "median") <- med
  knitr::kable(out, row.names = FALSE)
}

hr_table <- function(model) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    dplyr::transmute(term, HR = round(estimate, 2),
                     CI_95 = sprintf("%.2f - %.2f", conf.low, conf.high),
                     p_value = fmt_p(p.value), estimate, conf.low, conf.high)
}

forest_plot <- function(hr, title) {
  labels <- c("rxLev" = "Levamisole vs observation", "rxLev+5FU" = "Lev + 5-FU vs observation",
              "node4>4 nodes" = ">4 positive nodes", "extent_grpSerosa or beyond" = "Extends to serosa or beyond",
              "obstruct" = "Bowel obstruction", "differ_grpPoor" = "Poorly differentiated",
              "differ_grpUnknown" = "Differentiation unknown", "sexMale" = "Male", "age10" = "Age (per 10 years)")
  hr$label <- factor(labels[hr$term], levels = rev(labels))
  ggplot(hr, aes(estimate, label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "#52514e") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, colour = "#2a78d6") +
    geom_point(size = 3, colour = "#2a78d6") +
    scale_x_log10() +
    labs(title = title, x = "Hazard ratio (log scale; <1 = lower risk)", y = NULL) +
    theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
}

pct  <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
pfmt <- function(chisq, df) format.pval(pchisq(chisq, df, lower.tail = FALSE), digits = 2, eps = 1e-6)
zph_df <- function(z) {
  data.frame(term = rownames(z$table), chisq = round(z$table[, "chisq"], 2), df = z$table[, "df"],
             p = signif(z$table[, "p"], 2), row.names = NULL)
}
