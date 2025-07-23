## Time lag analyze before GAM
library(mgcv)
library(ggplot2)
library(dplyr)
library(patchwork)

# import the data in NLME.xls and Z-score nomarilized
NLME_scaled <- NLME %>%
  mutate(NH4_total_supply = NH4_supply_zoo + NH4_supply_bacteria) %>%
  mutate(across(c(niche_overlap, NH4_inventory, NH4_total_supply, mu_AOA, mu_phyto, par), scale))

niche <- as.vector(NLME_scaled$niche_overlap)
NH4_inv <- as.vector(NLME_scaled$NH4_inventory)
NH4_sup <- as.vector(NLME_scaled$NH4_total_supply)
mu_aoa <- as.vector(NLME_scaled$mu_AOA)
mu_phyto <- as.vector(NLME_scaled$mu_phyto)
par <- as.vector(NLME_scaled$par)

# CCF analyze: ccf(x, y) is checking the lag relationship of variable x with respect to y
# if lag > 0, x is prior to y, which means x maybe the driven factor of y
# ccf(driver, response)
par(mfrow = c(3, 2))  # 3x2 subplot panel structure

ccf(NH4_inv, niche, lag.max = 50, main = "CCF: NH₄ Inventory → Niche Overlap")
ccf(NH4_sup, niche, lag.max = 50, main = "CCF: NH₄ Total Supply → Niche Overlap")
ccf(mu_aoa, niche, lag.max = 50, main = "CCF: μ_AOA → Niche Overlap")
ccf(mu_phyto, niche, lag.max = 50, main = "CCF: μ_phyto → Niche Overlap")
ccf(par, niche, lag.max = 50, main = "CCF: PAR → Niche Overlap")

extract_ccf_peak <- function(x, y, lag.max = 50) {
  cc <- ccf(x, y, plot = FALSE, lag.max = lag.max)
  lag_at_max <- cc$lag[which.max(abs(cc$acf))]
  max_corr <- cc$acf[which.max(abs(cc$acf))]
  return(data.frame(Lag = lag_at_max, Correlation = max_corr))
}

ccf_table <- rbind(
  extract_ccf_peak(NH4_inv, niche),
  extract_ccf_peak(NH4_sup, niche),
  extract_ccf_peak(mu_aoa, niche),
  extract_ccf_peak(mu_phyto, niche),
  extract_ccf_peak(par, niche)
)
rownames(ccf_table) <- c("NH4_inventory", "NH4_total_supply", "mu_AOA", "mu_phyto", "PAR")
ccf_table

# CCF plot (figure S2)
make_ccf_plot <- function(driver, response, name, color = "#82B29A") {
  ccf_out <- ccf(driver, response, lag.max = 50, plot = FALSE)
  df <- data.frame(
    lag = ccf_out$lag,
    correlation = ccf_out$acf
  )
  n <- length(response)

  ggplot(df, aes(x = lag, y = correlation)) +
    geom_bar(stat = "identity", fill = color, width = 1) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(2, -2) / sqrt(n), linetype = "dashed", color = "gray40") +
    labs(
      title = name,
      x = "Lag (days)",
      y = "Cross-correlation"
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold")
    )
}

p1 <- make_ccf_plot(mu_phyto, niche, "μ_Phyto")
p2 <- make_ccf_plot(mu_aoa, niche, "μ_AOA")
p3 <- make_ccf_plot(NH4_inv, niche, "NH₄⁺ Inventory")
p4 <- make_ccf_plot(NH4_sup, niche, "NH₄⁺ Total Supply")
p5 <- make_ccf_plot(par, niche, "PAR")

(p1 | p2 | p5) / (p3 | p4 | patchwork::plot_spacer())

## Time lagged GAM analyze and figure plot (figure 3)
library(ggplot2)
library(patchwork)
library(ggtext)
plot_gam_effect <- function(varname, fit_model, data, x_label_html, show_y = TRUE, color = "#F2CC8E") {
  x_seq <- seq(min(data[[varname]], na.rm = TRUE),
               max(data[[varname]], na.rm = TRUE), length.out = 200)
  
  newdata <- data[1, , drop = FALSE]
  newdata <- newdata[rep(1, 200), ]
  newdata[[varname]] <- x_seq
  
  for (v in c("mu_AOA_lag14", "mu_phyto_lag4", "NH4_inventory", 
              "NH4_total_supply_lag34", "par_lag4")) {
    if (v != varname) {
      newdata[[v]] <- mean(data[[v]], na.rm = TRUE)
    }
  }
  newdata$phase <- factor("1", levels = levels(data$phase))
  
  pred <- predict(fit_model, newdata = newdata, se.fit = TRUE)
  newdata$fit <- pred$fit
  newdata$lower <- pred$fit - 1.96 * pred$se.fit
  newdata$upper <- pred$fit + 1.96 * pred$se.fit
  
  p <- ggplot(newdata, aes_string(x = varname, y = "fit")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = color, alpha = 0.2) +
    geom_line(color = color, size = 1.2) +
    labs(
      x = x_label_html,
      y = if (show_y) expression(bold(atop("Predicted", "Niche Overlap"))) else NULL
    ) +
    scale_y_continuous(
      limits = c(0.1, 0.45),
      breaks = seq(0.1, 0.45, by = 0.1),
      expand = c(0, 0)
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.title.x = element_markdown(size = 12, face = "bold"),  # 用于解析 <b>
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold")
    )
  
  return(p)
}

p1 <- plot_gam_effect("mu_phyto_lag4", fit_gam_lagged_all, NLME_lagged_all, "<b>&mu;_Phyto<sub>standardized</sub></b>", TRUE)
p2 <- plot_gam_effect("mu_AOA_lag14", fit_gam_lagged_all, NLME_lagged_all, "<b>&mu;_AOA<sub>standardized</sub></b>", FALSE)
p3 <- plot_gam_effect("par_lag4", fit_gam_lagged_all, NLME_lagged_all, "<b>PAR<sub>standardized</sub></b>", FALSE)
p4 <- plot_gam_effect("NH4_inventory", fit_gam_lagged_all, NLME_lagged_all, "<b>NH<sub>4</sub><sup>+</sup> Inventory<sub>standardized</sub></b>", TRUE)
p5 <- plot_gam_effect("NH4_total_supply_lag34", fit_gam_lagged_all, NLME_lagged_all, "<b>NH<sub>4</sub><sup>+</sup> Supply<sub>standardized</sub></b>", FALSE)

(p1 | p2 | p5) / (p3 | p4 | patchwork::plot_spacer())