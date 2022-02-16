# Comparison between FG and FRED-MD datasets
pkgs <- c("tidyverse", "gridExtra")
void <- lapply(pkgs, library, character.only = TRUE)

# Load data ####
FG_data$int_ix <- c(5, 96, 75, 106)
int_vars_fred <- c("INDPRO", "CPIAUCSL", "FEDFUNDS", "EXSZUSx")
FRED_light$int_ix <- sapply(int_vars_fred, function(x) which(names(FRED_light$df)==x))

# IRF estimation #####
# SVAR: FG data
fg_var <- est_ar(obj = FG_data$df[,FG_data$int_ix],
                 method = 'ols',
                 p.max = 9,
                 ic = "max",
                 mean_estimate = "intercept")
fg_sirf <- fg_var$model$sys %>%
  pseries(lag.max = 50) %>%
  unclass() %>%
  irf_x(post_mat = fg_var$model$sigma_L) %>%
  finalize_irf(shock_size = 0.5,
               norm_id = c(3,3),
               trans_ix = FRED_light$trans_ix[FRED_light$int_ix],
               int_vars = 1:4)

# SVAR: FRED-MD data
fred_var <- est_ar(obj = FRED_light$df[, FRED_light$int_ix],
                   method = 'ols',
                   p.max = 9,
                   ic = "max",
                   mean_estimate = "intercept")
fred_sirf <- fred_var$model$sys %>%
  pseries(lag.max = 50) %>%
  unclass %>%
  irf_x(post_mat = fred_var$model$sigma_L) %>%
  finalize_irf(shock_size = 0.5,
               norm_id = c(3,3),
               trans_ix = FRED_light$trans_ix[FRED_light$int_ix],
               int_vars = 1:4)

colnames(fred_sirf) <- colnames(fg_sirf) <- int_vars_fred
var_irfs <- bind_rows(bind_cols(fred_sirf, lag = 0:50, Dataset = "FRED"),
                      bind_cols(fg_sirf, lag = 0:50, Dataset = "FG")) %>%
  pivot_longer(cols=-c(lag, Dataset))

# DFM-FGLR: FG data
fg_sirf <- DfmRawImp(X = FG_data$df, q = 4, r = 16, k = 2, h = 50) %>%
  chol_ident(int_vars = FG_data$int_ix, est_type = "fglr") %>%
  finalize_irf(shock_size = 0.5,
               norm_id = c(FG_data$int_ix[3], 3),
               int_vars = FG_data$int_ix,
               trans_ix = FG_data$trans_ix)
# DFM-FGLR: FRED data
fred_sirf <- DfmRawImp(X = FRED_light$df, q = 4, r = 16, k = 2, h = 50) %>%
  chol_ident(int_vars = FRED_light$int_ix, est_type = "fglr") %>%
  finalize_irf(shock_size = 0.5,
               norm_id = c(FRED_light$int_ix[3],3),
               int_vars = FRED_light$int_ix,
               trans_ix = FRED_light$trans_ix)
colnames(fred_sirf) <- colnames(fg_sirf) <- int_vars_fred
dfm_irfs <- bind_rows(bind_cols(fred_sirf, lag = 0:50, Dataset = "FRED"),
                      bind_cols(fg_sirf, lag = 0:50, Dataset = "FG")) %>%
  pivot_longer(cols=-c(lag, Dataset))

# Draw IRFs ####
plot_all <- list()
irf_tbl <- bind_rows(bind_cols(var_irfs, model="var"), bind_cols(dfm_irfs, model="dfm"))
for(jj in seq_along(int_vars_fred)){
  plot_all[[jj]] <- irf_tbl %>%
    filter(model == "var", name==int_vars_fred[jj]) %>%
    ggplot(aes(x=lag, y = value, col = Dataset)) +
    geom_line(aes(linetype=Dataset)) +
    scale_color_manual(values=c('blue','red')) +
    labs(title = ifelse(jj==1, "VAR", ""), x="", y=int_vars_fred[jj]) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))

  plot_all[[jj+length(int_vars_fred)]] <- irf_tbl %>%
    filter(model == "dfm", name==int_vars_fred[jj]) %>%
    ggplot(aes(x=lag, y = value, col = Dataset)) +
    geom_line(aes(linetype=Dataset)) +
    scale_color_manual(values=c('blue','red')) +
    labs(title= ifelse(jj==1, "DFM",""), x="", y="") +
    if(jj==length(int_vars_fred)){
      theme(legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            plot.title = element_text(hjust = 0.5))
    } else { theme(legend.position = "none",
                   plot.title = element_text(hjust = 0.5)) }
}

marrangeGrob(plot_all, nrow = 4, ncol = 2, top = NULL)

# Compare ACFs ####
light_acf <- get_acf_qnt(FRED_light$df) %>% round(2)
heavy_acf <- get_acf_qnt(FRED_heavy$df) %>% round(2)
rownames(light_acf) <- rownames(heavy_acf) <- gsub("%", "", rownames(light_acf))

rbind("Percentile" = c("Lag", rep("", 7)),
      rbind(Light = paste(1:8), light_acf),
      rbind(Heavy = paste(1:8), heavy_acf)) -> acf_tbl
colnames(acf_tbl) <- rep("", 8)
print(acf_tbl, quote = FALSE)
