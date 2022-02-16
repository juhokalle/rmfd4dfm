# Preliminaries ####
pkgs <- c("gridExtra") # for nice plots
void <- lapply(pkgs, library, character.only = TRUE)
set.seed(20220111) # for initial values and bootstrapping

# Load data ####
int_vars_fred <- c("INDPRO", "CPIAUCSL", "FEDFUNDS", "EXSZUSx")
FRED_heavy$int_ix <- sapply(int_vars_fred, function(x) which(names(FRED_heavy$df)==x))

# Estimate the number of static factors ####
baing_cr <- baingcriterion(df_fred$df %>% as.matrix, rmax = 25)$IC[1:2]
abc_cr <- replicate(50, abc_crit(df_fred$df %>% as.matrix, kmax = 25))
get_mode <- function(x) unique(x)[which.max(colSums(sapply(unique(x), function(z) x %in% z)))]
r_hats <- c("ICp1" = baing_cr[1],
            "ICp2" = baing_cr[2],
            "ABCp1" = get_mode(unlist(abc_cr[1,])),
            "ABCp2" = get_mode(unlist(abc_cr[2,]))
            )

# Model estimation ####
est_obj0 <- do_everything_rmfd(df = FRED_heavy,
                               r = 8,
                               h = 50,
                               nrep = 0,
                               conv_crit = 1e-3,
                               ci = 0.8,
                               init_rep = 0,
                               verbose = TRUE)
# DFM-FGLR
est_fglr <- do_everything_fglr(df = FRED_heavy,
                               r = 8,
                               k = 2,
                               h = 50,
                               nrep = 500,
                               ci = 0.8)

# SVAR
svar_irf <- do_everything_svar(df = FRED_heavy,
                               p = 9,
                               nrep = 500,
                               h = 50,
                               ci = 0.8)
# IRF plots #####
p1 <- plot_irfs(est_obj$irf, int_vars_fred, "D-DFM")
p2 <- plot_irfs(est_fglr$irf, int_vars_fred, "S-DFM", label_y = FALSE)
p3 <- plot_irfs(svar_irf, int_vars_fred, "SVAR", label_y = FALSE)
plot1 <- marrangeGrob(c(p1,p2,p3), nrow = 4, ncol = 3, top = NULL)
ggsave("plot1.pdf", plot1, width = 15, height = 20, units = "cm")
