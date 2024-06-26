#--------------------------------------------------------------#
# Prediction models for hospitalization in CKD
# External validation of models
# Roemer J. Janse - 2023/03/22
#--------------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",     # Data manipulation
               "magrittr",  # Efficient piping
               "writexl",   # Writing to .xlsx
               "rms",       # C-statistic
               "cowplot",   # Data viz. add-ons
               "tableone",  # Table one
               "ggplot2"    # Data visualization
)

# Clean Global Environment
rm(list = ls())

# Path based on directionary
if(grepl("rjjanse", getwd())){
    # Set path
    path <- "C:/Users/rjjanse/OneDrive - LUMC/Research/Projects/10. EV_Hosp/"
} else {
    # Set path
    path <- "C:/Users/rjjan/OneDrive/Bureaublad/10. EV_Hosp/"
}

# Set working directionary
setwd(paste0(path, "codes/dataframes/"))

# Width and height for figures
wd = 7; hg = 7

# 1. General validation functions ----
## 1.1. Calibration plot ----
validate <- function(observed, predicted, df, continuous = FALSE, model, axis_lab = "probability", dynamic_histogram_label = FALSE, histogram_label = 0.8,
                     annotation = "", deciles = TRUE, zoom = FALSE){
    # Get observed data
    o <- df[[deparse(substitute(observed))]]
    
    # Get predicted data
    p <- df[[deparse(substitute(predicted))]]
    
    # Combine data into data frame
    dat_tmp <- data.frame(o = o,
                          p = p)
    
    if(deciles){
        # Set seed for jitter
        set.seed(1)
        
        # Create deciles data
        dat_dec <- dat_tmp %>%
            # Create deciles
            mutate(dec = cut(p, breaks = abs(as.numeric(sort(jitter(quantile(p, probs = seq(0, 1, 0.1)), amount = 0.00001)))), include.lowest = TRUE)) %>%
            # Sort for grouping
            arrange(dec) %>%
            # Group per decile
            group_by(dec) %>%
            # Get mean outcome and probability
            mutate(# Outcome
                   out_prop = mean(o),
                   # Predicted
                   pred_prop = mean(p),
                   # Check number of individuals
                   nmax = max(row_number())) %>%
            # Keep one row per decile
            slice(1L) %>%
            # Ungroup again
            ungroup() %>%
            # Keep only relevant columns
            dplyr::select(dec, out_prop, pred_prop, nmax)
    }
    
    # Set max plot limits for continuous
    if(continuous){
        # Upper limit
        scale_max <- pmax(max(p), max(o))
        
        # Scale
        brks <- as.vector(round(quantile(0:scale_max, probs = seq(0, 1, 0.2))))
    }
    
    # Set max plot limits for binary
    else {
        # Upper limit
        scale_max <- 1
        
        # Breaks
        brks <- seq(0, 1, 0.2)
    }
    
    # Dynamic event label in histogram
    if(dynamic_histogram_label){
        # Get location
        loc <- as.numeric(names(sort(table(dat_tmp[["p"]]))[1]))
    }
    
    # Otherwise take specified
    else{
        # Set location
        loc <- histogram_label
    }
    
    # Create base scatter plot
    plot_cal <- ggplot(dat_tmp, aes(x = p, y = o)) +
        # Geometries
        geom_abline(colour = "black", linewidth = 2, alpha = 0.33) +
        #geom_point(inherit.aes = FALSE, data = dat.plot.dec, aes(x = avg_risk, y = avg_rate)) +
        geom_smooth(colour = "darkred", fill = "darkred", method = "loess", formula = y ~ x) +
        geom_point(alpha = 0.25)
    
    # Add deciles of binary
    if(deciles) plot_cal <- plot_cal + geom_point(inherit.aes = FALSE, data = dat_dec, aes(x = pred_prop, y = out_prop), shape = 0)
    
    # Finish plot
    plot_cal <- plot_cal +
        # Scaling
        scale_x_continuous(breaks = brks, name = paste0("Predicted ", axis_lab)) +
        scale_y_continuous(breaks = brks, name = paste0("Observed ", axis_lab)) +
        # Transformations
        coord_cartesian(xlim = c(0, scale_max), ylim = c(0, scale_max)) +
        # Aesthetics
        theme(panel.background = element_blank(),
              plot.background = element_blank(),
              panel.grid = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.line = element_blank()) +
        panel_border(colour = "black", size = 1)
    
    # If continuous, change x axis theme
    if(continuous){
        # Change theme
        plot_cal <- plot_cal +
            # Aesthetics
            theme(axis.ticks.x = element_line(),
                  axis.text.x = element_text(),
                  axis.title.x = element_text())
    }
    
    # If zoom, make zoom plot
    if(length(zoom) > 1){
        # Make zoom plot
        plot_zoom <- plot_cal + 
            # Scaling
            scale_y_continuous(breaks = seq(0, zoom[1], zoom[1] / 5), name = waiver()) +
            scale_x_continuous(breaks = seq(0, zoom[1], zoom[1] / 5), name = waiver()) +
            # Transformatiosn
            coord_cartesian(xlim = c(0, zoom[1]), ylim = c(0, zoom[1])) + 
            # Labels
            ggtitle("Close-up") +
            xlab("") + 
            ylab("") + 
            # Aesthetics
            theme(panel.border = element_rect(colour = "black"),
                  axis.ticks.x = element_line(),
                  axis.text.x = element_text(),
                  panel.background = element_blank(),
                  plot.background = element_blank(),
                  panel.grid = element_blank(),
                  plot.title = element_text(face = "bold", hjust = 0.5))
        
        # Add zoom plot
        plot_cal <- plot_cal + annotation_custom(ggplotGrob(plot_zoom), xmin = zoom[2], xmax = scale_max, ymax = zoom[3])
    }
    
    # Check annotation
    if(length(annotation) > 1){
        plot_cal <- plot_cal + annotate("text", x = as.numeric(annotation[2]), y = as.numeric(annotation[3]), label = annotation[1], fontface = "bold", size = 8)
    }
    
    # Only create histogram plot if not continuous
    if(!continuous){
        # Create histogram of predicted probabilities
        # Warning that two rows with missing values for geom_bar() are removed can be ignored; warning originates from the scale_x_continuous argument limits, but histograms do not actually change
        plot_his <- ggplot(data = dat_tmp) +
            # Geometries
            geom_histogram(data = filter(dat_tmp, o == 0), aes(x = p), binwidth = 0.001, colour = "black", fill = "black") +
            geom_histogram(data = filter(dat_tmp, o == 1), aes(x = p, y = -after_stat(count)), binwidth = 0.001, colour = "black", fill = "black") +
            annotate("text", x = loc, y = (sort(table(dat_tmp[["p"]]))[[1]] + max(table(dat_tmp[["p"]]))) / 2, label = "No event", hjust = 0) +
            annotate("text", x = loc, y = -(sort(table(dat_tmp[["p"]]))[[1]] + max(table(dat_tmp[["p"]]))) / 2, label = "Event", hjust = 0) +
            # Scaling
            scale_x_continuous(breaks = brks, limits = c(0, scale_max)) +
            # Transformations
            coord_cartesian(xlim = c(0, scale_max), ylim = c(-max(table(round(dat_tmp[["p"]], 3)) ), max(table(round(dat_tmp[["p"]], 3))))) +
            # Labels
            xlab(paste0("Predicted ", axis_lab)) + ylab("") +
            # Aesthetics
            theme(panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background = element_rect(fill = "transparent", colour = NA),
                  axis.line = element_blank(),
                  panel.grid = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank()) +
            panel_border(colour = "black", size = 1)
    
        # Combine plots
        plot <- suppressWarnings(plot_grid(plot_cal, plot_his, align = c("hv"), axis = c("tblr"), ncol = 1, rel_heights = c(3, 1)))
    }
        
    # If continuous, make calibration plot the final plot
    if(continuous){
        # Final plot is calibration plot
        plot <- plot_cal
    }
    
    ## Calculate calibration statistics
    # Calibration-in-the-large
    citl <- format(round(mean(o) / mean(p), 3), nsmall = 3)
    
    # Calibration slope
    cslope <- format(round(glm(o ~ p, family = model)[["coefficients"]][["p"]], 3), nsmall = 3)
    
    # C-statistic
    if(!continuous){
        # Calculate C statistic
        c <- rcorr.cens(p, o)[["C Index"]]
        
        # Calculate confidence interval around C statistic
        ci <- quantile(do.call("c", lapply(1:500, \(x){
            # Set seed
            set.seed(x)
            
            # Get random sample numbers
            nrs <- sample.int(length(p), replace = TRUE)
            
            # Get samlpes
            samp_p <- p[nrs]; samp_o <- o[nrs]
            
            # Calculate statistic
            c_bootstrap <- rcorr.cens(samp_p, samp_o)[["C Index"]]
            
            # Return statistic
            return(c_bootstrap)
        })), probs = c(0.025, 0.975))
        
        # Pretty up C statistic
        c <- paste0(format(round(c, 2), nsmall = 2), " (",
                    format(round(ci[[1]], 2), nsmall = 2), "-",
                    format(round(ci[[2]], 2), nsmall = 2), ")")
    
        # Output statistics
        print.table(tibble("Calibration-in-the-large" = citl, "Calibration slope" = cslope, "C statistic" = c))
    }
    
    if(continuous){
        # Output statistics
        print.table(tibble("Calibration-in-the-large" = citl, "Calibration slope" = cslope))
    }

    # Return plot
    return(plot)
}
 
# 2. Fried et al. - Charlson ----
## 2.1. NECOSAD HD ----
# Load data
load("dat_necosad_hd_hosp.Rdata")

# Create validation data
dat_val <- dat_necosad_hd_hosp %>%
    # Set Charlson comorbidity index to 2 if 0 or 1
    mutate(charl = ifelse(charl < 2, 2, charl)) %>%
    # Only first visit
    filter(visit == 1) %>%
    # No missing outcome
    filter(missing_hosp == 0) %>%
    # At least one comorbidity
    filter(com >= 1) %>%
    # Estimate linear predictor based on model
    mutate(lp = cci45 * 0.626 + cci67 * 0.718 + ccig7 * 0.784 * + alb * -0.223 + ifelse(black == 1, 0.270, 0)) 

# Save single imputation for tables
dat_fried_cci_necosad_hd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor (log years accounts for rate instead of count): 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(hosp_rate ~ offset(lp), family = "poisson", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted rates
dat_val <- dat_val %>% 
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted rate
        pred_rate = exp(a + lp)) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_rate = mean(pred_rate)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(hosp_rate, pred_rate, dat_val, continuous = TRUE, model = "poisson", axis_lab = "rate (hospitalizations/year)", zoom = c(20, 195, 195),
         annotation = c("A", -1, 487))

# Save plot
ggsave("fried_cci_necosad_hd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

## 2.2. NECOSAD PD ----
# Load data
load("dat_necosad_pd_hosp.Rdata")

# Create validation data
dat_val <- dat_necosad_pd_hosp %>%
    # Set Charlson comorbidity index to 2 if 0 or 1
    mutate(charl = ifelse(charl < 2, 2, charl)) %>%
    # Only first visit
    filter(visit == 1) %>%
    # No missing outcome
    filter(missing_hosp == 0) %>%
    # At least one comorbidity
    filter(com >= 1) %>%
    # Estimate linear predictor based on model
    mutate(lp = cci45 * 0.626 + cci67 * 0.718 + ccig7 * 0.784 * + alb * -0.223 + ifelse(black == 1, 0.270, 0)) 

# Save single imputation for tables
dat_fried_cci_necosad_pd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor (log years accounts for rate instead of count): 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(hosp_rate ~ offset(lp), family = "poisson", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted rates
dat_val <- dat_val %>% 
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted rate
        pred_rate = exp(a + lp)) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_rate = mean(pred_rate)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(hosp_rate, pred_rate, dat_val, continuous = TRUE, model = poisson, axis_lab = "rate (hospitalizations/year)", zoom = c(20, 45, 45),
         annotation = c("B", -1, 112))

# Save plot
ggsave("fried_cci_necosad_pd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

## 2.3. EQUAL ----
# Load data
load("dat_equal_hosp.Rdata")

# Get relevant data
dat_val <- dat_equal_hosp %>% 
    # At least one comorbidity
    filter(com >= 1) %>%
    # Only first visit
    filter(visit == 1) %>%
    # No missing outcome
    filter(missing_hosp == 0) %>%
    # Estimate linear predictor based on model
    mutate(lp = cci45 * 0.626 + cci67 * 0.718 + ccig7 * 0.784 * + albumin * -0.223 + ifelse(race_short == 2, 0.270, 0)) 

# Save single imputation for tables
dat_fried_cci_equal <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor (log years accounts for rate instead of count): 
# https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression)
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(hosp_rate ~ offset(lp), family = "poisson", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted rates
dat_val <- dat_val %>% 
    # Create new variables
    mutate(# Intercept
           a = a,
           # Predicted rate
           pred_rate = exp(a + lp)) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_rate = mean(pred_rate)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()
                   
# Perform validation
validate(hosp_rate, pred_rate, dat_val, continuous = TRUE, model = poisson, axis_lab = "rate (hospitalizations/year)", zoom = c(9, 48, 48), 
         annotation = c("C", -1, 120.975))

# Save plot
ggsave("fried_cci_equal.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

# 3. Fried et al. - Davies ----
## 3.1. NECOSAD HD ----
# Load data
load("dat_necosad_hd_hosp.Rdata")

# Create validation data
dat_val <- dat_necosad_hd_hosp %>%
    # Only first visit
    filter(visit == 1) %>%
    # No missing outcome
    filter(missing_hosp == 0) %>%
    # At least one comorbidity
    filter(com >= 1) %>%
    # Estimate linear predictor based on model
    mutate(lp = dav1 * 0.3351 + dav2 * 0.875 + davg2 * 1.037 * + alb * -0.21 + ifelse(black == 1, 0.307, 0) + age * -0.007) 

# Save single imputation for tables
dat_fried_davies_necosad_hd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor (log years accounts for rate instead of count): 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(hosp_rate ~ offset(lp), family = "poisson", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted rates
dat_val <- dat_val %>% 
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted rate
        pred_rate = exp(a + lp)) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_rate = mean(pred_rate)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(hosp_rate, pred_rate, dat_val, continuous = TRUE, model = "poisson", axis_lab = "rate (hospitalizations/year)", zoom = c(20, 195, 195),
         annotation = c("A", -1, 487))

# Save plot
ggsave("fried_dav_necosad_hd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

## 3.2. NECOSAD PD ----
# Load data
load("dat_necosad_pd_hosp.Rdata")

# Create validation data
dat_val <- dat_necosad_pd_hosp %>%
    # Only first visit
    filter(visit == 1) %>%
    # No missing outcome
    filter(missing_hosp == 0) %>%
    # At least one comorbidity
    filter(com >= 1) %>%
    # Estimate linear predictor based on model
    mutate(lp = dav1 * 0.3351 + dav2 * 0.875 + davg2 * 1.037 * + alb * -0.21 + ifelse(black == 1, 0.307, 0) + age * -0.007) 

# Save single imputation for tables
dat_fried_davies_necosad_pd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor (log years accounts for rate instead of count): 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(hosp_rate ~ offset(lp), family = "poisson", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted rates
dat_val <- dat_val %>% 
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted rate
        pred_rate = exp(a + lp)) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_rate = mean(pred_rate)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(hosp_rate, pred_rate, dat_val, continuous = TRUE, model = poisson, axis_lab = "rate (hospitalizations/year)", zoom = c(20, 45, 45),
         annotation = c("B", -1, 112))

# Save plot
ggsave("fried_dav_necosad_pd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

# 4. Flythe et al. - Admission ----
## 4.1. NECOSAD HD ----
# Load data
load("dat_necosad_hd_los.Rdata")

# Get relevant data
dat_val <- dat_necosad_hd_los %>%
    # No missing readmission
    filter(missing_readmi == 0) %>%
    # Calculate linear predictor
    mutate(lp = mal * 1.012 + prev_hosp_yr_ge3 * 0.806 + med_cou * 0.732 + sas * 0.56 + non_unc * 0.747 + bp_sys_adm_se110 * 0.756 + bp_sys_adm_ge176 * 0.03) 

# Save single imputation for tables
dat_flythe_necosad_hd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor: 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(readmi ~ offset(lp), family = "binomial", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted probability of readmission
dat_val <- dat_val %>%
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted probability
        pred_prob = 1 / (1 + exp(-(a + lp)))) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_prob = mean(pred_prob)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(readmi, pred_prob, dat_val, model = binomial, continuous = FALSE, histogram_label = 0.3, annotation = c("A", 0, 0.975))

# Save plot
ggsave("flythe_necosad_hd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

## 4.2. NECOSAD PD ----
# Load data
load("dat_necosad_pd_los.Rdata")

# Get relevant data
dat_val <- dat_necosad_pd_los %>%
    # No missing readmission
    filter(missing_readmi == 0) %>%
    # Calculate linear predictor
    mutate(lp = mal * 1.012 + prev_hosp_yr_ge3 * 0.806 + med_cou * 0.732 + sas * 0.56 + non_unc * 0.747 + bp_sys_adm_se110 * 0.756 + bp_sys_adm_ge176 * 0.03) 

# Save single imputation for tables
dat_flythe_necosad_pd <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor: 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(readmi ~ offset(lp), family = "binomial", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted probability of readmission
dat_val <- dat_val %>%
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted probability
        pred_prob = 1 / (1 + exp(-(a + lp)))) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_prob = mean(pred_prob)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(readmi, pred_prob, dat_val, model = binomial, continuous = FALSE, histogram_label = 0.3, annotation = c("B", 0, 0.975))

# Save plot
ggsave("flythe_necosad_pd.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

## 4.3. EQUAL ----
# Load data
load("dat_equal_los.Rdata")

# Get relevant data
dat_val <- dat_equal_los %>%
    # No missing readmission
    filter(missing_readmi == 0) %>%
    # Calculate linear predictor
    mutate(lp = mal * 1.012 + prev_hosp_yr_ge3 * 0.806 + med_at * 0.732 + sas * 0.56 + non_unc * 0.747 + bp_sys_adm_se110 * 0.756 + bp_sys_adm_ge176 * 0.03) 

# Save single imputation for tables
dat_flythe_equal <- filter(dat_val, .imp == 1)

# Re-estimate model intercept with offset on linear predictor: 
# Do this per imputation
a <- mean(do.call("c", lapply(1:20, \(x){glm(readmi ~ offset(lp), family = "binomial", data = filter(dat_val, .imp == x))[["coefficients"]][["(Intercept)"]]})))

# Calculate predicted probability of readmission
dat_val <- dat_val %>%
    # Create new variables
    mutate(# Intercept
        a = a,
        # Predicted probability
        pred_prob = 1 / (1 + exp(-(a + lp)))) %>%
    # Sort for grouping
    arrange(studynr) %>%
    # Group on studynr
    group_by(studynr) %>%
    # Take mean of predicted value
    mutate(pred_prob = mean(pred_prob)) %>%
    # Keep one row per person
    slice(1L) %>%
    # Ungroup again
    ungroup()

# Perform validation
validate(readmi, pred_prob, dat_val, model = binomial, continuous = FALSE, histogram_label = 0.8, annotation = c("C", 0, 0.975))

# Save plot
ggsave("flythe_equal.png", path = paste0(path, "/figures/validations"), width = wd, height = hg, dpi = 600)

# 5. Create comparison tables for each model ----
## 5.1. Fried et al. - Charlson ----
# Combine data
dat_tab <- rbind(# NECOSAD HD
                 mutate(dat_fried_cci_necosad_hd, group = "nechd") %>%
                     # Keep only relevant variables
                     select(group, age, dm, female, black, charl, alb, bmi),
                 # NECOSAD PD
                 mutate(dat_fried_cci_necosad_pd, group = "necpd") %>%
                     # Keep only relevant variables
                     select(group, age, dm, female, black, charl, alb, bmi),
                 # EQUAL
                 mutate(dat_fried_cci_equal, group = "equal") %>%
                     # Keep only relevant variables
                     select(group, age, dm, female, race_short, charlson_baseline, albumin, bmi) %>%
                     # Rename variables
                     rename(black = race_short, charl = charlson_baseline, alb = albumin) %>%
                     # Change black
                     mutate(black = ifelse(black == 2, 1, 0))) %>%
    # Add Charlson comorbidity categories
    mutate(charl_cat = case_match(charl,
                                  2:3 ~ "2-3",
                                  4:5 ~ "4-5",
                                  6:7 ~ "6-7",
                                  .default = ">7"),
           # Change to factor
           charl_cat = factor(charl_cat, levels = c("2-3", "4-5", "6-7", ">7")))

# Create table
tab_fried_cci <- dat_tab %>%
    # Create table frame
    # Note that we get IQR here for CCI but we need range as this is what the paper uses
    CreateTableOne(vars = c("age", "female", "black", "dm", "alb", "charl", "charl_cat", "bmi"), strata = "group",
                   factorVars = c("female", "black", "dm", "charl_cat")) %>%
    # Print table one
    print(printToggle = FALSE, nonnormal = "charl") %>%
    # Change to matrix
    as.matrix() %>%
    # Keep only columns of interest
    extract(, c(2, 3, 1)) %>%
    # Change column names
    set_colnames(c("NECOSAD HD", "NECOSAD PD", "EQUAL")) %>%
    # Change to data frame
    as.data.frame() %>%
    # Add column with row names
    mutate(var = rownames(.)) %>%
    # Put in front
    relocate(var, .before = "NECOSAD HD")

# Set IQR for CCI to range
tab_fried_cci[7, 2:4] <- c(# NECOSAD HD
                           paste0(summary(filter(dat_tab, group == "nechd")[["charl"]])[[3]], " [",
                                  summary(filter(dat_tab, group == "nechd")[["charl"]])[[1]], ", ",
                                  summary(filter(dat_tab, group == "nechd")[["charl"]])[[6]], "]"),
                           # NECOSAD PD
                           paste0(summary(filter(dat_tab, group == "necpd")[["charl"]])[[3]], " [",
                                  summary(filter(dat_tab, group == "necpd")[["charl"]])[[1]], ", ",
                                  summary(filter(dat_tab, group == "necpd")[["charl"]])[[6]], "]"),
                           # EQUAL
                           paste0(summary(filter(dat_tab, group == "equal")[["charl"]])[[3]], " [",
                                  summary(filter(dat_tab, group == "equal")[["charl"]])[[1]], ", ",
                                  summary(filter(dat_tab, group == "equal")[["charl"]])[[6]], "]")) 

# Save table to Excel
write_xlsx(tab_fried_cci, path = paste0(gsub("dataframes", "results", getwd()), "/tab_fried_cci.xlsx"))

## 5.2. Fried et al. - Davies ----
# Combine data
dat_tab <- rbind(# NECOSAD HD
                 mutate(dat_fried_davies_necosad_hd, group = "nechd") %>%
                     # Keep only relevant variables
                     select(group, age, dm, female, black, dav_tot, alb, bmi),
                 # NECOSAD PD
                 mutate(dat_fried_davies_necosad_pd, group = "necpd") %>%
                     # Keep only relevant variables
                     select(group, age, dm, female, black, dav_tot, alb, bmi)) %>%
    # Add Davies categories
    mutate(dav_cat = case_match(dav_tot,
                                0 ~ "0",
                                1 ~ "1",
                                2 ~ "2",
                                .default = ">2"),
           # Change to factor
           dav_cat = factor(dav_cat, levels = c("0", "1", "2", ">2")))

# Create table
tab_fried_dav <- dat_tab %>%
    # Create table frame
    # Note that we get IQR here for CCI but we need range as this is what the paper uses
    CreateTableOne(vars = c("age", "female", "black", "dm", "alb", "dav_tot", "dav_cat", "bmi"), strata = "group",
                   factorVars = c("female", "black", "dm", "dav_cat")) %>%
    # Print table one
    print(printToggle = FALSE, nonnormal = "dav_tot") %>%
    # Change to matrix
    as.matrix() %>%
    # Keep only columns of interest
    extract(, 1:2) %>%
    # Change column names
    set_colnames(c("NECOSAD HD", "NECOSAD PD")) %>%
    # Change to data frame
    as.data.frame() %>%
    # Add column with row names
    mutate(var = rownames(.)) %>%
    # Put in front
    relocate(var, .before = "NECOSAD HD")

# Set IQR for CCI to range
tab_fried_dav[7, 2:3] <- c(# NECOSAD HD
                           paste0(summary(filter(dat_tab, group == "nechd")[["dav_tot"]])[[3]], " [",
                                  summary(filter(dat_tab, group == "nechd")[["dav_tot"]])[[1]], ", ",
                                  summary(filter(dat_tab, group == "nechd")[["dav_tot"]])[[6]], "]"),
                           # NECOSAD PD
                           paste0(summary(filter(dat_tab, group == "necpd")[["dav_tot"]])[[3]], " [",
                                  summary(filter(dat_tab, group == "necpd")[["dav_tot"]])[[1]], ", ",
                                  summary(filter(dat_tab, group == "necpd")[["dav_tot"]])[[6]], "]")) 

# Save table to Excel
write_xlsx(tab_fried_dav, path = paste0(gsub("dataframes", "results", getwd()), "/tab_fried_dav.xlsx"))

## 5.3. Flythe et al. - Admission ----
# Combine data
dat_tab <- rbind(# NECOSAD HD
    mutate(dat_flythe_necosad_hd, group = "nechd") %>%
        # Keep only relevant variables
        select(group, age, female, black, dm, hf, cad, hypert, mal, prev_hosp_yr_ge3, med_cou, sas, bp_sys_adm_se110, bp_sys_adm_ge176) %>%
        # Add empty cerebrovascular disease (CVD) and antithrombotics
        mutate(cvd = NA, med_at = NA),
    # NECOSAD PD
    mutate(dat_flythe_necosad_pd, group = "necpd") %>%
        # Keep only relevant variables
        select(group, age, female, black, dm, hf, cad, hypert, mal, prev_hosp_yr_ge3, med_cou, sas, bp_sys_adm_se110, bp_sys_adm_ge176) %>%
        # Add empty cerebrovascular disease (CVD) and antithrombotics
        mutate(cvd = NA, med_at = NA),
    # EQUAL
    mutate(dat_flythe_equal, group = "equal") %>%
        # Keep only relevant variables
        select(group, age, female, race_short, dm, hf, cad, hypert, mal, cvd, prev_hosp_yr_ge3, med_at, sas, bp_sys_adm_se110, bp_sys_adm_ge176) %>%
        # Rename variables
        rename(black = race_short) %>%
        # Change black and add empty coumarins
        mutate(black = ifelse(black == 2, 1, 0), med_cou = NA)) %>%
    # Add age categories and categorize admission blood pressure
    mutate(# Create age categories
           age_cat = case_when(age <= 49 ~ "<=49",
                               age > 49 & age < 60 ~ "50-59",
                               age > 60 & age < 70 ~ "60-69",
                               age >= 70 ~ ">=70"),
           # Change to factor
           age_cat = factor(age_cat, levels = c("<=49", "50-59", "60-69", ">=70")),
           # Create categories for blood pressure at admission
           bp_sys_adm = case_when(bp_sys_adm_se110 == 1 ~ "<=110",
                                  bp_sys_adm_ge176 == 1 ~ ">=176",
                                  .default = "111-175"),
           # Change to factor
           bp_sys_adm = factor(bp_sys_adm, levels = c("<=110", "111-175", ">=176")))

# Set vars
vars = c("age_cat", "female", "black", "dm", "hf", "cad", "hypert", "cvd", "mal", "prev_hosp_yr_ge3", "med_at", "med_cou", "sas", "bp_sys_adm")

# Create table
tab_flythe <- dat_tab %>%
    # Create table frame
    CreateTableOne(vars = vars, strata = "group", factorVars = vars) %>%
    # Print table one
    print(printToggle = FALSE) %>%
    # Change to matrix
    as.matrix() %>%
    # Keep only columns of interest
    extract(, c(2, 3, 1)) %>%
    # Change column names
    set_colnames(c("NECOSAD HD", "NECOSAD PD", "EQUAL")) %>%
    # Change to data frame
    as.data.frame() %>%
    # Add column with row names
    mutate(var = rownames(.)) %>%
    # Put in front
    relocate(var, .before = "NECOSAD HD")

# Save table to Excel
write_xlsx(tab_flythe, path = paste0(gsub("dataframes", "results", getwd()), "/tab_flythe.xlsx"))
