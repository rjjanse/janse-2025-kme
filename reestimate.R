# Load packages
library("dplyr",           # Data manipulation
        "magrittr",        # Efficient piping
        "broom",           # Converting statistical objects
        "stringr",         # Working with strings
        "rms",             # Checking assumptions
        "splines",         # Splines
        "cowplot",         # Data viz. add-ons
        "survival",        # Modelling
        "ggplot2",         # Data visualization
        "intsurv"          # Package for calculating C-statistic
)


# Create function
reestimate <- function(data,                     # Data for development and internal validation
                       imputed = TRUE,           # Is the development data imputed? If TRUE, a column .imp should be present indicating the imputation
                       event_label = "dynamic",  # In the histogram with predicted risks per, should the location of the event be dynamic or at a location between 0 and 1
                       annotation = "",          # Should an annotation be added to the plot?
                       save_location = "~/",     # Where should plots be saved?
                       save_name_prefix = ""     # Add-on to add to name of saved plots   
){
    
    # Loading functions
    source("https://raw.githubusercontent.com/rjjanse/hosps/main/3x-functions-20231129.R")
    
    # Saving plots
    quick_ggsave <- function(name, plot){ggsave(name, plot, path = path, width = 5, height = 5, dpi = 600)}
    
    # Create function to get linearity plot
    lin_check <- function(var){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(dat_tmp, .imp == 1)[[var]], y = filter(dat_tmp, .imp == 1)[["tte"]], 
                                                         noprint = TRUE, event = as.numeric(filter(dat_tmp, .imp == 1)[["event"]]), 
                                                         model = "cox", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame()
        
        # Get x-axis label
        xlabel <- ifelse(var == "age", "Age", ifelse(var == "bmi", "BMI", "eGFR"))
        
        # Get plot
        p <- ggplot(dat_plot, aes(x = x, y = V3)) +
            # Geometries
            geom_ribbon(aes(ymin = V4, ymax = V5), alpha = 0.2) +
            geom_line() +
            # Scaling
            scale_y_continuous(expand = expansion(add = 1)) +
            # Labelling
            xlab(xlabel) +
            ylab("Log relative hazard") +
            # Aesthetics
            theme_bw() +
            theme(panel.grid = element_blank())
        
        # Return plot
        return(p)
    }
    
    # Create function to draw proportional hazard plots
    ph_plot <- function(var){
        # Create ylabel
        ylabel = case_when(var == "cardiovasc" ~ "cardiovascular disease",
                           var == "dm" ~ "diabetes mellitus",
                           var == "mal" ~ "malignancy",
                           var == "egfr" ~ "eGFR",
                           var == "bmi" ~ "BMI",
                           .default = var)
        
        # Get data
        dat_plot <- ph_check %>%
            # Keep only relevant variables
            dplyr::select(V1, all_of(var)) %>%
            # Rename second variable
            rename(y = 2)
        
        # Draw plot
        p <- ggplot(dat_plot, aes(x = V1, y = y)) +
            # Geometries
            geom_point(alpha = 0.2) +
            geom_smooth(formula = y ~ x, method = "loess", colour = "black") +
            # Scaling
            scale_x_continuous(breaks = seq(0, max(dat_plot[["V1"]]), 100)) +
            # Labels
            xlab("Time (days)") +
            ylab(bquote(beta * "(t) for " * .(ylabel))) +
            # Aesthetics
            theme_bw() +
            theme(panel.grid = element_blank())
        
        # Return plot
        return(p)
    }
    
    # Save path
    path <- save_location
    
    # If data is not imputed, add artificial imputation column
    if(!imputed) dat_tmp <- mutate(data, .imp = 1) else dat_tmp <- data
    
    ## Check amount of events
    # Get single imputation
    dat_check <- filter(dat_tmp, .imp == 1) %>%
        # Single indicator for censoring or death
        mutate(event2 = ifelse(event == "hospitalization", 1, 0))
    
    # Check lowest event count
    lowest_event_count <- min(table(dat_check[["event2"]]))
    
    # Stop if events or non-events too low (<150)
    if(lowest_event_count < 150) stop(paste0("Not enough events or non-events to re-estimate model"))
    
    ## Check assumptions
    # Collinearity
    coll <- cor(filter(dat_tmp, .imp == 1)[, c("age", "cardiovasc", "dm", "mal", "egfr", "bmi")], method = "spearman") %>%
        # Round correlations
        round(3)
    
    # Save linearity plots for each continuous variable
    invisible(lapply(c("age", "egfr", "bmi"), \(x) quick_ggsave(paste0(save_name_prefix, "lin_", x, ".png"), lin_check(x))))
    
    # Proportional hazards for all predictors
    # Fit full model
    fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + as.factor(cardiovasc) + as.factor(dm) + as.factor(mal) + 
                     egfr + bmi, weight = fgwt, data = finegray(Surv(tte, event) ~ ., data = filter(dat_tmp, .imp == 1)))
    
    # Get data for assumption check
    ph_check <- 
        # Bind data together
        cbind(# Time (x-axis)
            cox.zph(fit)[["time"]],
            # Y per variable
            cox.zph(fit)[["y"]]) %>%
        # Change to data frame
        as.data.frame() %>%
        # Remove factor from column names
        set_colnames(gsub("as.factor\\(|\\)", "", colnames(.)))
    
    # Get proportional hazard plots
    invisible(lapply(c("age", "cardiovasc", "dm", "mal", "egfr", "bmi"), \(x){quick_ggsave(paste0(save_name_prefix, "ph_", x, ".png"), ph_plot(x))}))
    
    # Develop model
    model_vars <- dev(dat_tmp, "Surv(tte, event) ~ age + as.factor(cardiovasc) + as.factor(dm) + as.factor(mal) + egfr + bmi", "fine-gray")
    
    # Get calibration plot
    plot_cal <- quiet(pred(filter(dat_tmp, .imp %in% 1:5), "fine-gray", event, tte, lpsamp = lpsamp(dat_tmp)) %>%
                          # Validate
                          validate(observed, pred, lp, "fine-gray", time, histogram_label = 0.3, deciles = TRUE))
    
    # Add annotation to calibration plot
    plot_cal <- plot_cal + annotate("text", x = 0.15, y = 0.925, label = annotation, colour = "black", fontface = "bold", size = 10)
    
    # Save calibration plot
    quick_ggsave(paste0(save_name_prefix, "calibration.png"), plot_cal)
    
    # Get model performance statistics
    model_perf <- quiet(capture.output(pred(filter(dat_tmp, .imp %in% 1:5), "fine-gray", event, tte, lpsamp = lpsamp(dat_tmp)) %>%
                                           # Validate
                                           validate(observed, pred, lp, "fine-gray", time, histogram_label = event_label, deciles = TRUE)))
    
    # Clean model performance statistics
    model_stats <- str_extract_all(model_perf[[2]], "\\d.\\d{2,3}")[[1]]
    
    # Calibration in the large
    citl <- model_stats[[1]]
    
    # Calbration slope
    cslope <- model_stats[[2]]
    
    # C-statistic
    c <- paste0(model_stats[[3]], " (", model_stats[[4]], "-", model_stats[[5]], ")")
    
    ## Print results
    # Model
    model_print <- model_vars %>%
        # Tranpose
        t() %>%
        # To data frame
        as.data.frame() %>%
        # Add column name
        set_colnames("") %>%
        # Set row names
        set_rownames(c("Baseline hazard", "Age", "Cardiovascular disease", "Diabetes mellitus", "Malignancy", "eGFR", "BMI"))
    
    # Print message for model
    message("")
    message(paste0("The plots for linearity and proportional hazards are available in ", path.expand(save_location)))
    message("The correlation matrix for collinearity is:")
    print(coll)
    message("")
    
    message("The baseline hazard and coefficients of the re-estimated model are: ")
    print.data.frame(model_print)
    message("")
    
    # Formula for calculating individual risks
    message("To calculate individual risks of hospitalization, do the following:")
    message("1. Calculate each individual's linear predictor:")
    message(paste0("LP = age (years) * ", format(round(as.numeric(model_vars[["age"]]), 3), nsmall = 3), 
                   " + presence of cardiovascular disease * ", format(round(as.numeric(model_vars[["as.factor(cardiovasc)1"]]), 3), nsmall = 3),
                   " + presence of diabetes mellitus * ", format(round(as.numeric(model_vars[["as.factor(dm)1"]]), 3), nsmall = 3), 
                   " + presence of malignancy * ", format(round(as.numeric(model_vars[["as.factor(mal)1"]]), 3), nsmall = 3), 
                   " + eGFR (mL/min/1.73m2) * ", format(round(as.numeric(model_vars[["egfr"]]), 3), nsmall = 3), 
                   " + BMI (kg/m2) * ",  format(round(as.numeric(model_vars[["bmi"]]), 3), nsmall = 3)))
    message("")
    
    message("2. Transform the linear predictor into the risk:")
    message(paste0("1 - exp(", format(round(as.numeric(model_vars[["bh"]]), 3), nsmall = 3), " ^ exp(LP - ", format(round(lpsamp(dat_tmp), 3), nsmall = 3), "))"))
    message("")
    
    message("The validation metrics are:")
    print.table(tibble("Calibration-in-the-large" = citl, "Calibration slope" = cslope, "C statistic" = c))
}
