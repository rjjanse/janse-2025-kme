#--------------------------------------------------------------#
# Prediction models for hospitalization in CKD
# Creating functions for development script
# Roemer J. Janse - 2023/11/29
#--------------------------------------------------------------#

# Function for quick save
quick_ggsave <- function(name, plot){ggsave(name, plot, path = paste0(path, "/figures"), width = 5, height = 5, dpi = 600)}

# Function for quiet output
quiet <- function(x){ 
    # Put data in a temporary file
    sink(tempfile()) 
    
    # Sink again on end of function
    on.exit(sink()) 
    
    # Forcefully make output invisible
    invisible(force(x)) 
} 

# Function for linearity check
lin_check <- function(var, annotation, model){
    # Get event
    if(is.factor(dat_dev[["event"]])) ev <- as.numeric(filter(dat_dev, .imp == 1)[["event"]]) - 1 
    else ev <- as.numeric(filter(dat_dev, .imp == 1)[["event"]])
    
    # Fine-Gray model
    if(model == "fine-gray"){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(dat_dev, .imp == 1)[[var]], y = filter(dat_dev, .imp == 1)[["tte"]], noprint = TRUE,
                                                         event = ev, model = "cox", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Cox and AFT model
    if(model %in% c("cox", "aft")){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(dat_dev, .imp == 1)[[var]], y = filter(dat_dev, .imp == 1)[["tte"]], noprint = TRUE,
                                                         event = ev == 1, model = "cox", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Logistic model
    if(model == "logistic"){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(dat_dev, .imp == 1)[[var]], y = ev == 1, 
                                                         noprint = TRUE, model = "logistic", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Get x-axis label
    xlabel <- ifelse(var == "age", "Age (years)", ifelse(var == "bmi", expression("BMI (kg/m" ^ 2 * ")"), 
                                                         expression("eGFR (mL/min/1.73m" ^ 2 * ")")))
    
    # Get plot
    p <- ggplot(dat_plot, aes(x = x, y = V3)) +
        # Geometries
        geom_ribbon(aes(ymin = V4, ymax = V5), alpha = 0.2) +
        geom_line() +
        # Scaling
        scale_y_continuous(limits = c(min(dat_plot[["V4"]]) - 1, max(dat_plot[["V5"]]) + 1), 
                           name = ifelse(model == "logistic", "Log odds", "Log relative hazard")) +
        # Labelling
        xlab(xlabel) +
        annotate("text", x = as.numeric(annotation[2]), y = as.numeric(annotation[3]), label = annotation[1], fontface = "bold", size = 8) +
        # Aesthetics
        theme_bw() +
        theme(panel.grid = element_blank())
    
    # Return plot
    return(p)
}

# Function for proportional hazards plot
ph_plot <- function(var, ylabel, annotation){
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
        ylab(bquote(beta[t] *" for " * .(ylabel))) +
        annotate("text", x = as.numeric(annotation[2]), y = as.numeric(annotation[3]), label = annotation[1], fontface = "bold", size = 8) +
        # Aesthetics
        theme_bw() +
        theme(panel.grid = element_blank())
    
    # Return plot
    return(p)
}

# Function to create model per imputation
develop <- function(df, imp, formula, model, aft_dist){
    # Data
    dat_mod_tmp <- filter(df, .imp == imp) 
    
    # Fine-Gray model fit
    if(model == "fine-gray"){
        # Get left-hand-side of formula
        lhs <- str_extract(formula, ".*(?=~)")
        
        # Get right hand side of formula
        rhs <- str_extract(formula, "(?<=~).*")
        
        # Fit model
        fit <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs)), x = TRUE, 
                     weight = fgwt, data = finegray(as.formula(paste0(lhs, "~ .")), data = dat_mod_tmp))
    }
    
    # Cox model fit
    if(model == "cox") fit <- coxph(as.formula(formula), x = TRUE, data = dat_mod_tmp)
    
    # Accelerated failure time model fit
    if(model == "aft") fit <- survreg(as.formula(formula), data = dat_mod_tmp, dist = aft_dist)

    # Logistic model fit
    if(model == "logistic") fit <- glm(as.formula(formula), family = "binomial", data = dat_mod_tmp)
    
    # If survival model, get baseline hazard and add to coefficients
    if(model %in% c("cox", "fine-gray")){
        # Get baseline hazard for 1 year
        bh <- basehaz(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Remove everything after a year
            filter(time <= 365.25) %>%
            # Keep only last observation
            slice_tail(n = 1L) %>%
            # Keep baseline hazard
            extract2("hazard")
        
        # Get coefficients
        coeffs <- fit[["coefficients"]] %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose
            t() %>%
            # Add baseline hazard
            bind_cols(as.data.frame(bh), .)
    }
    
    # If logistic or accelerated failure time model, get intercept
    if(model %in% c("logistic", "aft")){
        # Get coefficients
        coeffs <- fit[["coefficients"]] %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose
            t() %>%
            # Rename intercept
            set_colnames(c("a", colnames(.)[2:length(colnames(.))]))
    }
    
    # Return model data
    return(coeffs)
}

# Function for model development
dev <- function(df, formula, model, aft_dist = NULL){
    # Fit model for each imputation and add together the final models
    mod <- do.call("rbind", lapply(unique(df[[".imp"]]), \(x) develop(df, x, formula, model, aft_dist)))
    
    # Get final model fit by taking means of each model
    model_vars <- summary(mod) %>% 
        # Change to data frame
        as.data.frame() %>% 
        # Keep only mean
        filter(grepl("Mean", Freq)) %>% 
        # Change parameters to numeric
        mutate(param = as.numeric(gsub("\\D{4} *:", "", Freq))) %>%
        # Rename Var2
        rename(var = Var2) %>%
        # Keep only columns of interest
        dplyr::select(var, param) %>%
        # Transpose
        t() %>%
        # Set column names
        set_colnames(gsub(" ", "", .[1, ])) %>%
        # Change to data frame
        as.data.frame() %>%
        # Keep only second row
        filter(row_number() == 2) 
    
    # Return model vars
    return(model_vars)
}

# C-statistic function
cstat <- function(df, model){
    # Wolbers' adaptation for competing events
    if(model == "fine-gray") dat_tmp <- mutate(df, tim = ifelse(obs >= 2, Inf, tim)) else dat_tmp <- df
    
    # For logistic models
    if(model == "logistic"){
        # All predictions of individuals with event
        events <- dat_tmp %>%
            # Keep only events
            filter(obs == 1) %>%
            # Extract predictions
            extract2("prd")
        
        # All predictions of individuals without event
        no_events <- dat_tmp %>%
            # Keep only events
            filter(obs == 0) %>%
            # Extract predictions
            extract2("prd")
        
        # Compare all events with all non-events and calculate C-statistic
        cstatistic <- sum((tmp <- do.call("c", lapply(events, \(x) ifelse(x > no_events, 1, ifelse(x == no_events, 0.5, 0)))))) / length(tmp); rm(tmp)
    }
    
    # For survival models (Harrell's C-statistic)
    if(model %in% c("cox", "fine-gray")) cstatistic <- cIndex(dat_tmp[["tim"]], dat_tmp[["obs_ncr"]], dat_tmp[["prd"]])[["index"]]
    
    # For AFT models, use aft_status obs
    if(model == "aft") cstatistic <- cIndex(dat_tmp[["tim"]], dat_tmp[["aft_status"]], dat_tmp[["prd"]])[["index"]]
    
    # Return C-statistic
    return(cstatistic)
}

# Function for validation
validate <- function(.data,                                     # Data
                     observed,                                  # Observed outcome
                     predicted,                                 # Predicted outcome
                     lp,                                        # Linear predictor
                     model,                                     # Regression model used to create the prediction model,
                     time = NULL,                               # Time variable (only relevant for Cox/AFT/Fine-Gray)
                     aft_dist = "lognormal",                    # Distribution used for the AFT model
                     # Output
                     print_stratified = FALSE,                  # Print stratified C-statistics (only relevant for survival model)
                     plot = TRUE,                               # Should a calibration plot be made
                     deciles = TRUE,                            # Should deciles be added in calibration plot
                     # Calibration plot details
                     unit = "probability",                      # Unit of prediction for axes of plot
                     annotation = c("", 0, 1),                  # Annotation to add to plot as c("annotation", x, y)
                     smooth_colour = "darkred",                 # Colour of smoother
                     histogram_label = NULL                     # Location of event / no-event label in probabilities histogram
){
    ## Prepare data
    # Observed
    obs <- .data[[deparse(substitute(observed))]]
    
    # Set observed to numeric if it is factor
    if(is.factor(obs)) obs <- as.numeric(obs) - 1
    
    # Observed without competing risks
    if(model == "fine-gray") obs_ncr <- ifelse(obs == 1, 1, 0) else obs_ncr <- obs
    
    # Predicted
    prd <- .data[[deparse(substitute(predicted))]]
    
    # Linear predictors
    lps <- .data[[deparse(substitute(lp))]] 
    
    # Time
    if(!(deparse(substitute(time)) == "NULL")) tim <- .data[[deparse(substitute(time))]] else tim <- NA
    
    # If model is AFT, observed should be time to event and status should be stored elsewhere
    if(model == "aft") {
        # Store status elsewhere
        aft_status <- obs
        
        # Store EFT
        obs <- obs_ncr <- tim
    }
    
    # Else empty aft_status variable
    else aft_status <- NA
    
    # Create data to work with
    dat <- tibble(obs = obs,
                  obs_ncr = obs_ncr,
                  lps = lps,
                  prd = prd, 
                  tim = tim,
                  aft_status = aft_status)
    
    # Get total number of individuals
    n <- nrow(dat)
    
    ## If plotting, prepare data and make plot
    if(plot){
        # Calculate deciles
        if(deciles){
            # Create deciles data
            dat_dec <- dat %>%
                # Create decile
                mutate(dec = cut(prd, breaks = as.numeric(quantile(prd, probs = seq(0, 1, 0.1))), include.lowest = TRUE)) %>%
                # Sort for grouping
                arrange(dec) %>%
                # Group per decile
                group_by(dec) %>%
                # Get mean outcome and probability
                mutate(# Outcome
                       out_prop = mean(obs_ncr),
                       # Predicted
                       pred_prop = mean(prd),
                       # Check number of individuals
                       nmax = max(row_number())) %>%
                # Keep one row per decile
                slice(1L) %>%
                # Ungroup again
                ungroup() %>%
                # Keep only relevant columns
                select(dec, out_prop, pred_prop, nmax)
        }
        
        # Calculate pseudo-observations
        if(model == "fine-gray"){
            # Create empty id column
            dat <- mutate(dat, ids = 1:nrow(dat))
            
            # Get A-J estimate for outcome in total population
            aj_tot <- survfit(Surv(tim, obs, type = "mstate") ~ 1, data = dat) %>%
                # Change fit list to dataframe
                tidy() %>%
                # Keep only events
                filter(state == 1) %>%
                # Keep only last observation, assuming this is the time point of interest
                slice_tail(n = 1L) %>%
                # Keep baseline hazard
                extract2("estimate")
            
            # Get A-J estimates for excluding each individual with jackknife
            aj_in <- do.call("c", lapply(dat[["ids"]], \(x){
                # Create new data
                dat_tmp <- filter(dat, ids != x)
                
                # Calculate A-J estimate
                aj <- survfit(Surv(tim, obs, type = "mstate") ~ 1, data = dat_tmp) %>%
                    # Change fit list to dataframe
                    tidy() %>%
                    # Keep only events
                    filter(state == 1) %>%
                    # Keep only last observation, assuming this is the time point of interest
                    slice_tail(n = 1L) %>%
                    # Keep baseline hazard
                    extract2("estimate")
                
                # Return estimate
                return(aj)
            }))
            
            # Get jackknife A-J estimate per individual
            dat <- dat %>%
                # Calculate new variables
                mutate(# Total A-J estimate
                    aj_t = aj_tot,
                    # Jackknife A-J estimate
                    aj_i = aj_in,
                    # Pseudovalue for the outcome
                    aj_o = (n * aj_t) - ((n - 1) * aj_i))
        }
        
        # Set y variable based on model
        if(model == "fine-gray") dat <- mutate(dat, y = aj_o) else dat <- mutate(dat, y = obs)
        
        ## Define characteristics of the plot
        # Upper limits of axes
        if(model %in% c("poisson", "linear", "aft")) xmax <- ymax <- max(obs) else xmax <- ymax <- 1
        
        # Lower limits of axes
        if(model %in% c("poisson", "linear", "aft")) xmin <- ymin <- pmin(0, min(obs)) else xmin <- ymin <- 0
        
        # Breaks of axes
        brks <- seq(xmin, xmax, (xmax - xmin) / 10)
        
        # Labels of x-axis
        if(model %in% c("poisson", "linear", "aft")) xlab <- paste0("Predicted ", unit) else xlab <- "Predicted probability"
        
        # Label of y-axis
        if(model %in% c("poisson", "linear", "aft")) ylab <- paste0("Observed ", unit) else ylab <- "Observed probability"
        
        # Determine location of event label in probabilities histogram automatically if not specified, based on least predicted probability
        if(is.null(histogram_label)) histogram_label <- as.numeric(names(sort(table(round(dat[["prd"]], 3))))[[1]])
        
        ## Create plot
        # Create base scatter plot
        plot_cal <- ggplot(dat, aes(x = prd, y = y)) +
            # Geometries
            geom_abline(colour = "black", linewidth = 2, alpha = 0.33) +
            geom_point(alpha = 0.25) +
            geom_smooth(colour = smooth_colour, fill = smooth_colour, method = "loess", formula = y ~ x)
        
        # If AFT model, overwrite base scatter plot to add colouring for different statuses
        if(model == "aft"){
            # Create base scatter plot
            plot_cal <- ggplot(dat, aes(x = prd, y = y)) +
                # Geometries
                geom_abline(colour = "black", linewidth = 2, alpha = 0.33) +
                geom_point(alpha = 0.25, mapping = aes(colour = factor(aft_status, levels = c(0, 1), labels = c("No-event", "Event")))) +
                geom_smooth(colour = smooth_colour, fill = smooth_colour, method = "loess", formula = y ~ x) + 
                # Scaling
                scale_colour_manual(values = c("darkorange", "darkred")) +
                # Aesthetics
                theme(legend.position = "bottom",
                      legend.title = element_blank())
        }
        
        # Add deciles
        if(deciles) plot_cal <- plot_cal + geom_point(inherit.aes = FALSE, data = dat_dec, aes(x = pred_prop, y = out_prop), shape = 0)
        
        # Finish calibration plot
        plot_cal <- plot_cal +
            # Scaling
            scale_y_continuous(breaks = brks, name = ylab) +
            scale_x_continuous(breaks = brks, name = xlab) +
            # Transformations
            coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            # Labels
            annotate("text", x = as.numeric(annotation[[2]]), y = as.numeric(annotation[[3]]), 
                     label = annotation[[1]], fontface = "bold", size = 10) +
            # Aesthetics
            theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank()) +
            panel_border(colour = "black", size = 1) 
        
        # Create probability histogram for non-continuous models
        if(!(model %in% c("poisson", "linear", "aft"))){
            # Turn off x-axis for calibration plot
            plot_cal <- plot_cal + theme(axis.ticks.x = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.title.x = element_blank())
            
            # Create histogram of predicted probabilities
            plot_his <- ggplot(data = dat) +
                # Geometries
                geom_histogram(data = filter(dat, obs_ncr == 1), aes(x = prd), binwidth = 0.001, colour = "black", fill = "black") +
                geom_histogram(data = filter(dat, obs_ncr == 0), aes(x = prd, y = -after_stat(count)), binwidth = 0.001, colour = "black", 
                               fill = "black") +
                annotate("text", x = histogram_label, y = max(table(round(dat[["prd"]], 3))) / 2, label = "Event", hjust = 0) +
                annotate("text", x = histogram_label, y = -max(table(round(dat[["prd"]], 3))) / 2, label = "No event", hjust = 0) +
                # Scaling
                scale_x_continuous(breaks = brks, name = xlab, limits = c(xmin, xmax)) +
                # Transformations
                coord_cartesian(xlim = c(xmin, xmax), ylim = c(-max(table(round(dat[["prd"]], 3))), max(table(round(dat[["prd"]], 3))))) +
                # Aesthetics
                theme(panel.background = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
                panel_border(colour = "black", size = 1)
            
            # Combine plots
            plot_cal <- suppressWarnings(plot_grid(plot_cal, plot_his, align = c("hv"), axis = c("tblr"), ncol = 1, rel_heights = c(3, 1)))
        }
    }
    
    ## Calculate performance measures
    # Get proportion of outcome for non-continuous models
    if(!(model %in% c("linear", "poisson", "aft", "cox", "fine-gray"))) prop_out <- prop.table(table(dat[["obs"]]))[["1"]]
    
    # Get mean of outcome for continous models
    if(model %in% c("linear", "poisson")) prop_out <- mean(obs)
    
    # Calibration-in-the-large
    if(model %in% c("linear", "poisson", "logistic")) citl <- format(round(prop_out - mean(prd), 3), nsmall = 3)
    
    # Calibation-in-the-large with cumulative incidence function for survival models
    if(model %in% c("cox", "fine-gray")){
        citl <- 
            (# Fit cumulative incidence function
             survfit(Surv(tim, obs, type = "mstate") ~ 1, data = dat) %>%
                 # Change fit list to dataframe
                 tidy() %>%
                 # Keep only events
                 filter(state == 1) %>%
                 # Keep only last observation, assuming this is the time point of interest
                 slice_tail(n = 1L) %>%
                 # Keep baseline hazard
                 extract2("estimate")) %>%
            # Divide by mean predicted risk
            subtract(mean(dat[["prd"]])) %>%
            # Round
            round(digits = 3) %>%
            # Format
            format(nsmall = 3)
    }

    ## Calculate calibration slope
    # Generalized linear model
    if(!(model %in% c("cox", "fine-gray", "aft"))) cslope <- format(round(glm(y ~ lps, family = ifelse(model == "logistic", "binomial", model), 
                                                                            data = dat)[["coefficients"]][["lps"]], 3), nsmall = 3)
    
    # Cox model
    if(model == "cox") cslope <- format(round(coxph(Surv(tim, obs) ~ lps, data = dat)[["coefficients"]][["lps"]], 3), nsmall = 3)
    
    # Fine-Gray model
    if(model == "fine-gray") cslope <- format(round(coxph(Surv(fgstart, fgstop, fgstatus) ~ lps, weight = fgwt, 
                                                          data = finegray(Surv(tim, as.factor(obs)) ~ ., data = dat))[["coefficients"]][["lps"]], 
                                                    3), nsmall = 3)
    
    # AFT model
    if(model == "aft") cslope <- format(round(survreg(Surv(tim, aft_status) ~ lps, data = dat, dist = aft_dist)[["coefficients"]][["lps"]], 3), nsmall = 3)

    # C-statistic
    if(!(model %in% c("linear", "poisson"))){
        # Calculate C statistic
        c <- cstat(dat, model)
        
        # Calculate confidence interval around C statistic
        ci <- quantile(do.call("c", lapply(1:500, \(x){
            # Set seed
            set.seed(x)
            
            # Get random sample numbers
            nrs <- sample.int(length(prd), replace = TRUE)
            
            # Get samples
            dat_samp <- dat[nrs, ] %>%
                # Change studynr
                mutate(studynr = 1:nrow(dat))
            
            # Calculate statistic
            c_bootstrap <- cstat(dat_samp, model)
            
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
    
    # Output statistics
    if(model %in% c("linear", "poisson")) print.table(tibble("Calibration-in-the-large" = citl, "Calibration slope" = cslope))
    
    # Return plot
    if(plot) return(plot_cal)
}

# Create function to get model information
model_inf <- function(){
    # Get factor info if factors were used
    if(length(names(model_vars[grepl("factor\\(", names(model_vars))])) > 0){
        # Get information on factor variables
        factor_info <- names(model_vars[grepl("factor\\(", names(model_vars))]) %>%
            # Change into data frame
            as.data.frame() %>%
            # Arrange for grouping
            arrange(.) %>%
            # Group per variable
            group_by(.) %>%
            # Create new columns
            mutate(# Variable
                var = str_replace_all(., "(as.)+factor\\(|\\)\\d*", ""),
                # Levels
                levels = as.numeric(str_extract(., "(?<=\\))\\d+")),
                # Rows per group
                rows = max(row_number())) %>%
            # Ungroup again
            ungroup()
        
        # Add information to spline_info
        factor_info <- mutate(factor_info, 
                              # Coefficient per level
                              coef = as.numeric(as.vector(model_vars[grepl("factor", names(model_vars))])),
                              # Add type
                              type = "factor") %>%
            # Keep only relevant variables
            select(var, levels, coef, type) %>%
            # Create empty lower and upper level columns
            mutate(lower_level = -Inf, upper_level = Inf)
    }
    
    # Else, create empty data
    else factor_info <- tibble(var = as.character(),
                               coef = as.numeric(),
                               lower_level = as.numeric(), 
                               upper_level = as.numeric(),
                               type = as.character())
    
    # Get spline info if splines were used
    if(length(names(model_vars[grepl("knot", names(model_vars))])) > 0){
        # Get information on splined variables
        spline_info <- names(model_vars[grepl("knot", names(model_vars))]) %>%
            # Change into data frame
            as.data.frame() %>%
            # Arrange for grouping
            arrange(.) %>%
            # Group per variable
            group_by(.) %>%
            # Create new columns
            mutate(# Variable
                var = str_replace_all(., ".*(?<!,knot=c)\\(|\\)\\d|,knot.*", ""),
                # Knots
                knots = paste0("-Inf,", str_extract(., "(?<=[\\(|=])\\d+(,\\d*)*"), ",Inf"),
                # Levels
                levels = as.numeric(str_extract(., "(?<=\\))\\d+")),
                # Rows per group
                rows = max(row_number())) %>%
            # Ungroup again
            ungroup()
        
        # Matrix of spline levels
        spline_levels <- str_split(spline_info[["knots"]], ",", simplify = TRUE)
        
        # Add information to spline_info
        spline_info <- mutate(spline_info, 
                              # Lower level of section
                              lower_level = as.numeric(spline_levels[rows[[1]], levels]),
                              # Upper level of section
                              upper_level = as.numeric(spline_levels[rows[[1]], levels + 1]),
                              # Coefficient per level
                              coef = as.numeric(as.vector(model_vars[grepl("knot", names(model_vars))])),
                              # Add type
                              type = "spline") %>%
            # Keep only relevant variables
            select(var, lower_level, upper_level, coef, type) %>%
            # Create empty levels variable
            mutate(levels = NA)
    }
    
    # Else, create empty data
    else spline_info <- tibble(var = as.character(),
                               coef = as.numeric(),
                               lower_level = as.numeric(), 
                               upper_level = as.numeric(),
                               type = as.character())
    
    # Create data for all variables in the model, adding onto the splined variables
    model_info <- tibble(var = names(model_vars[!grepl("knot|factor\\(", names(model_vars))][-1]),
                         coef = as.numeric(t(model_vars[!grepl("knot|factor\\(", names(model_vars))][-1])),
                         lower_level = -Inf, 
                         upper_level = Inf,
                         levels = NA,
                         type = "as_is") %>%
        # Clean up names
        mutate(var = str_replace_all(var, ".*(?<!,knot=c)\\(|\\)\\d", "")) %>%
        # Add splined and factor variables
        bind_rows(spline_info, factor_info)
    
    # Return data
    return(model_info)
}

# Create function to calculate sample linear predictor
lpsamp <- function(df){
    # Get model information
    model_info <- model_inf()
    
    # Calculate LP fractions
    lp_fracs <- bind_cols(lapply(1:nrow(model_info), \(x){
        # Get variable of interest
        var <- model_info[x, "var"][[1]]
        
        # Get coefficient of interest
        coef <- model_info[x, "coef"][[1]]
        
        # Get lower level of section
        lower_level <- model_info[x, "lower_level"][[1]]
        
        # Get upper level of section
        upper_level <- model_info[x, "upper_level"][[1]]
        
        # Get type of variable
        type <- model_info[x, "type"][[1]]
        
        # Get level of variable
        level <- model_info[x, "levels"][[1]]
        
        # Determine whether variable is logical
        if(length(table(df[[var]])) == 2 | type == "factor") logic <- TRUE else logic <- FALSE
        
        # Calculate fraction of LP
        lp_frac <- df %>%
            # Select relevant variables
            select(.imp, all_of(var)) %>%
            # Rename varying 'var' to value
            rename(value = 2) %>%
            # Arrange for grouping
            arrange(.imp) %>%
            # Group per imputation
            group_by(.imp) %>%
            # Calculate new variables
            mutate(# Set value based on conditions
                value = case_when(type == "as_is" ~ value,                                                # No changes for as_is variables
                                  type == "spline" & value >= lower_level & value < upper_level ~ value,  # No changes if value is in spline section
                                  type == "factor" & as.character(value) == as.character(level) ~ 1,      # Factor to 1 if same level (this represents dummy variable)
                                  type == "factor" & as.character(value) != as.character(level) ~ 0,      # Factor to 0 if not the same level
                                  .default = NA),                                                         # Otherwise, set value to missing
                # Get mean of value
                mean_value = mean(value, na.rm = TRUE),
                # Round mean value if logical is true
                mean_value = ifelse(logic, round(mean_value), mean_value),
                # Multiply mean value with coefficient
                lp = mean_value * coef) %>%
            # Keep one row per group
            slice(1L) %>%
            # Ungroup again
            ungroup() %>%
            # Select only relevant variables
            select(lp) %>%
            # Rename to make name unique
            set_colnames(paste0("frac_", x))
        
        # Return data
        return(lp_frac)
    })) 
    
    # Get and return linear predictor of sample
    return(mean(rowSums(lp_fracs, na.rm = TRUE)))
}

# Create function to get predicted risks
pred <- function(df, model, observed, time = NULL, lpsamp = NULL, aft_dist = NULL){
    # Get model information
    model_info <- model_inf()
    
    # Calculate LP fractions
    lp_fracs <- bind_cols(lapply(1:nrow(model_info), \(x){
        # Get variable of interest
        var <- model_info[x, "var"][[1]]
        
        # Get coefficient of interest
        coef <- model_info[x, "coef"][[1]]
        
        # Get lower level of section
        lower_level <- model_info[x, "lower_level"][[1]]
        
        # Get upper level of section
        upper_level <- model_info[x, "upper_level"][[1]]
        
        # Get type of variable
        type <- model_info[x, "type"][[1]]
        
        # Get level of variable
        level <- model_info[x, "levels"][[1]]
        
        # Calculate fraction of LP
        lp_frac <- df %>%
            # Select relevant variables
            select(studynr, .imp, all_of(var)) %>%
            # Rename varying 'var' to value
            rename(value = 3) %>%
            # Arrange for grouping
            arrange(studynr, .imp) %>%
            # Group per imputation
            group_by(studynr, .imp) %>%
            # Calculate fraction
            mutate(lp_frac = case_when(type == "as_is" ~ value * coef,                                                # As_is variables
                                       type == "spline" & value >= lower_level & value < upper_level ~ value * coef,  # Splines
                                       type == "factor" & as.character(value) == as.character(level) ~ coef,          # Factors
                                       .default = 0)) %>%  
            # Ungroup again
            ungroup() %>%
            # Keep only LP fraction
            select(lp_frac) %>%
            # Rename to make name unique
            set_colnames(paste0("frac_", x))
        
        # Return data
        return(lp_frac)
    })) 
    
    # Get observed event
    obs = arrange(df, studynr, .imp)[[deparse(substitute(observed))]]
    
    # Get time (for survival models)
    tim = arrange(df, studynr, .imp)[[deparse(substitute(time))]]
    
    # Get final linear predictors
    dat_tmp <- df %>%
        # Arrange on imputation
        arrange(studynr, .imp) %>%
        # Keep only imputation
        select(studynr, .imp) %>%
        # Add variables of interest
        mutate(# Linear predictor
            lp = rowSums(lp_fracs, na.rm = TRUE),
            # Observed event
            observed = obs,
            # Time (for survival models)
            time = tim)
    
    # For Cox/FIne-Gray models
    if(model %in% c("cox", "fine-gray")){
        ## Calculate individual risks
        # The predicted risks from predict.coxph() from {survival} do not correspond to our calculated probabilities, because they are
        # the survival at the observed follow-up time for that individual. Instead, predictCox() from {riskRegression} gives predictions
        # at a certain timepoint: predictCox(fit, newdata = df, times = 365, type = "survival")[["survival"]]
        # https://community.rstudio.com/t/predict-probability-of-default-in-cox-ph-model-and-get-baseline-hazard/118953 
        dat_tmp <- dat_tmp %>%
            # Calculate new variables
            mutate(# Centered linear predictors
                lp = lp - lpsamp,
                # Predicted probability
                prob = 1 - exp(-as.numeric(model_vars[["bh"]])) ^ exp(lp)) %>%
            # Ungroup again
            ungroup()
    }
    
    # For AFT models (currently only for lognormal and Weibull distributions)
    if(model == "aft"){
        # Calculate individual times
        dat_tmp <- dat_tmp %>%
            # Calculate expected time to failure
            mutate(prob = survreg.distributions[[aft_dist]][["itrans"]](as.numeric(model_vars[["a"]]) + lp))
    }
    
    # For logistic models
    if(model == "logistic"){
        # Calculate individual risks
        dat_tmp <- dat_tmp %>%
            # Calculate predicted probability
            mutate(prob = 1 / (1 + exp(-(as.numeric(model_vars[["a"]]) + lp))))
    }
    
    # Take mean of predicted value
    dat_tmp <- dat_tmp %>%
        # Sort for grouping
        arrange(studynr) %>%
        # Group on person
        group_by(studynr) %>%
        # Calculate final predicted probability
        mutate(pred = mean(prob)) %>%
        # Keep one row per person
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Drop irrelevant columns
        dplyr::select(-prob, -.imp)
    
    # Return data
    return(dat_tmp)
}
