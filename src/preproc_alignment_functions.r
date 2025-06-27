
# Function to handle missing values with time dependency for both continuous and binary variables
impute_missing_time_dependent.3 <- function(ts, variable_types, num_imputations = 10) {
    
    time <- ts$time
    id <- ts$id
    features <- ts
    features_df <- as.data.frame(features)
    
    imputation_methods <- c("","",sapply(variable_types, function(type) {
        if (type == "binary") {
            return("logreg")  # Logistic regression for binary variables
        } else {
            return("rf")
        }
    }))
    
    print("imputation_methods");print(imputation_methods)
    print("features_df)");print(head(features_df))
    print("apply(features_df,2,function(x) sum(!is.na(x)))");print(apply(features_df,2,function(x) sum(!is.na(x))))
    print("dim(features_df)");print(dim(features_df))
    print("dim(unique(features_df))");print(dim(unique(features_df)))
    
    ##se imputa todo revuelto
    imputed <- suppressWarnings(mice(features_df, m = num_imputations,
                                     method = imputation_methods, maxit = 3,
                                     seed = 123, group = "id"))
    
    print("imput'o")
    
                                        ## Combine results (locally weighted regression for continuous, majority vote for binary)
    completed_imputations <- lapply(1:num_imputations, function(i) {
        complete(imputed, action = i)
    })
    
    print("extrajo")
    print("class(completed_imputations)");print(class(completed_imputations))
    
    ids =  unique(completed_imputations[[1]]$id)
    
    features_df.init =  features_df
    completed_imputations.init =  completed_imputations
    
    print("class completed_imputations.init");print(class(completed_imputations.init[[1]]))
    print("head completed_imputations.init");print(head(completed_imputations.init[[1]]))
    print("a cada id sacarle su imputacion combinada")
    
    features_df_list = lapply(ids, function(id.pp){
        features_df = subset(features_df.init, id == id.pp)
        if( dim(features_df)[1] <=2) print(features_df) 
        completed_imputations = lapply(completed_imputations.init, function(m) subset(m, id == id.pp ) )
        features. = ts[, -grep("time|id", names(features_df))]
        
        for (j in 1:ncol(features.)) {
            if (variable_types[j] == "binary") {
                ## Majority vote for binary variables
                binary_columns <- sapply(completed_imputations, function(df) df[, j])
                features_df[, j] <- apply(binary_columns, 1, function(values) {
                    if (all(is.na(values))) {
                        return(NA)
                    }
                    return(as.numeric(names(sort(table(values), decreasing = TRUE)[1])))
                })
            } else {
                ## Locally weighted regression for continuous variables
                
                imputed_columns <- lapply(completed_imputations, function(df) {
                    
                    df = df[,-1]
                    
                    na_indices <- is.na(features_df[, j+2])
                    if (any(na_indices)) {
                        
                        loess_fit <- loess(df[, j+1] ~ as.numeric(df$time) , na.action = na.omit)
                        df[na_indices, j+1] <- predict(loess_fit, newdata = as.numeric(df$time)[na_indices])
                    }
                    return(df[, j+1])
                })
                features_df[, j+2] <- Reduce("+", imputed_columns) / num_imputations
            }
        }
        features_df
    })
    print("Recombine time with imputed features")
     imputed_ts <- do.call(rbind,features_df_list)
    print("completed imputation!")
    return(imputed_ts)
}



process_ts_list = function(ts_list, num_imputations = 5,                           
                           variable_types =  c("continuos","continuos", "binary", "continuos","continuos")
                           ) {
    ts_list = as.data.frame(ts_list)
    imputed_list <- lapply(ts_list, function(ts) {
        impute_missing_time_dependent(ts, num_imputations, variable_types = variable_types)
  })
    return(imputed_list)
}



# Function to interpolate time series to uniform time steps, handling continuous and binary variables
interpolate_to_uniform_time <- function(ts_list, variable_types, step_size = 1/2) {
  interpolated_list <- lapply(ts_list, function(ts) {
    # Extract time and features
      original_time <- ts[, "time"]
      id <- ts[1, "id"]
    features <- ts[, -grep("time|id", colnames(ts))]

    # Generate uniform time steps that are multiples of step_size
    start_time <- round(min(original_time) / step_size) * step_size
    end_time <- round(max(original_time) / step_size) * step_size
    uniform_time <- seq(from = start_time, to = end_time, by = step_size)
      
    # Interpolate each feature based on type
    interpolated_features <- sapply(1:ncol(features), function(i) {
      if (variable_types[i] == "binary") {
        # For binary variables, use nearest neighbor interpolation
          return(approx(original_time, features[, i], xout = uniform_time, method = "constant", rule = 2)$y)
      } else {
        # For continuous variables, use linear interpolation
          return(approx(original_time, features[, i], xout = uniform_time, rule = 2)$y)
      }
    })

    # Combine uniform time and interpolated features
      interpolated_ts <- as.data.frame(cbind(uniform_time, interpolated_features))
      names(interpolated_ts) = colnames(ts)[-grep("id",colnames(ts))]
      interpolated_ts$id = id
      return(interpolated_ts)
      
  })
  return(interpolated_list)
}


## Consensus Alignment Function for Multivariate Time Series
consensus_alignment <- function(series_list, max_iter = 10, tol = 1e-6) {
    ## Step 1: Pad all sequences to the same length
    ## max_length <- max(sapply(series_list, nrow))
    
    ii.longest = which.max(unlist(lapply(series_list, nrow)))
    
    max.time = max(unlist(lapply(series_list, function(x) x[,"time"]) ))
    min.time = min(unlist(lapply(series_list, function(x) x[,"time"] )))
    
    print("min.time");print(min.time)
    print("max.time");print(max.time)
    
    ref.init = series_list[[ii.longest]]        
    ref.init = ref.init
    ref.init = extend_time_series(data = ref.init, min_time = min.time, max_time = max.time, time_step = ref.init[2,"time"] - ref.init[1,"time"])
    
    series_list_carried  = lapply(series_list, function(ts) extend_time_series(data = ts, min_time = min.time, max_time = max.time, time_step = ts[2,"time"] - ts[1,"time"]))
    
    ref = ref.init    
    
    ## Iterative alignment and refinement
    for (iter in 1:max_iter) {
        ##print("iter: ");print(iter)
        ## Step 3: Align each sequence to the reference
        aligned_list <- lapply(series_list_carried, function(ts) {
            
            ref = as.matrix(ref)
            ts = as.matrix(ts)
            
            alignment <- dtw(ref[,-1], ts[,-1], keep = T,
                             step.pattern = rabinerJuangStepPattern(4, "d", smoothed = TRUE),
                             dist.method = "Euclidean"
                             )  
            ## print("paso align step")
            
            ## Warp sequence based on alignment path
            aligned <- matrix(NA, nrow = dim(ref)[1], ncol = ncol(ts))  # Initialize aligned matrix
            aligned[alignment$index1, ] <- ts[alignment$index2, ]  # Map alignment

            aligned[,1] =  ref[,"time"]
            
            ##print("aligned")  ;print(aligned)
            ##print("ts")  ;
            ##print(class(aligned));print(dim(aligned))
            
            colnames(aligned) = colnames(ts)
            
            ##print("aligned")  ;print(aligned)
          
          return(aligned)
          
    })
        
    ## Step 4: Update consensus reference
        new_ref <- NA*ref
        ##print("new_ref");print(new_ref)
        
    for (i in 1:nrow(ref)) {
        for (j in 1:(ncol(ref))) {
   ##        print(c(i,j))
        vals <- sapply(aligned_list, function(x) x[i, j])  # Collect column values
        valid_vals <- vals[!is.na(vals)]  # Filter valid values
        if (length(valid_vals) > 0) {
          new_ref[i, j] <- median(valid_vals)  # Compute mean
        }
      }
    }
      
    ## Step 5: Check for convergence
    diff <- sum((new_ref - ref)^2, na.rm = TRUE)
    if (diff < tol) break
    ref <- new_ref
  }

##print("class(aligned_list)");print((aligned_list))
  
aligned_list = lapply(aligned_list, function(x) {x = as.data.frame(x);x$time = ref$time; x } )
    
  # Return final consensus reference and aligned sequences
  return(list(reference = ref, aligned = aligned_list, ref.init = ref.init))
}

extend_time_series <- function(data, min_time, max_time, time_step) {
    ## Create full time sequence
    print("data");print(head(data,5))
    full_time <- seq(from = min_time, to = max_time, by = time_step)
  
    ## Merge with the full time range
    extended_data <- merge(data.frame(time = full_time), data, by = "time", all.x = TRUE)
    
    ## Carry backward (fill NA at the beginning)
    extended_data[is.na(extended_data)] <- NA
    extended_data <- apply(extended_data, 2, function(x) {
        x[1:which(!is.na(x))[1] - 1] <- x[which(!is.na(x))[1]]
        x
    })
    
    ## Carry forward (fill NA at the end)
  extended_data <- apply(extended_data, 2, function(x) {
    x[which(!is.na(x))[length(which(!is.na(x)))]:length(x)] <- x[which(!is.na(x))[length(which(!is.na(x)))]]
    x
  })

  ## Return result as a data frame
  extended_data <- as.data.frame(extended_data)
  return(extended_data)
}



plot_feature_all.2 =
function(ts_list, colapsar = F, ref.df)  {
    
  combined_df <- do.call(rbind, lapply(seq_along(ts_list), function(i) {
    ts <- as.data.frame(ts_list[[i]])
    ##    data.frame(time = ts[, 1], ts[,-1], id = as.factor(i))
        data.frame(time = ts$time, ts[,-grep("time",names(ts))], id = as.factor(i))
  }))

    combined_df =  melt(combined_df, id.vars = c("time","id"))


    if(colapsar == T){p <- ggplot(combined_df, aes(x = time, y = value, color = id)) +
                             geom_line(size = 0.5) +##geom_point(size = 1) +
                             labs(title = "",
                                  x = "Time", y = "Value", color = "ID") +
                             theme_minimal() ##+ facet_grid(variable~., scale = "free")}else{
    
                         p <- ggplot(combined_df, aes(x = time, y = value, color = id)) +
                             geom_line(size = 0.5) +##geom_point(size = 1) +
                             labs(title = "",
                                  x = "Time", y = "Value", color = "ID") +
                             theme_minimal() ##+ facet_grid(id~variable, scale = "free") 
                     }
 
    return(

        p + geom_line(data = ref.df, aes(x =  time, y = value), col = "black")+ scale_y_log10() + facet_wrap(~variable, ncol = 3, scale = "free")  + theme(legend.position = "none")


           )
     
}
 

## Function to compute bootstrap-based significance of DTW alignment contribution per time point
compute_bootstrap_confidence_dtw <- function(aligned_list, consensus, n_boot = 1000) {
    
    ## Function to compute DTW alignment cost contribution per time point
    compute_dtw_contributions <- function(ts) {
        dtw_res <- dtw(consensus, ts, step.pattern = rabinerJuangStepPattern(4, "d", smoothed = TRUE), keep = TRUE, method = "Euclidean")
        
        ## Extract warping path
        path_consensus <- dtw_res$index1  ## Time indices in consensus
        path_individual <- dtw_res$index2 ## Time indices in individual series
    
        ## Extract cumulative cost matrix
        cost_matrix <- dtw_res$costMatrix
        total_cost <- dtw_res$distance  # Total DTW cost
    
        ## Compute stepwise cost contributions
        step_costs <- numeric(nrow(consensus))  # Initialize vector for contributions
        stepwise_differences <- diff(cost_matrix[cbind(path_consensus, path_individual)], differences = 1)
        
        ## Assign cost differences to corresponding time points
        step_costs[path_consensus[-1]] <- stepwise_differences
        ## Normalize by total DTW cost
        contribution <- abs(step_costs) / total_cost
    
    return(contribution)
  }
  
    ## Stack aligned sequences into a 3D array: (time × variables × individuals)
    aligned_array <- array(unlist(aligned_list), dim = c(nrow(consensus), ncol(consensus), length(aligned_list)))
    
    ## Compute observed DTW contributions per time point
  observed_contributions <- apply(aligned_array, 3, compute_dtw_contributions)
    
    x11(title = "e");par(mfrow = c(1,2));matplot(log2(observed_contributions), cex = 0.5)
    print("dim(observed_contributions)");print(dim(observed_contributions))
    
    costs.per.sequence =  apply(observed_contributions,2, sum)
    qq.costs.per.sequence =  quantile(costs.per.sequence)
    print("min(observed_contributions)");print(min(observed_contributions))
    print(qq.costs.per.sequence)
    ii = which(costs.per.sequence<=qq.costs.per.sequence["25%"])    
    print("ii");print(ii)
    observed_contributions = observed_contributions[,ii]
    aligned_array = aligned_array[,,ii]
    x11("cost contibutions before the bt");matplot(log2(observed_contributions), cex = 0.5)
    mean_observed_contribution <- rowMeans(observed_contributions, na.rm = TRUE)
    
    ## Bootstrap resampling
    boot_samples <- replicate(n_boot, {
        boot_indices <- sample(seq_len(dim(aligned_array)[3]), replace = TRUE)
        boot_sample <- aligned_array[, , boot_indices, drop = FALSE]
        rowMeans(apply(boot_sample, 3, compute_dtw_contributions), na.rm = TRUE)
  })
    
    ## Compute 95% confidence intervals
    ci_lower <- apply(boot_samples, 1, quantile, probs = 0.025, na.rm = TRUE)
    ci_upper <- apply(boot_samples, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    ## Compute p-values for significance testing
    p_values <- rowMeans(boot_samples >= mean_observed_contribution)  # Right-tailed test: high contribution = significant

  return(list(
    observed_contributions = mean_observed_contribution,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_values = p_values
  ))
}


compute_dtw_discrepancy_matrix.level <-  function(discrepancy_matrix, similar= FALSE, threshold = "50%"){
    ##matrix :  times x vars
    discrepancy_matrix.normalized = abs(discrepancy_matrix[,-grep("time", names(discrepancy_matrix))])
    print("head(discrepancy_matrix.normalized)");print(head(discrepancy_matrix.normalized))
    print("tail(discrepancy_matrix.normalized)");print(tail(discrepancy_matrix.normalized))
    qq_discrepancy_matrix.normalized = quantile(as.matrix(discrepancy_matrix.normalized), seq(0,1,by = 0.05), na.rm = T)
    print("qq_discrepancy_matrix.normalized");print(qq_discrepancy_matrix.normalized)
  
    if(similar == TRUE){
        ## the best aggreements
        row.pick =  qq_discrepancy_matrix.normalized[which(names(qq_discrepancy_matrix.normalized) == threshold)]
        print("row.pick");print(row.pick)
        accepted_level = discrepancy_matrix.normalized<=row.pick      
    }else{
        ## the largest discrepancies
        row.pick =  qq_discrepancy_matrix.normalized[which(names(qq_discrepancy_matrix.normalized) == threshold )]##    
        print("row.pick");print(row.pick)
        accepted_level = discrepancy_matrix.normalized>=row.pick
    }
    indices = accepted_level
    discrepancy_matrix.normalized*ifelse(indices,indices +0, NA)
}


compute_dtw_discrepancy_matrix <- function(T_ref, T_list) {
    ## Ensure all time series have the same dimensions
    num_timepoints <- nrow(T_ref)
    num_variables <- ncol(T_ref)
    num_series <- length(T_list)
    
    ## Initialize discrepancy matrix
    dtw_discrepancy_matrix <- matrix(0, nrow = num_timepoints, ncol = num_variables)
    
    ## Compute DTW discrepancy for each variable separately
    discrepancy_values <- map(T_list, function(T_i) {
        map(1:num_variables, function(v) {
            dtw_res <- dtw(T_ref[, v], T_i[, v], step.pattern = rabinerJuangStepPattern(6, "c"))
      
            ## Extract warping path indices
            ref_indices <- dtw_res$index1
            seq_indices <- dtw_res$index2
            
            T_ref[T_ref == 1 ] = 2  
            T_i[T_i == 1] = 2
            
            T_ref[T_ref == 0 ] = 1  
            T_i[T_i == 0] = 1
            
            abs_diff <- (T_ref[ref_indices, v] - T_i[seq_indices, v])/(T_ref[ref_indices, v])
            
            discrepancy_vector <- numeric(num_timepoints)
            
            ## Accumulate discrepancies at reference indices
            for (i in seq_along(ref_indices)) {
                discrepancy_vector[ref_indices[i]] <- discrepancy_vector[ref_indices[i]] + abs_diff[i]
            }
            
      ## Normalize by occurrence count
      index_counts <- table(ref_indices)
            discrepancy_vector[as.numeric(names(index_counts))] <- 
                discrepancy_vector[as.numeric(names(index_counts))] / as.numeric(index_counts)
            
            discrepancy_vector
        }) %>% reduce(cbind)  # Combine variable-wise results into a matrix
    })
    
    
    print("discrepancy_values"); print(class(discrepancy_values)); print(length(discrepancy_values))
    
    ## Aggregate discrepancy across all sequences
    dtw_discrepancy_matrix <- reduce(discrepancy_values, `+`) / num_series
    
    return(dtw_discrepancy_matrix)
}


compute_median_variance <- function(group1, group2, relative = TRUE) {
    ## print( "Ensure lists are not empty")
  if (length(group1) == 0 || length(group2) == 0) {
      stop("Both groups must contain at least one individual.")
  }
  
 ##  print("Extract matrix dimensions (assuming all matrices have the same dimensions")
  time_points <- nrow(group1[[1]])
  variables <- ncol(group1[[1]])
  
    ##print("Function to compute median difference and variance at each (t, v)")
  compute_stats <- function(t, v) {
    values1 <- map_dbl(group1, ~ .x[t, v])  # Extract time-variable values from group1
    values2 <- map_dbl(group2, ~ .x[t, v])  # Extract time-variable values from group2

    ##print("Remove NAs")
    values1 <- values1[!is.na(values1)]
    values2 <- values2[!is.na(values2)]
    
    if (length(values1) == 0 || length(values2) == 0) {
      warning(paste("No valid comparison at time", t, "variable", v, "(NaN returned)"))
      return(c(NaN, NaN))
    }
    
    if(relative == TRUE) {median_diff <- (median(values1, na.rm = T) - median(values2, na.rm = T)) / median(values1, na.rm = T)}else{median_diff <- (median(values1, na.rm = T) - median(values2, na.rm = T))}
        
    diffs <- outer(values1, values2, FUN = "-")  # Pairwise differences
    var_diff <- var(as.vector(diffs), na.rm = TRUE)
    
    sign_diff =  sum(as.vector(sign(diffs)), na.rm = TRUE)/length(diffs)

    out = list(median_diff = median_diff, var_diff = var_diff, sign_diff =  sign_diff, pair.comparison.availability = length(diffs)  )
    
    return(out)
    
  }
  
  results <- lapply(1:variables, function(v) lapply(1:time_points, function(t) {
      compute_stats(t,v)} ) )

median_diffs = do.call(cbind,lapply(1:length(results),function(ii) unlist(lapply(lapply(lapply(results[[ii]], function(x) lapply(x, function(y) lapply(y, function(z) z))), function(w) w), function(z) z$median_diff))))

    var_diffs = do.call(cbind,lapply(1:length(results),function(ii) unlist(lapply(lapply(lapply(results[[ii]], function(x) lapply(x, function(y) lapply(y, function(z) z))), function(w) w), function(z) z$var_diff))))

    sign_diff = do.call(cbind,lapply(1:length(results),function(ii) unlist(lapply(lapply(lapply(results[[ii]], function(x) lapply(x, function(y) lapply(y, function(z) z))), function(w) w), function(z) z$sign_diff)))) 

    pair.comparison.availability = do.call(cbind,lapply(1:length(results),function(ii) unlist(lapply(lapply(lapply(results[[ii]], function(x) lapply(x, function(y) lapply(y, function(z) z))), function(w) w), function(z) z$pair.comparison.availability)))) 
    
    list(median_diff = as.data.frame(median_diffs), variance = as.data.frame(var_diffs), sign_diff = as.data.frame(sign_diff), pair.comparison.availability = as.data.frame(pair.comparison.availability))

}

remove_flat_ends <- function(data, threshold = 3, not.these = c("aghbs", "antihbc", "antihbs", "antihcv")) {
  ## Identify relevant columns (excluding time and not.these)
  value_cols <- setdiff(names(data), c("time", not.these))
  
    ## Function to find positions of flat starts/ends (returns logical vector)
  flat_mask <- function(x) {
    mask <- rep(FALSE, length(x))
    r <- rle(x)
    # Flat start
    if (r$lengths[1] >= threshold) {
      mask[1:r$lengths[1]] <- TRUE
    }
    ## Flat end
    if (r$lengths[length(r$lengths)] >= threshold) {
        mask[(length(x) - r$lengths[length(r$lengths)] + 1):length(x)] <- TRUE
    }
    return(mask)
  }
    
    ## Compute flat masks for each column
  masks <- lapply(data[, value_cols], flat_mask)

  ## Count number of TRUEs (i.e., flat ends) in each mask
    flat_counts <- sapply(masks, sum)

  ## Identify the column with the fewest flat positions
  ref_col <- names(which.min(flat_counts))
    ref_mask <- masks[[ref_col]]

    ## Apply NA based on reference mask
    for (col in value_cols) {
        data[ref_mask, col] <- NA
    }
    
    ## For not.these columns, set to NA where all value_cols are NA
    row_all_na <- apply(data[, value_cols], 1, function(x) all(is.na(x)))
    data[row_all_na, not.these] <- NA
    
    return(data)
}


# Custom function to compute DTW distance while treating NAs as gaps
compute_dtw_dist <- function(ts_list) {
  n <- length(ts_list)
  dist_matrix <- matrix(NA, n, n)
  indices <- combn(n, 2, simplify = FALSE)  
  dtw_distances <- map_dbl(indices, ~ {
    i <- .x[1]; j <- .x[2]
    ## Extract time series
    ts_i <- ts_list[[i]]
    ts_j <- ts_list[[j]]
    ## Ensure they are matrices
    if (!is.matrix(ts_i)) ts_i <- as.matrix(ts_i)
    if (!is.matrix(ts_j)) ts_j <- as.matrix(ts_j)
    ## Identify valid (non-NA) indices for each variable
    valid_idx <- complete.cases(ts_i) & complete.cases(ts_j)
    ##print("valid_idx");print(valid_idx)
    ## If no common valid points, return a large distance
    if (sum(valid_idx) == 0) {
      return(Inf)
    }
    ##Compute DTW on valid data points only
    dtw_res <- dtw(ts_i[valid_idx, , drop = FALSE], ts_j[valid_idx, , drop = FALSE],
    keep = T,step.pattern = rabinerJuangStepPattern(4, "d", smoothed = TRUE), dist.method = "Euclidean")
    ##print("dtw_res$distance");print(dtw_res$distance)  
    dtw_res$distance
  })
  ##print("class(dtw_distances)");print(class(dtw_distances))
  dist_matrix[lower.tri(dist_matrix)] <- dtw_distances
  dist_matrix[upper.tri(dist_matrix)] <- t(dist_matrix)[upper.tri(dist_matrix)]  # Fill upper triangle without propagating NAs
 ## print("dist_matrix");print(dist_matrix)
  diag(dist_matrix) <- 0  # Zero diagonal
  as.dist(dist_matrix)  # Convert to dist object
}


compute_dtw_dist.alternative <- function(ts_list) {
  n <- length(ts_list)
  dist_matrix <- matrix(NA, n, n)
  indices <- combn(n, 2, simplify = FALSE)  
  dtw_distances <- map_dbl(indices, ~ {
    i <- .x[1]; j <- .x[2]
    ## Extract time series
    ts_i <- ts_list[[i]]
    ts_j <- ts_list[[j]]
    ## Ensure they are matrices
    if (!is.matrix(ts_i)) ts_i <- as.matrix(ts_i)
    if (!is.matrix(ts_j)) ts_j <- as.matrix(ts_j)
    ## Identify valid (non-NA) indices for each variable
    valid_idx <- complete.cases(ts_i) & complete.cases(ts_j)
    ##print("valid_idx");print(valid_idx)
    ## If no common valid points, return a large distance
    if (sum(valid_idx) == 0) {
      return(Inf)
    }
    
dtw_res <- dtw(ts_i[valid_idx, , drop = FALSE], ts_j[valid_idx, , drop = FALSE],
  step.pattern = symmetric1,  ## use symmetric1 or 'rigid' for identity-like match
  window.type = "none",       ## no window constraint
  distance.only = TRUE,
  keep = FALSE,
  open.begin = FALSE,
  open.end = FALSE
)
    dtw_res$distance
  })

  dist_matrix[lower.tri(dist_matrix)] <- dtw_distances
  dist_matrix[upper.tri(dist_matrix)] <- t(dist_matrix)[upper.tri(dist_matrix)]  # Fill upper triangle without propagating NAs
  diag(dist_matrix) <- 0  # Zero diagonal
  as.dist(dist_matrix)  # Convert to dist object
}


##Alignment Function for Multivariate Time Series
alignment_to_reference <- function(series_list, reference.ts, max_iter = 10, tol = 1e-6) {
    ## Step 1: Pad all sequences to the same length
    ## max_length <- max(sapply(series_list, nrow))
    
    max.time = max(unlist(lapply(series_list, function(x) x[,"time"]) ))
    min.time = min(unlist(lapply(series_list, function(x) x[,"time"] )))

    max.time = max(reference.ts$time,max.time)
    min.time = min(reference.ts$time,min.time)
    
    print("min.time");print(min.time)
    print("max.time");print(max.time)
    
    ref.init = reference.ts
    
    ##  print("head(ref.init)");print(head(ref.init))
    ##  print("series_list[[1]]");print(head(series_list[[8]]))
    
    series_list_carried  = lapply(series_list, function(ts) extend_time_series(data = ts, min_time = min.time, max_time = max.time, time_step = ts[2,"time"] - ts[1,"time"]))

    message("carried the series")
    
    ref = ref.init    
    
    alignment.list = lapply( 1:length(series_list_carried), function(ii){
        ts = series_list_carried[[ii]]
        alignment <- dtw(ref[,-1], ts[,-1], keep = T,step.pattern = rabinerJuangStepPattern(4, "d", smoothed = TRUE), dist.method = "Euclidean")
        
        ##print("head(ref)");print(head(ref))
        ##print("head(ts)");print(head(ts))
        ##message("alineó")
        aligned <- as.data.frame(matrix(NA, nrow = nrow(ref), ncol = ncol(ref)))  # Initialize aligned matrix

        aligned[alignment$index1, ] <- ts[alignment$index2, ]  # Map alignment
        
        colnames(aligned)  = colnames(alignment.al.allpat.list[[ii]])
        
        aligned$time = consensus.reference.df.hbv.hcc.short$time
        
        aligned$id =  names(alignment.al.allpat.list)[[ii]]
        
        aligned
        
    }
)

    }
 
