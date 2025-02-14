# Performing calculations on the data

#' @title Outlier detection in triplicate data
#' 
#' @description This function identifies outlier wells in triplicate data.
#'              If 3 wells are present:
#'                - Find the two closest values.
#'                - Mark the third one as an outlier if it differs more than 1 cycle.
#'                - If all three are different, mark all three.
#'              If 2 wells are present:
#'                - If they differ by more than 1, highlight both.
#'              If 1 well is present:
#'                - Highlight it.
#' 
#' 
#' @param cq_data A vector of length 3 containing the Cq values for a given 
#'                sample-target pair. NA values should be included.
#' 
#' @return A logical vector indicating whether each well is an outlier. An NA
#'         value indicates that the well is missing and is always considered an
#'         outlier.
#' 
#' @export
#' 
detect_outlier_well <- function(cq_data) {
  if (length(cq_data) != 3) {
    stop("The input vector must have a length of 3.")
  }

  valid_values <- na.omit(cq_data)
  n <- length(valid_values)

  if (n <= 1) {
    return(rep(TRUE, 3))
  }

  if (n == 2) {
    if (abs(valid_values[1] - valid_values[2]) > 1) {
      return(rep(TRUE, 3))
    } else {
      # Only the NA value is considered TRUE
      return(is.na(cq_data))
    }
  }

  if (n == 3) {
    # Compute pairwise absolute differences
    diffs <- abs(outer(valid_values, valid_values, "-"))
    diag(diffs) <- NA  # Ignore self-comparisons
    closest_pair <- which(diffs == min(diffs, na.rm = TRUE), arr.ind = TRUE)[1, ]

    # Check if the closest pair differs by more than 1 cycle
    if (diffs[closest_pair[1], closest_pair[2]] > 1) {
      # If so, all three are outliers
      return(rep(TRUE, 3))
    } else {
      # The one left is the potential outlier
      potential_outlier <- setdiff(1:3, closest_pair)

      # Check if the potential outlier differs by more than 1 cycle from the closest pair
      if (diffs[closest_pair[1], potential_outlier] > 1 | diffs[closest_pair[2], potential_outlier] > 1) {
        # The one left is the outlier
        rslt <- rep(FALSE, 3)
        rslt[potential_outlier] <- TRUE
        return(rslt)
      } else {
        # All three are valid
        return(rep(FALSE, 3))
      }
    }
  }
}


#' @title Highlight outlier wells in triplicate data
#' 
#' @description This function highlights outlier wells in triplicate data.
#' 
#' @param data A validated data frame of all samples, targets, and wells.
#' 
#' @return A data frame with an additional column indicating whether each well
#'         is an outlier.
#' 
#' @seealso \code{\link{detect_outlier_well}}
#' 
#' @import dplyr
#' 
#' @export
#' 
highlight_outliers <- function(data) {
  # Group by Sample and Target and add an Outlier column
  # Add an outlier_sample column to the data frame and TRUE if any well in
  # the sample-target pair is an outlier
  data <- data %>%
    group_by(Sample, Target) %>%
    mutate(Outlier = detect_outlier_well(Cq)) %>%
    mutate(Outlier_no_na = Outlier & !is.na(Cq)) %>%
    mutate(Outlier_Sample = any(Outlier)) %>%
    mutate(Outlier_Sample_no_na = any(Outlier_no_na)) %>%
    ungroup()

  # Identify different categories of outliers
  na_wells <- data %>%
    filter(is.na(Cq)) %>%
    select(Sample, Target, Group, Well, Cq)
  
  outlier_wells <- data %>%
   group_by(Sample, Target) %>%
    filter(Outlier_Sample_no_na) %>%
    ungroup() %>%
    select(Sample, Target, Group, Well, Cq, Outlier)

  cat("The following wells must be excluded due to missing values:\n")
  print(na_wells)
  cat("\n")
  cat("The following wells are outliers due to Cq values:\n")
  print(outlier_wells)
  cat("\n")
  return(data)
}

#' @title Get invalid wells
#' 
#' @description This function returns list of two vectors, one containing the
#'              wells with missing values and the other containing the wells
#'              with outlier values.
#' 
#' @param data A validated data frame of all samples, targets, and wells. Must
#'             be processed by \code{\link{highlight_outliers}}.
#' 
#' @return A list of two vectors, one containing the wells with missing values
#'         and the other containing the wells with outlier values.
#' 
#' @export
#' 
#' @import dplyr
#' 
get_invalid_wells <- function(data) {
  na_wells <- data %>%
    filter(is.na(Cq)) %>%
    select(Sample, Target, Group, Well, Cq)
  
  outlier_wells <- data %>%
    filter(Outlier_no_na) %>%
    select(Sample, Target, Group, Well, Cq)
  
  return(list(na_wells = na_wells, outlier_wells = outlier_wells))
}


#' @title Exclude invalid wells
#' 
#' @description This function excludes invalid wells from the data.
#' 
#' @param data A validated data frame of all samples, targets, and wells. Must
#'             be processed by \code{\link{highlight_outliers}}.
#' 
#' @param exclude_wells A manually curated vector of wells to exclude. These
#'                      should be wells that are manually determined to be
#'                      invalid based on results from \code{\link{highlight_outliers}}
#'                      and \code{\link{get_invalid_wells}}. Only outlier wells
#'                      are necessary to include in this vector. NA wells will
#'                      be automatically excluded.
#' 
#' @return A data frame with invalid wells marked to be excluded.
#' 
#' @export
#' 
#' @import dplyr
#' 
exclude_invalid_wells <- function(data, exclude_wells) {
  # Find and alert excluded wells that are valid
  invalid_exclusion <- setdiff(exclude_wells, data$Well)
  if (length(invalid_exclusion) > 0) {
    warning(paste("The following wells are not present in the data and will be ignored:", paste(invalid_exclusion, collapse = ", ")))
  }

  data <- data %>%
    select(Well, Target, Group, Sample, Cq) %>%
    mutate(Excluded = is.na(Cq) | Well %in% exclude_wells)

  return(data)
}

#' @title Calculate the mean Cq value for each sample-target pair
#' 
#' @description This function calculates the mean Cq value for each sample-target
#'              pair after excluding invalid wells.
#' 
#' @param data A validated data frame of all samples, targets, and wells. Must
#'            be processed by \code{\link{highlight_outliers}} and
#'            \code{\link{exclude_invalid_wells}}.
#' 
#' @return A data frame with the mean Cq value for each sample-target pair.
#' 
#' @export
#' 
#' @import dplyr
#' 
calculate_mean_cq <- function(data) {
  # Calculate the mean Cq value for each sample-target pair
  data <- data %>%
    filter(!Excluded) %>%
    mutate(
        Sample = factor(Sample, levels = unique(Sample)),
        Target = factor(Target, levels = unique(Target))
    )
  if (length(unique(data$Group)) > 1) {
    data <- data %>%
      group_by(Target, Sample, Group) %>%
      summarise(Mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
      ungroup() %>%
      mutate(Sample = as.character(Sample), Target = as.character(Target)) %>%
      select(Target, Group, Sample, Mean_Cq)
  } else {
    data <- data %>%
      group_by(Target, Sample) %>%
      summarise(Mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
      ungroup() %>%
      mutate(Sample = as.character(Sample), Target = as.character(Target)) %>%
      select(Target, Sample, Mean_Cq)
  }
}

#' @title Calculate relative expression
#' 
#' @description This function calculates the relative expression of each sample
#'              compared to a reference sample for each target.
#' 
#' @param data A validated data frame of all samples, targets, and wells. Must
#'             be processed by \code{\link{highlight_outliers}} and
#'             \code{\link{exclude_invalid_wells}}.
#' 
#' @param refs One or more targets that serve as the housekeeping gene(s)
#' 
#' @return A list containing the raw data and the relative expression of each
#'         sample compared to the reference sample for each target.
#' 
#' @export
#' 
#' @import dplyr rlang
#' 
calculate_relative_expression <- function(data, refs) {
  expr_results <- list(raw_data = data)

  # Calculate the mean Cq value for each sample-target pair
  mean_cq <- calculate_mean_cq(data)

  # Check if any na values are present
  if (any(is.na(mean_cq$Mean_Cq))) {
    stop("NA values are present in the mean Cq data.")
  }

  # Check if any reference genes are missing
  if (any(!refs %in% mean_cq$Target)) {
    stop("One or more reference genes are missing from the data.")
  }

  if (length(refs) == 0) {
    stop("At least one reference gene must be provided.")
  }

  # Calculate the relative expression of each sample compared to the reference sample
  for (ref in refs) {
    ref_mean_cq <- mean_cq %>%
      filter(Target == ref) %>%
      select(Sample, Mean_Cq) %>%
      rename(!!paste0("Mean_Cq_", ref) := Mean_Cq)

    result <- mean_cq %>%
      filter(Target != ref) %>%
      left_join(ref_mean_cq, by = "Sample") %>%
      mutate(delta_Cq = Mean_Cq - !!sym(paste0("Mean_Cq_", ref))) %>%
      mutate(Relative_Expr = 2 ^ -delta_Cq) %>%
      mutate(Reference = ref) %>%
      select(Target, Group, Sample, Mean_Cq, Reference, !!sym(paste0("Mean_Cq_", ref)), delta_Cq, Relative_Expr)

    expr_results[[ref]] <- result
  }

  return(expr_results)
}

