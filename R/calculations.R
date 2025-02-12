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
    select(Sample, Target, Well, Cq)
  
  outlier_wells <- data %>%
   group_by(Sample, Target) %>%
    filter(Outlier_Sample_no_na) %>%
    ungroup() %>%
    select(Sample, Target, Well, Cq, Outlier)

  cat("The following wells must be excluded due to missing values:\n")
  print(na_wells)
  cat("\n")
  cat("The following wells are outliers due to Cq values:\n")
  print(outlier_wells)
  cat("\n")
  return(data)
}