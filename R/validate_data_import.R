# Validate data import

#' @title Validate data import
#' 
#' @param data Data to validate, a data frame or a path to a data file. The file
#'             must be a standard file for qPCR data. Each row must be a well.
#'             Columns should include the minimum required columns for qPCR data,
#'             including: Well, Target, Sample, Cq, and Group.
#' 
#' @return A data frame with the required columns for qPCR data. The column names
#'         are standardized to Well, Target, Sample, and Cq.
#' 
#' @importFrom readxl read_excel
#' @importFrom readr read_csv read_tsv
#' @importFrom tools file_ext
#' 
#' @export
#' 
import_and_validate <- function(data) {
  
  # If data is a path, read the data
  if (is.character(data)) {
    # Check extension to determine the file type
    extension <- tools::file_ext(data)
    if (extension %in% c("xlsx", "xls")) {
      data <- readxl::read_excel(data)
    } else if (extension == "csv") {
      data <- readr::read_csv(data)
    } else if (extension %in% c("txt", "tsv")) {
      data <- readr::read_tsv(data)
    } else {
      stop("File type not supported. Please provide a .xlsx, .xls, .csv, .txt, or .tsv file.")
    }
  }

  # Check if the data has the required columns and identify exact column names
  # for Well, Target, Sample, Cq, and Group
  well_col <- intersect(WELL, colnames(data))
  target_col <- intersect(TARGET, colnames(data))
  sample_col <- intersect(SAMPLE, colnames(data))
  cq_col <- intersect(CQ, colnames(data))
  group_col <- intersect(GROUP, colnames(data))

  # Check if the data has the required columns
  if (length(well_col) == 0) {
    stop("Data does not contain a column for Well.")
  } else if (length(well_col) > 1) {
    stop("Data contains multiple columns for Well. Please provide only one column. Available column names for well location: ", paste(WELL, collapse = ", "))
  }

  if (length(target_col) == 0) {
    stop("Data does not contain a column for Target.")
  } else if (length(target_col) > 1) {
    stop("Data contains multiple columns for Target. Please provide only one column. Available column names for target gene: ", paste(TARGET, collapse = ", "))
  }

  if (length(sample_col) == 0) {
    stop("Data does not contain a column for Sample.")
  } else if (length(sample_col) > 1) {
    stop("Data contains multiple columns for Sample. Please provide only one column. Available column names for sample: ", paste(SAMPLE, collapse = ", "))
  }

  if (length(cq_col) == 0) {
    stop("Data does not contain a column for Cq.")
  } else if (length(cq_col) > 1) {
    stop("Data contains multiple columns for Cq. Please provide only one column. Available column names for Cq value: ", paste(CQ, collapse = ", "))
  }

  if (length(group_col) == 0) {
    warning("Data does not contain a column for Group. Group information will not be available.")
    data$Group <- NA
  } else if (length(group_col) > 1) {
    stop("Data contains multiple columns for Group. Please provide only one column. Available column names for group: ", paste(GROUP, collapse = ", "))
  } else {
    colnames(data)[which(colnames(data) %in% group_col)] <- "Group"
    # When only one unique group is present, warn the user
    if (length(unique(data$Group)) == 1) {
      warning("Only one unique group is present in the data. Group information will not be used.")
    }
  }

  # Rename columns to standard names
  colnames(data)[which(colnames(data) %in% well_col)] <- "Well"
  colnames(data)[which(colnames(data) %in% target_col)] <- "Target"
  colnames(data)[which(colnames(data) %in% sample_col)] <- "Sample"
  colnames(data)[which(colnames(data) %in% cq_col)] <- "Cq"
  

  # Ensure Cq values are numeric and handle any NA values to NA
  data$Cq <- as.numeric(data$Cq)

  # Ensure that all sample-target pairs have exactly rows
  sample_target_pairs <- unique(data[, c("Sample", "Target")])
  for (i in 1:nrow(sample_target_pairs)) {
    pair <- sample_target_pairs[i, ]
    n <- sum(data$Sample == pair$Sample & data$Target == pair$Target)
    if (n != 3) {
      stop("Sample-Target pair ", pair$Sample, "-", pair$Target, " does not have exactly 3 rows.")
    }
  }
  
  return(data)
}