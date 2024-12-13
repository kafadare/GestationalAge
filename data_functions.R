## File of useful functions to load at top of script

find_diff_cols <- function(X) {
  replicate(nrow(X), names(
    Filter(function(z) length(z) > 1,
           lapply(X, unique))
  ), simplify = FALSE)
}

#Function to create custom size & step batches from df, sorted via given column
create_batch <- function(data, batch_size, method = split, batch_step = NULL){
  methods <- c("split", "sliding_window")
  method <- match.arg(method, methods)
  # Load or install necessary libraries
  required_pkg <- c("dplyr", "purrr")
  load_install_pkg(required_pkg)
  batches <- list()
  # Split the dataframe into regular batches
  if (method == "split"){
    batches <- df %>%
      group_by(batch = cumsum(row_number() %% batch_size == 1)) %>%
      ungroup() %>%
      nest()
  }
  #Split into sliding window batches
  else if (method == "sliding_window"){
    #check if batch step argument is given
    if (!is.null(batch_step)) {
      for (i in seq(1, ncol(data), batch_step)) {
        end <- min(i + batch_size - 1, ncol(data))
        batch <- data[, i:end]
        if(ncol(batch) == batch_size | (batch_size - ncol(batch)) < batch_step){
          batches[[length(batches) + 1]] <- batch
        }
      }
    } else {stop("Batch Step size is missing. Cannot perform sliding window analysis without step size.")}
    
  }
  return(batches)
}


save_list_elements <- function(.data, delimiter, file_extension, save_folder, file_name = "element_") {
  # Check if the save folder exists, if not, create it
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  # Use apply to save each element as a delimited file
  lapply(seq_along(.data), function(i) {
    filename <- file.path(save_folder, paste0(file_name, i, ".", file_extension))
    write.table(.data[[i]], file = filename, sep = delimiter, col.names = TRUE, row.names = FALSE)
  })
}

# # Function to save a batch to given directory or as a TEMP file
# save_batch <- function(batch_df, batch_index, batch_dir, file_extension = "csv", temp = FALSE) {
#   #library(VariantAnnotation)
#   
#   # Get directory from arguments
#   # Generate the file name with the specified file extension
#   batch_index <- 
#   batch_file_name <- paste0("batch.", batch_index, file_extension)
#   batch_file_path <- file.path(batch_dir, batch_file_name)
#   
#   #if fct specifies creating temp files
#   if (temp == TRUE){
#     batch_file_path <- tempfile(batch_file_name, batch_dir, fileext = file_extension)
#   }
#   #if saving non-temp files
#   else{
#     path <- file.path(batch_dir_name, batch_file_name) 
#   }
#   
#   # Save the batch dataframe with the specified file extension
#   if (file_extension == "csv") {
#     write.csv(batch_df$data, file = batch_file_path, row.names = FALSE)
#   } else if (file_extension == "vcf") {
#     writeVcf(batch_df$data, batch_file_path)
#   } else {
#     stop("Unsupported file extension. Use 'csv' or 'vcf'.")
#   }
# }

# Function to load or install required packages
load_install_pkg <- function(required_packages) {
  # Check if packages are installed, and install missing ones
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  }
  # Load all required packages
  lapply(required_packages, require, character.only = TRUE)
}

# `%notin%` <- function(x, y) {
#   !(x %in% y)
# }

load_data <- function(names = c("data"), ...){
  # Determine the number of arguments
  num_files <- length(as.list(...))
  num_names <- length(names)
  if (num_files != num_names){
    stop(paste0("You provided ", num_files, " files but only ", num_names, " names. There should be a name matching each file to be loaded"))
  }
  # Initialize an empty list to store the results
  out <- vector("list", length = num_files)
  names(out) <- names
  # Iterate through each argument
  for (i in seq_along(as.list(...))) {
    #get file type
    file_ext <- tools::file_ext(as.list(...)[[i]])
    if (file_ext %in% c("csv", "txt")) {
      # Read data from .csv or .txt file
      table <- read.csv(as.list(...)[[i]], header = TRUE, sep = ",")
    } else if (file_ext %in% c("tsv")) {
      # Read data from .csv or .txt file
      table <- read.csv(as.list(...)[[i]], header = TRUE, sep = "\t")
    }
    else {
      # Unsupported file type
      stop(paste0(as.list(...)[[i]]," is an unsupported file type. Only .csv, .txt are supported."))
    }
    # Store the result in the list
    out[[i]] <- table
  }
  return(out)
}
