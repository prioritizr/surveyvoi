# Initialization
## load packages
library(knitr)

## set variables
start_seed <- 100
max_seed <- 1000

# Preliminary processing
## purl R code
file_path <- tempfile(fileext = ".R")
print(file_path)
knitr::purl(input = "vignettes/surveyvoi.Rmd", output = file_path)

# Main processing
##  initialize loop variables
curr_seed <- start_seed
success <- FALSE

## main loop
while ((curr_seed <= max_seed) && (!success)) {
  ## print information
  message(paste0("starting seed: ", curr_seed))

  ## update seed in R script
  code <- readLines(file_path)
  idx <- which(startsWith(code, "seed <- "))
  code[idx] <- paste0("seed <- ", curr_seed)
  writeLines(code, file_path)

  ## try to run R script
  out <- capture.output(suppressMessages({
    result <- try(source(file_path), silent = TRUE)
  }))

  ## determine if successful
  if (!inherits(result, "try-error")) {
    success <- TRUE
  }

  ## increment seed if not successful
  if (!success) {
    curr_seed <- curr_seed + 1
  }
}

# Exports
## print seed
if (success) {
  message("success with seed:")
  message(curr_seed)
} else {
  message("Failed to find valid seed")
}
