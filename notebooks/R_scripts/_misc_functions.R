# Miscellanous functions
#' This script contains function that are used
#' across multiple scripts notebooks 

# ______________________________________________________________________________
#                     Utility Functions
# ______________________________________________________________________________

# inverse of logical statement
`%nin%` <- Negate(`%in%`)

#' Function for executing shell command in R with STDOUT and STEDRR control
#' specifying TRUE for stdout_path redirects the output of the command 
#' into an R object
shell_do <- function(command_string, stdout_path = "", stderr_path = "") {
  inputs <- unlist(stringr::str_split(command_string, " "))
  system2(
    command = inputs[1],
    args = inputs[-1],
    stdout = stdout_path,
    stderr = stderr_path
  )
}


# Function to run any shell command using slurm and future.batchtools
slurm_shell_do <- function(cmd,
                           jobname = glue("slurm-shell-{get_time()}"),
                           working_dir = wkdir,
                           template_path = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
                           memory = "1G",
                           ncpus = 1,
                           walltime = 3600) {
  require(magrittr)
  require(future)
  require(future.batchtools)
  # Initiate future.batchtools backend for parallel processing
  future::plan(
    future.batchtools::batchtools_slurm,
    template = template_path,
    resources = list(
      name = jobname,
      memory = memory,
      ncpus = ncpus,
      walltime = walltime
    )
  )
  job %<-% shell_do(cmd)
}


get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}

# function to count ongoing slurm jobs
count_slurm_jobs <- function(params = c("-u", "jboktor")) {
  queue <- system2(
    command = "squeue",
    args = params,
    stdout = TRUE
  )
  length(queue) - 1
}

#' sleep until a condition is true
wait_until <- function(conditional, interval = 2) {
  # Keep looping until the condition is met
  while (!conditional()) {
    Sys.sleep(interval)
  }
}

check_slurm_overload <- function(njobs = 9999) {
  wait_until(function() {
    count_slurm_jobs() < njobs
  })
}

#' Function to rapidly extract the number of files matching a search string
#' while integrating a list of sample names.
#' Input:
#' search_pattern   - a string with glue syntax {.}, not yet glued
#'                      to indicate placement of sampleIDs
#' name_list        - list of names
#' nworkders        - number of threads to use for
#'                    parallel processing
#' Output:
#' A named list of samples with the number of files
#' matching the search critera for each sample
search_file_n <- function(search_pattern, name_list, nworkers = 2, ...) {
  plan(multisession, workers = nworkers)
  file_n <- name_list %>%
    purrr::set_names() %>%
    purrr::map(~glue(search_pattern)) %>%
    furrr::future_map(
      ~ system2(
        command = "ls",
        args = c(., "|", "wc", "-l"),
        stdout = TRUE
      ) %>% as.numeric(),
      .progress = TRUE
    )
  return(file_n)
  }

# chunking function from https://stackoverflow.com/a/16275428
chunk_func <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

has_error_message <- function(stderr_file) {
  lines <- readLines(stderr_file)
  for (line in lines) {
    if (grepl("error|ERR|quota exceeded", line)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

seriate_matrix_rows <- function(mat,
                                seriate_method = "HC_average",
                                dist_method = "euclidean",
                                nthreads = 8) {
  order <- mat %>%
    parallelDist::parDist(method = dist_method, threads = nthreads) %>%
    # stats::dist(method = dist_method) %>%
    seriation::seriate(method = seriate_method) %>%
    seriation::get_order()
  ranked_order <- rownames(mat)[order]
  return(ranked_order)
}
