

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package onLoad and onAttach
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {

  # if (!"dbplyr" %in% loadedNamespaces()) {
  #   requireNamespace("dbplyr", quietly = TRUE)
  # }

  # Set package options
  IST_options <- list(
    IST.project = "IST_Project",
    IST.verbose = TRUE,
    IST.check_python = TRUE,
    IST.default_assay = "Spatial",
    IST.parallel_workers = 4,
    IST.temp_dir = tempdir()
  )

  # Set options if not already set
  for (i in seq_along(IST_options)) {
    op <- names(IST_options)[i]
    if (is.null(getOption(op))) {
      do.call(options, stats::setNames(list(IST_options[[i]]), op))
    }
  }

  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Get package version
  pkgVersion <- tryCatch(
    utils::packageVersion(pkgname),
    error = function(e) "unknown"
  )

  # Startup message
  packageStartupMessage(
    paste0(
      "\n",
      "==================================================\n",
      "  ISTools (Infectious Spatiotemporal Transcriptomics Tools)\n",
      "  Version ", pkgVersion, "\n",
      "==================================================\n",
      "\n",
      "Welcome to ISTools package for analyzing \n",
      "Infectious Spatiotemporal Transcriptomics.\n",
      "\n",
      "Key functionalities:\n",
      "  - Background correction (CorrectBackgroud)\n",
      "  - Spot detection (SpotDetect_Gene, SpotDetect_Geneset)\n",
      "  - Niche detection (NicheDetect_Lasso, NicheDetect_STS)\n",
      "  - Niche analysis (CalNicheComposition, CalNicheDEGs, CalNicheCellComm)\n",
      "  - Data conversion (annh5ad2rds, rds2annh5ad)\n",
      "\n",
      "For more information, see the package documentation.\n",
      "==================================================\n"
    )
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' IST package options
#'
#' Options that control the behavior of the IST package. These can be set via
#' \code{\link{options}}.
#'
#' @section Options:
#' \describe{
#'   \item{\code{IST.project}}{Default project name (default: "IST_Project")}
#'   \item{\code{IST.verbose}}{Whether to print verbose output (default: TRUE)}
#'   \item{\code{IST.check_python}}{Whether to check Python configuration (default: TRUE)}
#'   \item{\code{IST.default_assay}}{Default assay name (default: "Spatial")}
#'   \item{\code{IST.parallel_workers}}{Number of parallel workers (default: 4)}
#'   \item{\code{IST.temp_dir}}{Temporary directory (default: tempdir())}
#' }
#'
#' @name IST-options
#' @rdname IST-options
#'
#' @examples
#' \dontrun{
#' options(IST.verbose = FALSE)
#' options(IST.parallel_workers = 8)
#' getOption("IST.project")
#' }
NULL

#' Set IST package options
#'
#' Convenience function to set multiple IST package options at once.
#'
#' @param ... Named options to set. Valid options include:
#'   \itemize{
#'     \item \code{verbose}: Logical, print verbose output
#'     \item \code{parallel_workers}: Integer, number of parallel workers
#'     \item \code{project}: Character, default project name
#'     \item \code{default_assay}: Character, default assay name
#'     \item \code{check_python}: Logical, check Python configuration
#'   }
#'
#' @return NULL (invisibly)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Set multiple options at once
#' IST_options(verbose = FALSE, parallel_workers = 8)
#'
#' # Set project name
#' IST_options(project = "MyProject")
#' }
IST_options <- function(...) {
  args <- list(...)
  # Map user-friendly names to internal option names
  names_mapping <- c(
    "verbose" = "IST.verbose",
    "parallel_workers" = "IST.parallel_workers",
    "project" = "IST.project",
    "default_assay" = "IST.default_assay",
    "check_python" = "IST.check_python",
    "temp_dir" = "IST.temp_dir"
  )

  # Convert names
  new_args <- args
  for (user_name in names(args)) {
    if (user_name %in% names(names_mapping)) {
      internal_name <- names_mapping[user_name]
      names(new_args)[which(names(new_args) == user_name)] <- internal_name
    } else if (!grepl("^IST\\.", user_name)) {
      # If no prefix and not in mapping, add IST. prefix
      names(new_args)[which(names(new_args) == user_name)] <- paste0("IST.", user_name)
    }
  }

  do.call(options, new_args)
  invisible()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables and constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' IST package global variables
#'
#' Internal global variables and constants used throughout the IST package.
#'
#' @keywords internal
#' @noRd
.IST_globals <- new.env(parent = emptyenv())

.IST_globals$supported_formats <- c("StereoSeq", "Visium")
.IST_globals$supported_hosts <- c("human", "mouse", "unknown")
.IST_globals$supported_pathogens <- c("virus", "bacteria", "parasite", "unknown")
.IST_globals$score_methods <- c("AddModuleScore", "AUCell", "UCell", "MeanExp", "SumExp")
.IST_globals$blur_methods <- c("isoblur", "medianblur")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package dependencies and requirements
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check Python dependencies
#'
#' Internal function to verify Python packages.
#'
#' @param envname Character, name of the conda environment
#' @param packages Character vector, Python packages to install
#'
#' @return Logical, TRUE if successful
#'
#' @keywords internal
#' @noRd
.check_python_deps <- function(envname = "r-ist",
                               packages = c("scanpy", "squidpy", "matplotlib",
                                            "seaborn", "pandas", "numpy", "anndata")) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("Package 'reticulate' is required for Python integration")
    return(FALSE)
  }

  tryCatch({
    for (pkg in packages) {
      if (!reticulate::py_module_available(pkg)) {
        message("Python package '", pkg, "' not found. Attempting to install...")
        reticulate::py_install(pkg, envname = envname)
      }
    }
    return(TRUE)
  }, error = function(e) {
    warning("Failed to install Python packages: ", e$message)
    return(FALSE)
  })
}

#' Check BLAST installation
#'
#' Internal function to verify BLAST+ tools.
#'
#' @return Logical, TRUE if BLAST tools are available
#'
#' @keywords internal
#' @noRd
.has_blast <- function() {
  blast_tools <- c("blastp", "blastn", "makeblastdb")
  available <- sapply(blast_tools, function(tool) {
    nchar(Sys.which(tool)) > 0
  })

  if (all(available)) {
    return(TRUE)
  } else {
    missing <- names(available)[!available]
    warning("BLAST tools not found in PATH: ", paste(missing, collapse = ", "))
    return(FALSE)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package-wide utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get IST package information
#'
#' Returns information about the IST package including version, available functions,
#' and current package options.
#'
#' @return A list containing package information with the following components:
#' \itemize{
#'   \item \code{version}: Package version
#'   \item \code{functions}: List of available functions by category
#'   \item \code{supported_formats}: Supported spatial transcriptomics formats
#'   \item \code{supported_hosts}: Supported host species
#'   \item \code{supported_pathogens}: Supported pathogen types
#'   \item \code{score_methods}: Available scoring methods
#'   \item \code{blur_methods}: Available blurring methods
#'   \item \code{options}: Current package options
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' info <- IST_info()
#' print(info$version)
#' print(info$supported_formats)
#' }
IST_info <- function() {
  pkg_version <- tryCatch(
    utils::packageVersion("IST"),
    error = function(e) "unknown"
  )

  list(
    version = pkg_version,
    functions = list(
      conversion = c("as.IST", "as.Seurat.IST", "annh5ad2rds", "rds2annh5ad"),
      background = c("CorrectBackgroud"),
      spot_detection = c("SpotDetect_Gene", "SpotDetect_Geneset"),
      niche_detection = c("NicheDetect_Lasso", "NicheDetect_STS", "NicheDetect_Spot", "NicheExpand"),
      niche_analysis = c("CalNicheComposition", "CalNicheAggIndex", "CalNicheCoLoc",
                         "CalNicheDEGs", "CalNicheGeneCor", "CalNichePPI",
                         "CalNicheCellComm", "CalNicheGRN"),
      utility = c("GetInfo", "SetInfo", "GetMetaData", "AddMetaData",
                  "CreateSingleSampNiche", "CreateMultiSampNiche")
    ),
    supported_formats = .IST_globals$supported_formats,
    supported_hosts = .IST_globals$supported_hosts,
    supported_pathogens = .IST_globals$supported_pathogens,
    score_methods = .IST_globals$score_methods,
    blur_methods = .IST_globals$blur_methods,
    options = list(
      IST.project = getOption("IST.project"),
      IST.verbose = getOption("IST.verbose"),
      IST.default_assay = getOption("IST.default_assay"),
      IST.parallel_workers = getOption("IST.parallel_workers")
    )
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clean up on package unload
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onUnload <- function(libpath) {
  # Clean up temporary files
  temp_dir <- getOption("IST.temp_dir", tempdir())
  temp_files <- list.files(temp_dir, pattern = "^IST_.*", full.names = TRUE)
  unlink(temp_files, recursive = TRUE, force = TRUE)

  # Close any open connections
  while (sink.number() > 0) {
    sink()
  }

  invisible()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# R version compatibility
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check R version compatibility
#'
#' Internal function to verify R version meets minimum requirements.
#'
#' @param minimum_version Character, minimum required R version (default: "4.0.0")
#'
#' @return Logical, TRUE if compatible
#'
#' @keywords internal
#' @noRd
.check_r_version <- function(minimum_version = "4.0.0") {
  if (utils::compareVersion(as.character(getRversion()), minimum_version) < 0) {
    warning("This package requires R version ", minimum_version,
            " or higher. Current version: ", getRversion())
    return(FALSE)
  }
  return(TRUE)
}

# Run version check when package is loaded in interactive session
if (interactive()) {
  if (!.check_r_version()) {
    packageStartupMessage("Warning: IST package may not be fully compatible with this R version.")
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' IST package version
#'
#' @keywords internal
#' @noRd
IST_VERSION <- tryCatch(
  utils::packageVersion("IST"),
  error = function(e) "unknown"
)

#' IST package URL
#'
#' @keywords internal
#' @noRd
IST_URL <- "https://github.com/yourusername/IST"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# End of zzz.R
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
