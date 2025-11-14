# Robust installer for: exomePeak, exomePeak2, TRESS, MeTPeak
# Usage:
#   install_mystery_packages()
#   install_mystery_packages(update = TRUE)
#   install_mystery_packages(github_map = c(exomePeak2 = "username/exomePeak2"))
install_mystery_packages <- function(pkgs = c("exomePeak", "exomePeak2", "TRESS", "MeRIPtools","foreach","tidyverse","data.table"),
                                     update = FALSE,
                                     ask = FALSE,
                                     github_map = NULL,
                                     bioc_install_args = list()) {
  # pkgs: character vector of package names to ensure installed
  # update: if TRUE, force update (passes update = TRUE to BiocManager::install when used)
  # ask: whether to ask before installing Bioconductor packages (passed to BiocManager::install)
  # github_map: named character vector mapping package -> "user/repo" for GitHub installation fallback
  # bioc_install_args: additional named args passed to BiocManager::install()
  #
  # Returns: invisibly a named logical vector indicating which packages are installed (TRUE) or not (FALSE).

  stopifnot(is.character(pkgs), length(pkgs) >= 1)
  if (!is.null(github_map)) {
    if (!is.character(github_map) || is.null(names(github_map))) {
      stop("github_map must be a named character vector like c(pkgname = 'owner/repo').")
    }
  }

  installed_before <- rownames(utils::installed.packages())
  need_install <- setdiff(pkgs, installed_before)

  #install MeTPeak from source
  if (length(need_install) > 0 ) {
    if(length(grep(need_install, pattern="MeTPeak"))>0){
      file_path <- system.file("extdata", "MeTPeak-master", package = "mysterypackage")
      if(file_path != ""){
        install.packages(file_path,repos=NULL,type="source")
      }else{
        message("Required source package for MeTPeak is not found")
      }
    }
    #install remaining packages
    installed_before <- rownames(utils::installed.packages())
    need_install <- setdiff(pkgs, installed_before)
  }
  # #install exomePeak from source
  # if (length(need_install) > 0 ) {
  #   if(length(grep(need_install, pattern="exomePeak") )>0){
  #     file_path <- system.file("extdata", "exomePeak-master", package = "mysterypackage")
  #     if(file_path != ""){
  #       install.packages(file_path,repos=NULL,type="source")
  #     }else{
  #       message("Required source package for exomePeak is not found")
  #     }
  #   }
  #   #install remaining packages
  #   installed_before <- rownames(utils::installed.packages())
  #   need_install <- setdiff(pkgs, installed_before)
  # }
  # #install GenomicFeatures<1.61 from source
  # if (length(need_install) > 0 ) {
  #   if(length(grep(need_install, pattern="GenomicFeatures"))>0){
  #     file_path <- system.file("extdata", "GenomicFeatures_1.56.0.tar.gz", package = "mysterypackage")
  #     if(file_path != ""){
  #       install.packages(file_path,repos=NULL,type="source")
  #     }else{
  #       message("Required source package for GenomicFeatures_1.56 is not found")
  #     }
  #   }
  #   #install remaining packages
  #   installed_before <- rownames(utils::installed.packages())
  #   need_install <- setdiff(pkgs, installed_before)
  # }
  # Ensure BiocManager is available for Bioconductor packages and general robust install
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager from CRAN...")
    install.packages("BiocManager")
  }

  # First attempt: use BiocManager::install for all targets (BiocManager can install CRAN + Bioconductor)
  if (length(need_install) > 0 || update) {
    message("Attempting installation via BiocManager::install(): ",
            paste(pkgs, collapse = ", "))
    bioc_args <- c(list(pkgs, ask = ask, update = update), bioc_install_args)
    # Use do.call to forward additional args
    tryCatch({
      do.call(BiocManager::install, bioc_args)
    }, error = function(e) {
      message("BiocManager::install failed: ", conditionMessage(e))
    })
  }

  # Re-check installed status
  installed_now <- rownames(utils::installed.packages())
  still_missing <- setdiff(pkgs, installed_now)

  # Second attempt: try CRAN via install.packages() for remaining
  if (length(still_missing) > 0) {
    message("Attempting installation from CRAN for: ", paste(still_missing, collapse = ", "))
    tryCatch({
      utils::install.packages(still_missing)
    }, error = function(e) {
      message("install.packages failed: ", conditionMessage(e))
    })
  }

  # Re-check
  installed_now <- rownames(utils::installed.packages())
  still_missing <- setdiff(pkgs, installed_now)

  # Third attempt: GitHub fallback for packages with provided github_map
  if (length(still_missing) > 0 && !is.null(github_map)) {
    # Ensure remotes is present
    if (!requireNamespace("remotes", quietly = TRUE)) {
      message("Installing remotes from CRAN for GitHub installs...")
      install.packages("remotes")
    }
    for (pkg in still_missing) {
      if (!is.null(github_map[[pkg]])) {
        repo <- github_map[[pkg]]
        message(sprintf("Attempting to install %s from GitHub repo '%s' ...", pkg, repo))
        tryCatch({
          remotes::install_github(repo)
        }, error = function(e) {
          message("remotes::install_github failed for ", pkg, ": ", conditionMessage(e))
        })
      }
    }
  }

  # Final status
  final_installed <- pkgs %in% rownames(utils::installed.packages())
  names(final_installed) <- pkgs

  if (all(final_installed)) {
    message("All requested packages are installed.")
  } else {
    missing_pkgs <- names(final_installed)[!final_installed]
    message("The following packages are still missing: ", paste(missing_pkgs, collapse = ", "))
    message("If any are hosted on GitHub, provide a github_map like:\n  github_map = c(exomePeak2 = 'owner/exomePeak2')")
  }

  invisible(final_installed)
}

# Example:
# install_mystery_packages()
# If you know the GitHub repo for a package not on CRAN/Bioconductor:
# install_mystery_packages(github_map = c(exomePeak2 = "user/exomePeak2"))

