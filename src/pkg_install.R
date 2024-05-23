
env_pkgs <- .packages(all.available = T)
attached <- (.packages())
req_libs <- c(
  "jsonlite",
  "languageserver",
  "pandoc",
  "devtools",
  "NADA2",
  "ssizeRNA"
)

others <-c(
    "logger",
  "devtools",
  "lintr",
  "ICC.Sample.Size",
  "CVcalibration",
  "kableExtra",
  "dplyr",
  "tidyr",
  "pcr")

to_install <- req_libs[!req_libs %in% env_pkgs]
failed_pkgs <- c()

getOption("repos")
options(
  repos = c(
    CRANextra = "https://macos.rbind.io",
    CRANstudio = "https://cran.rstudio.com",
    CRAN = "https://CRAN.R-project.org"
  ),
  download.file.method = "wget"
)

if (!length(to_install)) {
  lapply(req_libs, library, character.only = TRUE)
} else {
  message(sprintf("installing missing package %s ...", to_install))
  tryCatch(
    {
      progress <- txtProgressBar(0, length(to_install), 3)
      lapply(seq_along(to_install), function(i) {
        pkg <- to_install[i]
        install.packages(pkg,
          dependencies = TRUE,
          INSTALL_opts = c("--no-lock")
        )
        setTxtProgressBar(progress, i)
      })
    },
    error = function(e) {
      print(e)
      failed_pkgs <<- c(failed_pkgs, pkg)
    },
    finally = {
      if (length(failed_pkgs)) warning(sprintf("Failed to install packages: %s", failed_pkgs))
      lapply(req_libs, library, character.only = TRUE)
    }
  )
}

BiocManager::install("ddCt")
BiocManager::install("RnaSeqSampleSize")
BiocManager::install("RNASeqPower")
BiocManager::install("ssizeRNA")
# github
devtools::install_github("hadley/emo")
