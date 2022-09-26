
env_pkgs <- .packages(all.available = T)
attached <- (.packages())
req_libs <- c(
  "emo",
  "logger",
  "devtools"
)
to_install <- req_libs[!req_libs %in% env_pkgs]
failed_pkgs <- c()


if (!length(to_install)) {
  lapply(req_libs, library, character.only = TRUE)
} else {
  logger::log_info(sprintf("installing missing package %s ...", to_install))
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

# github
# install.packages("devtools")
devtools::install_github("hadley/emo")
