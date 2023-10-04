.onAttach <- function(libname, pkgname) {
  version <- "1.1.0"
  art <- c(
    "   ___                _____  ___  ____ ", "  / _ \\___ ____  ___ / / _ \\/ _ \\/ __ \\",
    " / ___/ _ `/ _ \\/ -_) / ___/ , _/ /_/ /", "/_/   \\_,_/_//_/\\__/_/_/  /_/|_|\\____/ "
  )
  yellow <- crayon::yellow
  msg <- sprintf("PPP version: %s", version)
  packageStartupMessage("Welcome to the PPP package, by the BayesMendel Lab!")
  packageStartupMessage(msg)
  packageStartupMessage(paste(yellow(art), "\n"))
}
