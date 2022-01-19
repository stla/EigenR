dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "EigenR.dll", package = "EigenR")
  )
}

makedoc <- function(){
  roxygen2::roxygenise(load_code = roxygen2::load_installed)
}
