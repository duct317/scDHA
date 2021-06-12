library(testthat)
library(scDHA)

on_cran <- function() !identical(Sys.getenv("NOT_CRAN"), "true")

if((!torch::torch_is_installed()) & 
   (tolower(Sys.info()[["sysname"]]) != "windows") & 
   (!on_cran()))
{
  torch::install_torch()
}

Sys.setenv(IN_TEST_scDHA="true")

test_check("scDHA")

Sys.unsetenv("IN_TEST_scDHA")
