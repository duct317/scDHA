library(testthat)
library(scDHA)

if((!torch::torch_is_installed()) & (tolower(Sys.info()[["sysname"]]) != "windows"))
{
  torch::install_torch()
}

Sys.setenv(IN_TEST_scDHA="true")

test_check("scDHA")

Sys.unsetenv("IN_TEST_scDHA")
