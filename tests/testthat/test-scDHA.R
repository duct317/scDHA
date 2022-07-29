test_that("scDHA can run", {
  skip_on_os("windows")
  skip_on_cran()
  
  #Load example data (Goolam dataset)
  data('Goolam'); data <- t(Goolam$data); label <- as.character(Goolam$label)
  #Log transform the data 
  data <- log2(data + 1)
  #Generate clustering result, the input matrix has rows as samples and columns as genes
  result <- scDHA(data, ncores = 2, seed = 1)
  expect_is(result, 'list')
  expect_is(result$cluster, 'numeric')
  expect_equal(length(result$cluster), nrow(data))
})

test_that("scDHA.vis can run", {
  #Load preprocessing data
  data("Goolam_result")
  result <- Goolam_result
  #Generate 2D representation, the input is the output from scDHA function
  result <- scDHA.vis(result, ncores = 2, seed = 1)
  expect_is(result, 'list')
  expect_is(result$pred, 'matrix')
  expect_equal(length(result$pred), length(result$cluster)*2)
})


test_that("scDHA.pt can run", {
  #Load preprocessing data
  data("Goolam_result")
  result <- Goolam_result
  #Generate pseudo-time for each cell, the input is the output from scDHA function
  result <- scDHA.pt(result, start.point = 1, ncores = 2, seed = 1)
  expect_is(result, 'list')
  expect_is(result$pt, 'numeric')
  expect_equal(length(result$pt), length(result$cluster))
})

test_that("scDHA.class can run", {
  skip_on_os("windows")
  skip_on_cran()
  
  #Load example data (Goolam dataset)
  data('Goolam'); data <- t(Goolam$data); label <- as.character(Goolam$label)
  #Log transform the data 
  data <- log2(data + 1)
  #Split data into training and testing sets
  set.seed(1)
  idx <- sample.int(nrow(data), size = round(nrow(data)*0.75))
  train.x <- data[idx, ]; train.y <- label[idx]
  test.x <- data[-idx, ]; test.y <- label[-idx]
  #Predict the labels of cells in testing set
  prediction <- scDHA.class(train = train.x, train.label = train.y, test = test.x, 
                             ncores = 2, seed = 1)
  expect_is(prediction, 'character')
  expect_equal(length(prediction), nrow(test.x))
})