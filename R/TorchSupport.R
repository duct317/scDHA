#' @import torch
scDHA_dataset <- dataset(
  
  name = "scDHA_dataset",
  
  initialize = function(input) {
    self$data <- self$prepare_scDHA_data(input)
  },
  
  .getbatch = function(index) {
    
    x <- self$data[index, ]
    list(x)
  },
  
  .length = function() {
    self$data$size()[[1]]
  },
  
  prepare_scDHA_data = function(input) {
    input <- as.matrix(input)
    torch_tensor(input, dtype = torch_float())
  }
)

scDHA_AE <- nn_module(
  "scDHA_AE",
  initialize = function(original_dim, im_dim) {
    self$fc1 <- nn_linear(original_dim, im_dim)
    self$fc2 <- nn_linear(im_dim, original_dim)
    
    nn_init_zeros_(self$fc1$bias)
    nn_init_zeros_(self$fc2$bias)
    nn_init_xavier_uniform_manual(self$fc1$weight)
    nn_init_xavier_uniform_manual(self$fc2$weight)
  },
  forward = function(x) {
    x %>% 
      self$fc1() %>% 
      self$fc2() 
  }
)

sampling <- function(mu, var, epsilon_std, lat_dim)
{
  mu + torch_sqrt(var)*epsilon_std*torch_randn(c(1, lat_dim))
}

scDHA_VAE <- nn_module(
  "scDHA_VAE",
  initialize = function(original_dim, im_dim, lat_dim, epsilon_std, batch_norm=TRUE, zero_bias=TRUE) {
    self$batch_norm <- batch_norm
    self$epsilon_std <- epsilon_std
    self$lat_dim <- lat_dim
    if(batch_norm)
    {
      self$h <- nn_linear(original_dim, im_dim)
      self$bn <- nn_batch_norm1d(im_dim, momentum = 0.01, eps = 1e-3)
      
    } else {
      self$h <- nn_linear(original_dim, im_dim)
    } 
    
    self$mu <- nn_linear(im_dim, lat_dim)
    self$var <- nn_linear(im_dim, lat_dim)
    
    self$h1 <- nn_module_list(list(nn_linear(lat_dim, im_dim),
                                   nn_linear(lat_dim, im_dim)
                                    ))
    self$x_ <- nn_module_list(list(nn_linear(im_dim, original_dim),
                                   nn_linear(im_dim, original_dim)
                                    ))
    
    nn_init_xavier_uniform_manual(self$h$weight)
    nn_init_xavier_uniform_manual(self$mu$weight)
    nn_init_xavier_uniform_manual(self$var$weight)
    nn_init_xavier_uniform_manual(self$h1[[1]]$weight)
    nn_init_xavier_uniform_manual(self$h1[[2]]$weight)
    nn_init_xavier_uniform_manual(self$x_[[1]]$weight)
    nn_init_xavier_uniform_manual(self$x_[[2]]$weight)
    
    if(zero_bias)
    {
      nn_init_zeros_(self$h$bias)
      nn_init_zeros_(self$mu$bias)
      nn_init_zeros_(self$var$bias)
      nn_init_zeros_(self$h1[[1]]$bias)
      nn_init_zeros_(self$h1[[2]]$bias)
      nn_init_zeros_(self$x_[[1]]$bias)
      nn_init_zeros_(self$x_[[2]]$bias)
    }
  },
  forward = function(x) {
    if(self$batch_norm)
    {
      im <- x %>% 
        self$h() %>% self$bn()
      mu <- im %>% self$mu()
      var <- im %>% self$var() %>% nnf_softmax(dim = 2)
    } else {
      im <- x %>% 
        self$h() %>% nnf_selu()
      mu <- im %>% self$mu()
      var <- im %>% self$var() %>% nnf_softmax(dim = 2)
    }
    out <- list()
    out[[1]] <- mu
    out[[2]] <- var
    for (i in 1:2) {
      z <- sampling(mu, var, self$epsilon_std, self$lat_dim)
      out[[2+i]] <- z %>% 
        self$h1[[i]]() %>% nnf_selu() %>% self$x_[[i]]()
    }
    out
  },
  encode_mu = function(x)
  {
    if(self$batch_norm)
    {
      im <- x %>% 
        self$h() %>% self$bn()
      mu <- im %>% self$mu()
    } else {
      im <- x %>% 
        self$h() %>% nnf_selu()
      mu <- im %>% self$mu()
    }
    mu
  }
)

nn_init_xavier_uniform_manual <- function(tensor, gain = 1) {
  fans <- nn_init_calculate_fan_in_and_fan_out(tensor)
  fan_in <- fans[[1]]
  fan_out <- fans[[2]]
  std <- gain * sqrt(2.0 / (fan_in + fan_out))
  a <- sqrt(3.0) * std # Calculate uniform bounds from standard deviation
  nn_init_no_grad_uniform(tensor, -a, a)
}

nn_init_calculate_fan_in_and_fan_out <- function(tensor) {
  
  dimensions <- tensor$dim()
  num_input_fmaps <- tensor$size(2)
  num_output_fmaps <- tensor$size(1)
  receptive_field_size <- 1
  
  # if (dimensions > 2)
  #   receptive_field_size <- tensor[1,1,..]$numel()
  
  fan_in <- num_input_fmaps * receptive_field_size
  fan_out <- num_output_fmaps * receptive_field_size
  
  list(fan_in, fan_out)
}

nn_init_no_grad_uniform <- function(tensor, a, b) {
  with_no_grad({
    out <- tensor$uniform_(a, b)
  })
  out
}

state <- function(self) {
  attr(self, "state")
}

`state<-` <- function(self, value) {
  attr(self, "state") <- value
  self
}

optim_adamw <- optimizer(
  "optim_adamw", 
  initialize = function(params, lr=1e-3, betas=c(0.9, 0.999), eps=1e-8,
                        weight_decay=1e-2) {
    
    
    if (lr < 0)
      value_error("Invalid learning rate: {lr}")
    
    if (eps < 0)
      value_error("Invalid eps: {eps}")
    
    if (betas[[1]] < 0 || betas[[1]] > 1)
      value_error("Invalid beta parameter at index 1")
    
    if (betas[[2]] < 0 || betas[[2]] > 1)
      value_error("Invalid beta parameter at index 2")
    
    if (weight_decay < 0)
      value_error("Invalid weight decay value: {weight_decay}")
    
    defaults <- list(lr=lr, betas=betas, eps = eps, weight_decay = weight_decay)
    
    super$initialize(params, defaults)
  },
  
  step = function(closure = NULL) {
    loop_fun <- function(group, param, g, p) {
      
      grad <- param$grad
      
      # state initialization
      if (length(state(param)) == 0) {
        state(param) <- list()
        state(param)[["step"]] <- 0
        state(param)[["exp_avg"]] <- torch_zeros_like(param, memory_format=torch_preserve_format())
        state(param)[["exp_avg_sq"]] <- torch_zeros_like(param, memory_format=torch_preserve_format())
      }
      
      # Perform stepweight decay
      param$mul_(1 - group$lr*group$weight_decay)
      
      exp_avg <- state(param)[["exp_avg"]]
      exp_avg_sq <- state(param)[["exp_avg_sq"]]

      beta1 <- group$betas[[1]]
      beta2 <- group$betas[[2]]
      
      state(param)[["step"]] <- state(param)[["step"]] + 1
      bias_correction1 <- 1 - beta1 ^ state(param)[['step']]
      bias_correction2 <- 1 - beta2 ^ state(param)[['step']]
      
      # Decay the first and second moment running average coefficient
      exp_avg$mul_(beta1)$add_(grad, alpha=1 - beta1)
      exp_avg_sq$mul_(beta2)$addcmul_(grad, grad, value=1 - beta2)

      denom <- (exp_avg_sq$sqrt() / sqrt(bias_correction2))$add_(group$eps)

      step_size <- group$lr / bias_correction1
      
      param$addcdiv_(exp_avg, denom, value=-step_size)
    }
    private$step_helper(closure, loop_fun)
  }
)


cal.dis <- function(arg) {
  a <- arg
  a2 <- torch_square(a)
  a2sum <- torch_sum(a2, dim = 2, keepdim = TRUE)
  dis <- a2sum$t() + a2sum
  ab <- torch_matmul(a, a$t())
  final.dis <- dis - 2 * ab
  final.dis$fill_diagonal_(0)
  final.dis <- torch_sqrt(final.dis+1e-6)
  final.dis
}


scDHA_model_vis <- nn_module(
  "scDHA_model_vis",
  initialize = function(original_dim, out_dim, eps = 0.1) {
    self$eps <- eps
    self$h1 <- nn_linear(original_dim, 1024)
    self$h2 <- nn_linear(1024, 32)
    
    self$z <- nn_linear(32, 2)
    self$bn <- nn_batch_norm1d(2, momentum = 0.1, eps = 1e-3)
    
    self$x_ <- nn_linear(2, out_dim)
    
    nn_init_zeros_(self$h1$bias)
    nn_init_zeros_(self$h2$bias)
    nn_init_zeros_(self$z$bias)
    nn_init_zeros_(self$x_$bias)
    nn_init_xavier_uniform_manual(self$h1$weight)
    nn_init_xavier_uniform_manual(self$h2$weight)
    nn_init_xavier_uniform_manual(self$z$weight)
    nn_init_xavier_uniform_manual(self$x_$weight)
  },
  forward = function(x) {
    if(self$training)
    {
      x <- x + torch_randn_like(x)*self$eps
    }
    x <- x %>% self$h1() %>% nnf_elu() %>% nnf_dropout()
    x <- x %>% self$h2() %>% nnf_sigmoid()
    z <- x %>% self$z() %>% self$bn()
    
    z1 <- cal.dis(z)
    
    x_ <- z %>% self$x_() %>% nnf_softmax(dim = 2)
    
    list(z1, x_)
  },
  
  encode_latent = function(x) {
    x <- x %>% self$h1() %>% nnf_elu() %>% nnf_dropout()
    x <- x %>% self$h2() %>% nnf_sigmoid()
    z <- x %>% self$z() %>% self$bn()
    
    z
  }
)

check_grad_nan <- function(parameters) {
  if (is(parameters, "torch_tensor"))
    parameters <- list(parameters)
  parameters <- Filter(function(x) !is_undefined_tensor(x$grad), parameters)
  for (p in parameters) {
    if(p$grad$sum()$isnan()$item()) return(TRUE)
  }
  FALSE
}