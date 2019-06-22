#' gdp penalty function
#'
#' This function define the gdp penalty. The
#' formula is \code{gdp(x) = lambda * log( 1 + x/gamma)}.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}.
#'
#' @param x a non-negative numeric value
#' @param gamma hyper-tuning parameter
#' @param lambda tuning parameter
#'
#' @return the value of gdp penalty
#'
#' @examples
#' \dontrun{gdp(0:9, gamma = 1, lambda = 5)}
gdp <- function(x, gamma, lambda) {
    # gdp penalty for non-negative values x>=0
    
    # the domain of gdp() should be in [0,+inf)
    stopifnot(x >= 0)
    
    y <- lambda * log(1 + x/gamma)
    y
}

#' Super gradient of gdp penalty function
#'
#' This function define the super gradient of the gdp penalty.
#' The formula is \code{gdp_sg(x) = lambda/(gamma + x)}.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}
#'
#' @inheritParams gdp
#'
#' @return the super-gradient of the gdp penalty function
#'
#' @examples
#' \dontrun{gdp_sg(0:9,gamma=1,lambda=5)}
gdp_sg <- function(x, gamma, lambda) {
    # supergradient of gdp penalty x>=0
    
    # the domain of gdp_sg() should be in [0,+inf)
    stopifnot(x >= 0)
    
    y <- lambda/(gamma + x)
    y
}

#' Lq penalty function
#'
#' This function define the Lq penalty. The
#' formula is \code{Lq(x) =lambda*x^q}.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}
#'
#' @inheritParams gdp
#'
#' @return the value of the Lq penalty
#'
#' @examples
#' \dontrun{lq(0:9,gamma=0.5,lambda=5)}
lq <- function(x, gamma, lambda) {
    # lq penalty
    
    # the domain of lq() should be in [0,+inf)
    stopifnot(x >= 0)
    
    # Lq is non-differential at 0
    epsilon <- 0.5 * .Machine$double.eps
    y <- lambda * ((x + epsilon)^gamma)
    y
}

#' Super gradient of the Lq penalty function
#'
#' This function define the super gradient of the Lq penalty.
#' The formula is \code{lq_sg(x) = lambda * q * x^(p-1)}.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}
#'
#' @inheritParams lq
#'
#' @return the super-gradient of the Lq penalty function
#'
#' @examples
#' \dontrun{
#' lq_sg(0:9,gamma=1,lambda=5)
#' lq_sg(0:9,gamma=0.5,lambda=5)
#' }
lq_sg <- function(x, gamma, lambda) {
    # supergradient of lq penalty
    
    # the domain of lq() should be in [0,+inf)
    stopifnot(x >= 0)
    
    epsilon <- 0.5 * .Machine$double.eps
    y <- lambda * gamma * ((x + epsilon)^(gamma - 1))
    y
}

#' SCAD penalty function
#'
#' This function define the SCAD penalty.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}
#'
#' @inheritParams gdp
#'
#' @return the value of SCAD penalty
#'
#' @examples
#' \dontrun{scad(0:9,gamma=3.7,lambda=3)}
scad <- function(x, gamma, lambda) {
    
    # the domain of scad() should be in [0,+inf)
    stopifnot(x >= 0)
    
    y <- rep(1, length(x)) * (gamma + 1) * (lambda^2/2)
    
    ind <- (x <= lambda)
    y[ind] <- lambda * x[ind]
    
    ind <- (x > lambda & x <= gamma * lambda)
    y[ind] <- (-x[ind]^2 + 2 * gamma * lambda * x[ind] - lambda^2)/(2 * (gamma - 1))
    
    y
}

#' Super gradient for the SCAD penalty function
#'
#' This function define the supergradient of the SCAD penalty.
#' Details can be found in \url{https://arxiv.org/abs/1807.04982}
#'
#' @inheritParams scad
#'
#' @return the super gradient of the SCAD penalty
#'
#' @examples
#' \dontrun{scad_sg(0:9,gamma=3.7,lambda=3)}
scad_sg <- function(x, gamma, lambda) {
    # the domain of scad() should be in [0,+inf)
    stopifnot(x >= 0)
    
    y <- rep(0, length(x))
    ind <- (x <= lambda)
    y[ind] <- lambda
    
    ind <- (x > lambda & x <= gamma * lambda)
    y[ind] <- (gamma * lambda - x[ind])/(gamma - 1)
    
    y
}

