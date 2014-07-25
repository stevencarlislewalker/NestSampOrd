
mkLike <- function(Y) {
    local({
        Y <- Y
        n <- nrow(Y)
        m <- ncol(Y)
        thetaFun <- function(x, a, u, t) {
            theta <- c(x, a, u, t)
            setNames(theta, c(paste0("x", 1:n),
                              paste0("a", 1:m),
                              paste0("u", 1:m),
                              paste0("t", 1:m)))
        }
        thetaInv <- function(theta) {
            list(x = theta[1:n],
                 a = theta[(n+1):(n+m)],
                 u = theta[(n+m+1):(n+2*m)],
                 t = theta[(n+2*m+1):(n+3*m)])
        }
        etaFun <- function(theta, n, m) {
            with(thetaInv(theta), {
                matrix(a, n, m, byrow = TRUE) - 
                    matrix(t^(-2), n, m, byrow = TRUE) *
                        (outer(x, u, "-")^2)
            })
        }
        rPrior <- function() {
            thetaFun(rnorm(n),
                     rnorm(m),
                     rnorm(m),
                     abs(rnorm(m)))
        }
        dPrior <- function(theta) {
            -0.5*sum(dnorm(unlist(theta), log = TRUE))
        }
        rProp <- function() {
            thetaFun(rnorm(n, sd = 0.1),
                     rnorm(m, sd = 0.1),
                     rnorm(m, sd = 0.1),
                     abs(rnorm(m, sd = 0.1)))
        }
        theta <- thetaFun(scores(princomp(Y))[,1],
                          qlogis(apply(Y, 2, mean)),
                          rep(0, m), rep(1, m))
        dev <- binomial()$aic
        linkinv <- binomial()$linkinv
        function(theta) {
            dev(Y, 1, linkinv(etaFun(theta, n, m)), 1)
        }
    })
}


library(reo)
Y <- as.matrix(fish)
lfun <- mkLike(Y)
rho <- environment(lfun)
theta <- rho$theta
lfun(theta)
## thetaList <- rho$thetaInv(theta)
## thetaList$x[] <- rnorm(52)
## theta <- do.call(rho$thetaFun, thetaList)
## opt <- optim(theta, lfun, method = "BFGS",
##              control = list(maxit = 500,
##                  trace = 10))
## opt
## rho$thetaInv(opt$par)

theta0 <- replicate(100, rho$rPrior())
lfun0 <- apply(theta0, 2, lfun)
q0 <- quantile(lfun0, c(1/exp(1), 1-1/exp(1)))

lev <- 0


theta1 <- theta0
lfun0[wm <- which.min(abs(lfun0 - q0[2]))]
theta1[, 1] <- theta0[, wm]


for(i in 2:100) {
    if(runif(1) < 0.5) lev <- 1*ifelse(lev, 0, 1)
    print("iter"); print(i)
    print("lev"); print(lev)
    if(lev == 0) {
        theta1[, i] <- rho$rPrior()
    } else {
        prop <- theta1[, i-1] + rho$rProp()
        if(lfun(prop) > q0[1]) {
            aprob <- min(1, exp(rho$dPrior(prop) -
                                rho$dPrior(theta1[, i-1])))
        } else {
            aprob <- 0
        }
        uu <- runif(1)
        print("aprob, uu"); print(aprob); print(uu)
        if(uu < aprob) {
            theta1[, i] <- prop
        } else {
            theta1[, i] <- theta1[, i-1]
        }
    }
}
    
    
    

lfun1 <- apply(theta1, 2, lfun)
plot(1:100, lfun1, type = "l")
plot(1:100, apply(theta1, 2, rho$dPrior), type = "l")
