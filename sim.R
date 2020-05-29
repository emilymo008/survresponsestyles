# Simulation code, for running on multi-core machine.

library(lavaan)
library(parallel)
suppressMessages(library(foreach))
library(iterators)
library(doMC)
library(moments)

round2 <- function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

num_cores <- detectCores()
registerDoMC(num_cores)

set.seed(18)
newseeds <- NULL
for (i in 1:num_cores){
  tmp.ind <- floor(runif(1,2,627)) 
  next.rand <- .Random.seed[tmp.ind]
  newseeds <- c(newseeds,next.rand)
}

sims <- 1000
core.sims <- c(rep(times = num_cores-1, x = sims %/% num_cores), sims - (num_cores-1)*(sims %/% num_cores))
sum(core.sims) == sims


{
  shift.mat <- list()  #### this contains the no-response style tau cuts!
  shift.mat[[1]] <- c(-1.0, 0, 1.0)  # gran 4 --> index is 1
  shift.mat[[2]] <- c(-1.5, -0.5, 0.5, 1.5)  # gran 5 --> index is 2
  shift.mat[[3]] <- c(-1.6, -0.8, 0.0, 0.8, 1.6)  # gran 6 --> index is 3
  shift.mat[[4]] <- c(-2.0, -1.2, -0.4, 0.4, 1.2, 2.0)  # gran 7 --> index is 4
  shift.mat[[5]] <- c(-1.8, -1.4, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4, 1.8)  # gran 11 --> index is 5
}




#################################################################################################################################################

par.list <- foreach(i = 1:num_cores) %dopar% {
  set.seed(newseeds[i])
  count <- 0
  # note: We originally intended to look at percent contamination from 0 to 39 for granularity 4 to 7, and percent contamination from 0 to 60 for granularity 11. In final analysis, decided to limit it to 0 to 39 for all granularity levels, for the sake of having a completely balanced design.
  rows <- 11*(4*14 + 1*21)*core.sims[i]  # [11 = largest possible number of skews]*[(4 = gran 4 to 7)*(14 = number of variations of pct for gran 4 to 7) + (1 = gran 11)*(21 = number of variations of pct for gran 11)] 
  p <- 0.3  # actual pre-test effect
  t <- 0.5  # actual treatment effect
  est <- array(NA, dim = c(rows, 12))  # set up array to put results
  
  for (k in c(15)) {  # manifest items
    txt2.x <- paste('x', 1:k, sep = '')
    txt2.x <- paste(txt2.x, collapse = ' + ')
    txt2.y <- paste('y', 1:k, sep = '')
    txt2.y <- paste(txt2.y, collapse = ' + ')
    lmod <- paste("pre =~", txt2.x, '\n', "post =~", txt2.y, '\n', "post ~ trt + pre")  # make SEM model
    for (b in c(4, 5, 6, 7, 11)) {  # granularity levels
        if (b %in% c(4, 5, 6, 7)) {
          gind <- b-3  # index for shift.mat
        }
        else if (b == 11) {
          gind <- 5
        }
      
      for (n in c(250)) {  # n (number of respondents)
        for (g in c(0, -0.2, -0.4, -0.6, -0.8, -1, 0.2, 0.4, 0.6, 0.8, 1)) {  # indices for skews
          std.tc <- shift.mat[[gind]] + g  # no response style + skew
          cnt.tc <- rep(0, b-1) + g # response style + skew: make a vector of all 0s, then add skew
          for (j in c(0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36, 0.39,
                      0.42, 0.45, 0.48, 0.51, 0.54, 0.57, 0.60)) {  # percent contamination levels
            
            u <- round2(j*n, 0)  # number of data entries that get contaminated
            v <- n-u  # number of data entries that don't get contaminated
            
            for (h in 1:core.sims[i]) {  # number of simulations
              count <- count + 1
              
              # pre test:
              pre <- array(c(rnorm(n = n, mean = 0, sd = 1)), dim = c(n, 1))  # construct for each respondent
              
              lambda <- array(runif(min = 0.45, max = 0.85, n = k), dim = c(k, 1))  # kx1, factor loadings for each manifest item
              eps <- array(rnorm(mean = 0, sd = (1-lambda^2)^(1/2), n = n*k), dim = c(k, n))  # kxn 
              pre.star <- t(lambda %*% t(pre)) + t(eps)
              pre.itms <- as.data.frame(pre.star)  # pretest uncategorized (continuous) answers
              colnames(pre.itms) <- paste('x', 1:k, sep = '')
              if (j != 0) {
                cnt.pre.tc <- array(-2, dim=c(u, k))  # array for contaminated respondents
                std.pre.tc <- array(-2, dim=c(v, k))  # array for uncontaminated respondents
                # how this works: for each tau cut, adds 1 starting from lowest option if answer greater than tau cut
                for (ii in 1:(b-1)) {cnt.pre.tc <- cnt.pre.tc + (pre.star[(1:u),] > cnt.tc[ii])*1}  # build cont. array
                for (ii in 1:(b-1)) {std.pre.tc <- std.pre.tc + (pre.star[((u+1):n),] > std.tc[ii])*1}  # build uncont. array
                pre.tc <- as.data.frame(rbind(cnt.pre.tc, std.pre.tc))
              }
              else {  # if no contamination
                std.pre.tc <- array(-2, dim=c(n, k)) 
                for (ii in 1:(b-1)) {std.pre.tc <- std.pre.tc + (pre.star[1:n,] > std.tc[ii])*1}  # then all n respondents should be calculated with the no-RS tau cuts
                pre.tc <- as.data.frame(std.pre.tc) 
              }
              colnames(pre.tc) <- paste('x', 1:k, sep = '')
              
              # treatment and post test:
              trt <- t 
              x.trt <- sample(x = c(1, -1), size = n, prob= c(0.5, 0.5), replace = T)  # assign respondents to treatment / control
              xs <- data.frame(trt = x.trt)  # treatments in data frame form to pass to lavaan
              delta <- rnorm(mean = 0, sd = (0.66)^(1/2), n = n)
              post <- pre*p + x.trt*trt + delta
              eps <- array(rnorm(mean = 0, sd = (1-lambda^2)^(1/2), n = n*k), dim = c(k, n)) 
              post.star <- t(lambda %*% t(post)) + t(eps)
              post.itms <- as.data.frame(post.star)  # posttest uncategorized answers
              colnames(post.itms) <- paste('y', 1:k, sep = '')
              if (j != 0) {
                cnt.post.tc <- array(-2, dim=c(u, k))  # array for contaminated answers
                std.post.tc <- array(-2, dim=c(v, k))  # array for uncontaminated answers
                for (ii in 1:(b-1)) {cnt.post.tc <- cnt.post.tc + (post.star[(1:u),] > cnt.tc[ii])*1}  # build cont. array
                for (ii in 1:(b-1)) {std.post.tc <- std.post.tc + (post.star[((u+1):n),] > std.tc[ii])*1}  # build uncont. array
                post.tc <- as.data.frame(rbind(cnt.post.tc, std.post.tc))
              }
              else {
                std.post.tc <- array(-2, dim=c(n, k))
                for (ii in 1:(b-1)) {std.post.tc <- std.post.tc + (post.star[1:n,] > std.tc[ii])*1}
                post.tc <- as.data.frame(std.post.tc)
              }
              colnames(post.tc) <- paste('y', 1:k, sep = '')
              
              # input for lavaan to process
              input.df <- cbind(post.tc, xs, pre.tc) 
              
              # lavaan fit
              fit <- lavaan::sem(lmod, input.df)
              tryCatch(
                {
                  co <- standardizedsolution(fit)[which(standardizedsolution(fit)$lhs == 'post' & 
                                                          standardizedsolution(fit)$rhs == 'trt' | 
                                                          standardizedsolution(fit)$lhs == 'post' & 
                                                          standardizedsolution(fit)$rhs == 'pre'),] 
                  post.trt <- co$est[1]
                  post.trt.se <- co$se[1]
                  post.pre <- co$est[2]
                  post.pre.se <- co$se[2]
                  
                  chisq <- fitmeasures(fit, 'chisq')
                  df <- fitmeasures(fit, 'df')
                  cfi <- fitmeasures(fit, 'cfi')
                },
                error = function(e) {
                  post.trt <- NaN
                  post.trt.se <- NaN
                  post.pre <- NaN
                  post.pre.se <- NaN
                  chisq <- NaN
                  df <- NaN
                  cfi <- NaN
                  print(paste("error: ", e))
                }
              )
              
              # add estimates to data frame
              row <- array(c(
                as.numeric(k),  # manif
                as.numeric(b),  # gran
                as.numeric(g),  # shift
                as.numeric(n),  # n
                as.numeric(j),  # pct
                as.numeric(chisq),
                as.numeric(df),
                as.numeric(cfi),
                as.numeric(post.trt),
                as.numeric(post.pre),
                as.numeric(post.trt.se),
                as.numeric(post.pre.se))) 
              
              est[count,] <- row
            }  
          }
        }
      }
    }
  }
  
  est <- as.data.frame(est)
  colnames(est) <- c('manif', 
                     'gran',
                     'shift',
                     'n',
                     'pct',
                     'chisq', 
                     'df',
                     'cfi',
                     'post.trt',
                     'post.pre',
                     'post.trt.se',
                     'post.pre.se') 
  
  return(est)
}

mdf5 <- par.list[[1]]
for(ctr.i in 2:num_cores) {
  mdf5 <- rbind(mdf5, par.list[[ctr.i]])
}

write.csv(mdf5, 'mdf.05.18.csv')
