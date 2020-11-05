library(doParallel)
library(foreach)
progcombine <- function(nreps) {
  defcombine <- function(a, ...) c(a, list(...))
  pb <- txtProgressBar(min=1, max=nreps-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    flush.console()
    defcombine(...)
  }
}

progcombine_c <- function(nreps) {
  pb <- txtProgressBar(min=1, max=nreps-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    flush.console()
    c(...)
  }
}

progcombine_rbind <- function(nreps) {
  pb <- txtProgressBar(min=1, max=nreps-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    flush.console()
    rbind(...)
  }
}

progcombine_cbind <- function(nreps) {
  pb <- txtProgressBar(min=1, max=nreps-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    flush.console()
    cbind(...)
  }
}