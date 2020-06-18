plot.misfit <- function(x, se = TRUE){

  beta.hat <- x[["beta.hat"]]
  M <- length(x[["beta.hat"]])
  grid <- seq(from = 0,to = 1,length.out = M)
  plot(grid, beta.hat, type = 'l')

  if(se){
    if(is.list(x[["Cbeta"]])){
      Cbeta <- x[["Cbeta"]]["var.t"]
    }else{
      Cbeta <- x[["Cbeta"]]
    }
    me <- 1.96*sqrt(diag(Cbeta))
    lines(grid, beta.hat - me, lty = 2)
    lines(grid, beta.hat + me, lty = 2)
  }
}
