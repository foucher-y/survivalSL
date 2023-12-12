
auc <- function(sens, spec)
{
.tab.res <- data.frame(se=sens, sp=spec)
.tab.res <- .tab.res[!is.na(.tab.res$sp + .tab.res$se),]
.tab.res$sp1 <- 1-.tab.res$sp
.tab.res <- .tab.res[order(.tab.res$sp1, .tab.res$se),]
.tab.res <- rbind(c(0,1,0), .tab.res, c(1,0,1))

return( sum((.tab.res$sp1[2:length(.tab.res$sp1)] -
  .tab.res$sp1[1:(length(.tab.res$sp1)-1)]) * 0.5 * 
  (.tab.res$se[2:length(.tab.res$se)]+.tab.res$se[1:length(.tab.res$se)-1])) )
}


