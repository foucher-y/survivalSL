
plot.survrisca <- function(x, ..., col=1, lty=1, lwd=1, type="s", max.time=NULL, min.y=0, max.y=1, grid.lty=NULL){

col <- rep(col, ceiling(length(sort(unique(x$table.surv$variable)))/length(col)) )
lty <- rep(lty, ceiling(length(sort(unique(x$table.surv$variable)))/length(lty)) )
lwd <- rep(lwd, ceiling(length(sort(unique(x$table.surv$variable)))/length(lwd)) )

max.x <- min(max.time, max(x$table.surv$times))

plot.default(NULL, xlim=c(0, max.x), ylim=c(min.y, max.y), col=col, lty=lty, lwd=lwd, type=type, ...)

if(!is.null(grid.lty)) {grid(lty=grid.lty, col="grey90")}
        
j <- 1
	for (i in sort(unique(x$table.surv$variable))){
			lines(
			x$table.surv$times[x$table.surv$variable==i & x$table.surv$times<=max.x],
			x$table.surv$survival[x$table.surv$variable==i & x$table.surv$times<=max.x],
			col=col[j], lty=lty[j], lwd=lwd[j], type=type
			)
		j <- j + 1
	}
}
# 
# data(dataDIVAT2)
# 
# # adjusted Kaplan-Meier estimator by IPW
# Pr0 <- glm(ecd ~ 1, family = binomial(link="logit"), data=dataDIVAT2)$fitted.values[1]
# Pr1 <- glm(ecd ~ age + hla + retransplant, data=dataDIVAT2,
#            family=binomial(link = "logit"))$fitted.values
# W <- (dataDIVAT2$ecd==1) * (1/Pr1) + (dataDIVAT2$ecd==0) * (1)/(1-Pr1)
# res.akm <-ipw.survival(times=dataDIVAT2$times, failures=dataDIVAT2$failures,
#                        variable=dataDIVAT2$ecd, weights=W)
# 
# plot(res.akm, ylab="Confounder-adjusted survival",
#      xlab="Time post-transplantation (years)", col=c(1,2), grid.lty=1)
#  
# 




# 
# data(dataDIVAT2)
# 
# res <- survfit(Surv(times,failures) ~ ecd, data=dataDIVAT2)
# 
# couleurs <- c("royalblue", "#F8766D", "#00BA38")
# 
# # il faut commencer par un plot vide pour tracer le quadrillage avant les courbes sinon le quadrillage se retrouve dessus
# plot(NULL, xlim=c(0, max(res$time)), ylim=c(0,1), ylab="Confounder-adjusted survival", xlab="Time post-transplantation (years)")
# 
# grid(lty=1, col="grey90")
# 
# lines(res, col=couleurs)
# 
# for(i in 1:length(res$strata)){
#         x <- c(res[i]$time, rev(res[i]$time))
#         y <- c(res[i]$upper, rev(res[i]$lower))
#         polygon(x, y, col = adjustcolor(couleurs[i], alpha.f=0.4), border=FALSE)
# }
#  
#  
#  
 
 
 
 
# lines(res.km$times[res.km$variable==1], res.km$survival[res.km$variable==1],
# type="s",col=2,lty=2,lwd=2)
# lines(res.km$times[res.km$variable==0], res.km$survival[res.km$variable==0],
# type="s",col=1,lty=2,lwd=2)

# adjusted Kaplan-Meier estimator by IPW
#  Pr0 <- glm(ecd ~ 1, family = binomial(link="logit"), data=dataDIVAT2)$fitted.values[1]
#  Pr1 <- glm(ecd ~ age + hla + retransplant, data=dataDIVAT2,
#  family=binomial(link = "logit"))$fitted.values
#  W <- (dataDIVAT2$ecd==1) * (1/Pr1) + (dataDIVAT2$ecd==0) * (1)/(1-Pr1)
#  res.akm <-ipw.kaplan.meier(times=dataDIVAT2$times, failures=dataDIVAT2$failures,
#   variable=dataDIVAT2$ecd, weights=W)
# lines(res.akm$times[res.akm$variable==1], res.akm$survival[res.akm$variable==1],
#  type="s",col=2,lwd=2)
# lines(res.akm$times[res.akm$variable==0], res.akm$survival[res.akm$variable==0],
#  type="s",col=1,lwd=2)

# nb.risk1<-function(x) {sum(dataDIVAT2$times[dataDIVAT2$ecd==0]>x)}
# nb.risk2<-function(x) {sum(dataDIVAT2$times[dataDIVAT2$ecd==1]>x)}
# segments(x0=0, y0=0.1, x1=13, y1=0.1) 
# text(x=6, y=0.12, "number of at-risk patients", cex=0.8)
# tps <- seq(1,12,by=1)
# text(x=tps, y=rep(0.07,length(tps)), as.character(sapply(tps, FUN="nb.risk1")),
#  cex=0.8, col=1)
# text(x=tps, y=rep(0.02,length(tps)), as.character(sapply(tps, FUN="nb.risk2")),
#  cex=0.8, col=2)
# legend("topright", legend=c("Unadjusted estimator for SCD",
#  "Adjusted estimator for SCD", "Unadjusted estimator for ECD",
#  "Adjusted estimator for ECD"), col=c(1,1,2,2),
#  lty=c(2,1,2,1), lwd=2, cex=0.8)


