#--------------------------------------
# prediction of the DNA state of MGMT promoter
#  based on HM infinium platform
# 2014-09-09 modified by pbady
# 2015-05-07 modified by pbady
#--------------------------------------
data(MGMTSTP27, envir=environment())
.starsig <- function(x,cutpoints=NULL,symbols=NULL,...){
	if(is.null(cutpoints))
		cutpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
	if(is.null(symbols))
		symbols <- c("***", "**", "*", ".", " ")
	cut(x,breaks=cutpoints,label=symbols)	
}	
.roundsig <- function(x,digits=3,...){
	ifelse(x < 10^(-digits),paste("<0.",paste(rep(0,digits-1),collapse=""),"1",sep=""),round(x,digits))
}
.qqlinexy <- function (x,y,rev=FALSE,...){
	y <- quantile(y[!is.na(y)], c(0.25, 0.75))
	x <- quantile(x[!is.na(x)], c(0.25, 0.75))
	if (rev) {
		slope <- diff(x)/diff(y)
		int <- x[1] - slope * y[1]
	}
	else {
		slope <- diff(y)/diff(x)
		int <- y[1] - slope * x[1]
	}
	abline(int, slope, ...)
}
.logit <- function(x) log(x/(1-x))
.invlogit <- function(x) 1/(1+exp(-x))
# .sclass function based on the function 's.class' from ade4
.sclass <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
    cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
    clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
    xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 
        0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
    cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE,gbox=TRUE){
    f1 <- function(cl) {
        n <- length(cl)
        cl <- as.factor(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (unclass(cl) - 1)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        data.frame(x)
    }
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    dfxy <- data.frame(dfxy)
    if (!is.data.frame(dfxy)) 
        stop("Non convenient selection for dfxy")
    if (any(is.na(dfxy))) 
        stop("NA non implemented")
    if (!is.factor(fac)) 
        stop("factor expected for fac")
    dfdistri <- f1(fac) * wt
    coul = col
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
    cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
    if (nrow(dfxy) != nrow(dfdistri)) 
        stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (cpoint > 0) 
        for (i in 1:ncol(dfdistri)) {
            pch <- rep(pch, length = nrow(dfxy))
            points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, 
                i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
                cpoint, col = coul[i])
        }
    if (cstar > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                coul[i])
        }
    if (cellipse > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                cellipse = cellipse, axesell = axesell, coul[i])
        }
    if (clabel > 0) 
        scatterutil.eti(coox, cooy, label, clabel, coul = col)
    if(gbox)
        box()
    invisible(match.call())
}
# IC => normal, see Faraway (2006)
MGMTpredict <- function (x, level = 0.05, dispersion = FALSE, transpose = FALSE,ic.distrib="normal",cutoff=1,...){
  if (!inherits(x, "data.frame")) 
      stop("non convient object!")
  if(cutoff==1){
		  perf <- MGMTSTP27$perf1
	}else if(cutoff==2){
		  perf <- MGMTSTP27$perf2	
	}else{
	  stop("non convenient cut-off!")
	}
	if (!transpose) {
        if (!is.element("cg12981137", colnames(x))) 
            stop("the probe 'cg12981137' is missing!")
        if (!is.element("cg12434587", colnames(x))) 
            stop("the probe 'cg12434587' is missing!")
        data1 <- as.data.frame(x[, c("cg12981137", "cg12434587")])
    }
    else {
        if (!is.element("cg12981137", rownames(x))) 
            stop("the probe 'cg12981137' is missing!")
        if (!is.element("cg12434587", rownames(x))) 
            stop("the probe 'cg12434587' is missing!")
        data1 <- as.data.frame(t(x[c("cg12981137", "cg12434587"), 
            ]))
    }
    colnames(data1) <- c("cg12981137", "cg12434587")
    predmod <- predict(MGMTSTP27, newdata = data1, type = "link", 
        se.fit = TRUE)
    linkinv <- MGMTSTP27$family$linkinv
    df <- MGMTSTP27$df.residual
    pred <- predmod$fit
    if (dispersion) {
        sigma <- sum((MGMTSTP27$weights * MGMTSTP27$residuals^2)[MGMTSTP27$weights > 
            0])/df
    }
    else {
        sigma <- 1
    }
    sy <- predmod$se.fit * sqrt(sigma)

	if(ic.distrib=="student"){		
		lower <- pred - sy * qt(1 - level/2, df = df)
		upper <- pred + sy * qt(1 - level/2, df = df)
    }else if(ic.distrib=="normal"){
		lower <- pred - sy * qnorm(1-level/2)
		upper <- pred + sy * qnorm(1-level/2)
	}else stop("non convenient ic.distrib!")
	lower <- linkinv(lower)
    upper <- linkinv(upper)
    pred <- linkinv(pred)
    categ <- ifelse(pred >= perf$cut, "M", "U")
    if (is.null(rownames(data1))) 
        rownames(data1) <- 1:nrow(data1)
    res <- data.frame(sample = rownames(data1), cg12434587 = data1[, 
        "cg12434587"], cg12981137 = data1[, "cg12981137"], pred = pred, 
        lower = lower, upper = upper, state = categ)
    class(res) <- c("mgmt", class(res))
    return(res)
}
MGMTsim <- function(n=1000,newdata=NULL,maxIteration=100,tol=1e-3,method="each",expd.grid=FALSE,...){
	if(method=="each"){
  if(is.null(newdata)){
		gm1 <- gammaFitEM(MGMTSTP27$data[,c("cg12981137")],maxIteration=maxIteration,tol=tol,...)
		gm2 <- gammaFitEM(MGMTSTP27$data[,c("cg12434587")],maxIteration=maxIteration,tol=tol,...)
	}else{
		if(!all(is.element(c("cg12981137","cg12434587"),colnames(newdata))))
		        stop("the probes 'cg12981137'and/or 'cg12434587' are missing!")
		gm1 <- gammaFitEM(newdata[,c("cg12981137")],maxIteration=maxIteration,tol=tol,...)
		gm2 <- gammaFitEM(newdata[,c("cg12434587")],maxIteration=maxIteration,tol=tol,...)
	}
	simdata <- as.data.frame(cbind(.simprobe(n,gm1),.simprobe(n,gm2)))
  }else if(method=="all"){
    if(is.null(newdata)){
      gm3 <- gammaFitEM(unlist(MGMTSTP27$data),maxIteration=maxIteration,tol=tol,...)
    }else{
      if(!all(is.element(c("cg12981137","cg12434587"),colnames(newdata))))
        stop("the probes 'cg12981137'and/or 'cg12434587' are missing!")
      gm1 <- gammaFitEM(unlist(newdata),maxIteration=maxIteration,tol=tol,...)
    }
    simdata <- as.data.frame(cbind(.simprobe(n,gm3),.simprobe(n,gm3)))    
  }
  if(expd.grid)
		simdata <- expand.grid(simdata[,1],simdata[,2])
	colnames(simdata) <- c("cg12434587","cg12981137")
	return(simdata)
}
.simprobe <- function(n=1000, gammaFit = NULL, k = NULL, 
                      theta = NULL, shift = NULL,proportion = NULL, ...){
# require(lumi)
   if (!is.null(gammaFit)) {
        if (class(gammaFit) != "gammaFit") 
            stop("gammaFit should be an object of 'gammaFit' class!")
        k <- gammaFit$k
        theta <- gammaFit$theta
        shift <- gammaFit$shift
        proportion <- gammaFit$proportion
    }
    if (is.null(k) || is.null(theta) || is.null(shift) || is.null(proportion)) 
        stop("Information of parameters k, theta, shift and proportion is required!\n")
    n1 <- rbinom(n = 1, size = n, prob = proportion[1])
	y1 = shift[1] + rgamma(n1, shape = k[1], scale = theta[1])
    y2 =  shift[2] - rgamma(n-n1, shape = k[2], scale = theta[2])
	return(c(y1,y2))
}
MGMTqc.pop <- function(object,sim=FALSE,n=nrow(object),cutoff=1,which.plot=1:6,
                       mfrow=c(3,3),...){
    if (!inherits(object, "mgmt")) 
        stop("object 'mgmt' expected !")
	if(is.element("cg12981137",rownames(object)))
		stop("the probe 'cg12981137' is missing!")
	if(is.element("cg12434587",rownames(object)))
		stop("the probe 'cg12434587' is missing!")
	if(nrow(object) < 2)
		stop("non convenient dimension!")
	if(cutoff==1){
		perf <- MGMTSTP27$perf1
	}else if(cutoff==2){
		perf <- MGMTSTP27$perf2	
	}else{
		stop("non convenient cut-off!")
	}
	data1 <- as.data.frame(object[,c("cg12981137","cg12434587")])
	state1 <- object[,"state"]
	if(sim){
		simdata <- MGMTsim(n=n,...)
		simpred <- MGMTpredict(simdata)
		simstate <- factor(as.character(simpred$state),levels=c("U","M"))
		simlab <- "simulated"
	}else{
		simdata <- MGMTSTP27$data[,c("cg12981137","cg12434587")]
		simstate <- factor(ifelse(MGMTSTP27$fitted.values >= perf$cut,"M","U"),levels=c("U","M"))
		simlab <- "training"
	}
	cut1 <- perf[1]
	if(!is.null(mfrow))
	  par(mfrow=mfrow,mar=c(4.1,4.1,2.1,2.1))
	# Global 2D
	coord1 <- .logit(object[,"pred"])
	ord1 <- order(coord1)
	coord1 <- coord1[ord1]
	pred <- object[ord1,"pred"]
	lwr <- object[ord1,"lower"]
	upr <- object[ord1,"upper"]
  if(any(which.plot%in%1)){
  	plot(simdata,type="n",xlab=colnames(simdata)[1],ylab=colnames(simdata)[2],panel.first=c(grid(),abline(v=0),abline(h=0)))
      xykde = MASS::kde2d(simdata[,1], simdata[,2], lims = par("usr"))
      zlim = range(xykde$z, finite = TRUE)
      lev = seq(zlim[1], zlim[2], le = 8)
      lev = lev[2:7]
      contour(xykde, add = TRUE, lwd = 1, col = "gray", levels = lev,drawlabels = FALSE,lty=2)
  	points(data1[,1],data1[,2],pch=18,col=ifelse(state1=="M","red","blue"))
  	.sclass(simdata,fac=simstate,cstar=0,add.p=TRUE,gbox=FALSE,col=c("blue","red"),cpoint=0,clabel=2)
    title("Two-Dimensional Representation",adj=0,font.main=1,cex.main=1)
	}
	# global hist 
	d1 <- density(coord1)
	d0 <- density(na.omit(.logit(MGMTSTP27$fitted.values)))
	if(any(which.plot%in%2)){
  	hist(coord1,main="",xlab="logit(proba)",nclass=13,col="lightgrey",proba=TRUE)
  	lines(d1,col="black",lwd=2)
  	lines(d0,col="green3",lwd=2)
  	abline(v=.logit(cut1),col="darkgreen",lty=2)	
  	legend("topright",c("fitted",simlab),lty=1,lwd=2,col=c("black","green3"),cex=0.9,box.lty=0)
  	box()
  	title("Histogram for the MGMT score",adj=0,font.main=1,cex.main=1)
	}
	# global predict
	if(any(which.plot%in%3)){
  	plot(coord1,pred,ylim=c(0,1),xlab="logit(proba)",ylab="Probability",panel.first=c(grid()),type="n")
  	polygon(c(coord1,rev(coord1)),c(upr,rev(lwr)),col="lightgrey",border=FALSE)
  	lines(coord1,pred,lwd=2,col="black")
  	points(coord1,pred,pch=19,col="black")
  	abline(h=cut1,col="darkgreen",lty=2)
  	abline(v=.logit(cut1),col="darkgreen",lty=2)
  	abline(h=cut1,col="darkgreen",lty=2)
  	axis(4,at=c(0,1),c("U","M"),las=1)
  	legend("bottomright",c("fitted","cut-off","IC"),lty=c(1,2,NA),
  		col=c("black","darkgreen","lightgrey"),bg="white",
  		pch=c(19,NA,15),pt.cex=c(1,NA,1.5),cex=0.9,box.lty=0)
  	box()
  	title("Expected values and tolerance intervals",adj=0,font.main=1,cex.main=1)
	}
	# KS test + qqplot by prob
	if(any(which.plot%in%c(4,5,6))){
    for(k in c("cg12434587","cg12981137")){
  		tmp0 <- simdata[,k]
  		tmp1 <- data1[,k]
  		if(any(which.plot%in%4)){
    		ks1 <- ks.test(tmp1,tmp0)
    		qqplot(tmp0,tmp1,panel.first=c(grid()),xlab=paste(simlab,"values"),ylab="observed values")
    		abline(0,1,col="red",lty=3)
    		.qqlinexy(tmp0,tmp1,col="red")
    		legend("topleft",k,box.lty=0)
    		legend("bottomright",paste("D=",round(ks1$statistic,3),"; p=",.roundsig(ks1$p.value,3),sep=""),cex=0.9,box.lty=0)
    		title(paste("QQ-plot for ",k,sep=""),adj=0,font.main=1,cex.main=1)
  		}
  	# density plot
  	  if(any(which.plot%in%5)){
    		d1 <- density(tmp1)
    		d0 <- density(tmp0)
    		hist(tmp1,xlab=k,main="",nclass=13,col="lightgrey",proba=TRUE)
    		lines(d1,col="black",lwd=2)
    		lines(d0,col="green3",lwd=2)
    		legend("topright",c("observed",simlab),lty=1,lwd=2,col=c("black","green3"),cex=0.9,box.lty=0)
    		box()		
        title(paste("Histogram for ",k,sep=""),adj=0,font.main=1,cex.main=1)
  	  }
  	# marginal effect?
  	  if(any(which.plot%in%6)){
    		ord1 <- order(object[,k])
    		coord1 <- object[ord1,k]
    		pred <- object[ord1,"pred"]
    		lwr <- object[ord1,"lower"]
    		upr <- object[ord1,"upper"]
    		plot(coord1,pred,ylim=c(0,1),xlab=k,ylab="Probability",panel.first=c(grid()),type="n")
    		segments(coord1,upr,coord1,lwr,col="black")
    		points(coord1,pred,pch=19,col="black")
    		lines(lowess(pred~coord1),col="red",lwd=2)
    		abline(h=cut1,col="darkgreen",lty=2)
    		axis(4,at=c(0,1),c("U","M"),las=1)
    		legend("bottomright",c("fitted","cut-off","loess"),lty=c(1,2,1),
    			col=c("black","darkgreen","red"),bg="white",
    			pch=c(19,NA,NA),pt.cex=c(1,NA,NA),cex=0.9,box.lty=0)
    		box()
    	  title(paste("Expected values in function of ",k,sep=""),adj=0,font.main=1,cex.main=1)
  	  }
	  }
	}
	invisible(match.call())
}
MGMTqc.single <- function(object,nsample=NULL,sim=FALSE,n=nrow(object),cutoff=1,which.plot=1:4,
                          mfrow=c(2,3),...){
    if (!inherits(object, "mgmt")) 
        stop("object 'mgmt' expected !")
	if(is.element("cg12981137",rownames(object)))
		stop("the probe 'cg12981137' is missing!")
	if(is.element("cg12434587",rownames(object)))
		stop("the probe 'cg12434587' is missing!")
	if(is.null(nsample))
		nsample <- 1
	object <- object[nsample,]
	if(cutoff==1){
		perf <- MGMTSTP27$perf1
	}else if(cutoff==2){
		perf <- MGMTSTP27$perf2	
	}else{
		stop("non convenient cut-off!")
	}
	data1 <- as.data.frame(object[,c("cg12981137","cg12434587")])
	state1 <- object[,"state"]
	if(sim){
		simdata <- MGMTsim(n=n,...)
		simpred <- MGMTpredict(simdata)
		simstate <- factor(as.character(simpred$state),levels=c("U","M"))
		simlab <- "simulated"
	}else{
		simdata <- MGMTSTP27$data[,c("cg12981137","cg12434587")]
		simpred <- MGMTpredict(simdata)
		simstate <- factor(as.character(simpred$state),levels=c("U","M"))
		simlab <- "training"
	}
	cut1 <- perf[1]
  if(!is.null(mfrow))
	  par(mfrow=mfrow,mar=c(4.1,4.1,2.1,2.1))
	# Global 2D
  if(any(which.plot%in%1)){
  	plot(simdata,type="n",xlab=colnames(simdata)[1],ylab=colnames(simdata)[2],panel.first=c(grid(),abline(v=0),abline(h=0)))
      xykde = MASS::kde2d(simdata[,1], simdata[,2], lims = par("usr"))
      zlim = range(xykde$z, finite = TRUE)
      lev = seq(zlim[1], zlim[2], le = 8)
      lev = lev[2:7]
      contour(xykde, add = TRUE, lwd = 1, col = "gray", levels = lev,drawlabels = FALSE,lty=2)
  	points(data1[,1],data1[,2],pch=18,cex=1.5,col=ifelse(state1=="M","red","blue"))
  	.sclass(simdata,fac=simstate,cstar=0,add.p=TRUE,gbox=FALSE,col=c("blue","red"),cpoint=0,clabel=2)
  	title("Two-Dimensional Representation",adj=0,font.main=1,cex.main=1)
  }
	# global hist 
	coord1 <- .logit(object[1,"pred"])
	pred <- object[1,"pred"]
	lwr <- object[1,"lower"]
	upr <- object[1,"upper"]
	d0 <- density(na.omit(.logit(MGMTSTP27$fitted.values)))
  y0 <- max(d0$y)
	if(any(which.plot%in%2)){
    plot(d0,main="",xlab="logit(proba)",col="green3",lwd=2)
  	lines(c(coord1, coord1), c(2*y0/3, 0),col="black")
      points(coord1, 2*y0/3, pch = 18, cex = 2,col="black")
  	abline(v=.logit(cut1),col="darkgreen",lty=2)
  	legend("topright",c("fitted",simlab),lty=1,lwd=2,col=c("black","green3"),cex=0.9,box.lty=0)
  	box()
  	title(paste("Density plot for MGMT score",sep=""),adj=0,font.main=1,cex.main=1)
	}
	# global predict
	simcoord1 <- .logit(simpred[,"pred"])
	ord1 <- order(simcoord1)
	simcoord1 <- simcoord1[ord1]
	simpred1 <- simpred[ord1,"pred"]
	simupr1 <- simpred[ord1,"upper"]
	simlwr1 <- simpred[ord1,"lower"]
	if(any(which.plot%in%3)){	
  	plot(simcoord1,simpred1,ylim=c(0,1),xlab="logit(proba)",ylab="Probability",type="n",panel.first=c(grid()))
  	polygon(c(simcoord1,rev(simcoord1)),c(simupr1,rev(simlwr1)),col="lightgrey",border=FALSE)
    #polygon(c(.logit(upr),coord1,.logit(lwr),coord1),c(pred,upr,pred,lwr),col=grey(0.7))
  	lines(simcoord1,simpred1,lwd=2,col="green3")
  	segments(coord1,upr,coord1,lwr,col="black")
  	segments(.logit(upr),pred,.logit(lwr),pred,col="black")
    points(coord1,pred,pch=19,cex=1.5,col="black")
  	abline(h=cut1,col="darkgreen",lty=2)
  	abline(v=.logit(cut1),col="darkgreen",lty=2)
  	abline(h=cut1,col="darkgreen",lty=2)  
  	axis(4,at=c(0,1),c("U","M"),las=1)
  	legend("bottomright",c("fitted","training","cut-off","IC"),lty=c(1,1,2,NA),
  		col=c("black","green3","darkgreen","lightgrey"),bg="white",
  		pch=c(19,NA,NA,15),pt.cex=c(1.5,NA,NA,1.5),cex=0.9,box.lty=0)
    box()
  	title("Expected values and tolerance intervals",adj=0,font.main=1,cex.main=1)
	}
	# Density plot
	if(any(which.plot%in%4)){
  	for(k in c("cg12434587","cg12981137")){
  		tmp0 <- simdata[,k]
  		coord1 <- data1[1,k]
  	# density plot
  		d0 <- d0 <- density(tmp0)
  		y0 <- max(d0$y)
  		plot(d0,main="",xlab=k,col="green3",lwd=2)
  		lines(c(coord1, coord1), c(2*y0/3, 0),col="black")
  		points(coord1, 2*y0/3, pch = 18, cex = 2,col="black")
  		legend("topright",c("observed",simlab),lty=1,lwd=2,col=c("black","green3"),cex=0.9,box.lty=0)
  		box()		
  	  title(paste("Density plot for ",k,sep=""),adj=0,font.main=1,cex.main=1)
  	}
	}
	invisible(match.call())
}