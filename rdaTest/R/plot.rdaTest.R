'plot.rdaTest' <- 
	function(x, xax=1, yax=2, scaling=1, graph.type="notchosen",  
	plot.sites=TRUE, plot.spe=TRUE, plot.env=TRUE, binary=NULL, height=6, 
	width=6, xlim=NULL, ylim=NULL, ell=NULL, alpha=0.95, ell.axis=FALSE,  
	mul.spe=0.8, mul.env=0.90, mul.text=0.10, pos.ell=NULL, pos.bin=NULL,  
	pos.sites=1, label.sites=TRUE, label.spe=TRUE,label.env=TRUE, 
	label.ell=TRUE, col.env="blue", col.spe="red", col.site="black",
	col.ell="black", lty.env=1, lty.spe=1, lty.ell=1, lty.axis=2, cex=1, 
	cex.lab=1, cex.axis=1, lwd=1, len=0.1, saveplot=FALSE, path=NULL, 
	mar.perc=NULL, interior.mar.perc=0.05, select.spe=NULL, ...)
#
# License: GPL-2 
# Authors:
# Version 1.0: SŽbastien Durand Pierre Legendre, UniversitŽ de MontrŽal, 2005
# Version 1.2: Pierre Legendre, UniversitŽ de MontrŽal, April 2008
# Version 1.3: Pierre Legendre, UniversitŽ de MontrŽal, October 2010
{
    #### Internal function
	larger.frame <- function(mat, percent=0.05, xinv, yinv)
	# Produce an object plot 10% larger than strictly necessary
	{
	if(xinv<0) mat[,1] <- -mat[,1]
	if(yinv<0) mat[,2] <- -mat[,2]
	range.mat = apply(mat,2,range)
	z <- apply(range.mat, 2, function(x) x[2]-x[1])
	range.mat[1,]=range.mat[1,]-z*percent
	range.mat[2,]=range.mat[2,]+z*percent
	range.mat
	}
	####

# Save plot
	if(saveplot==TRUE){
		if(is.null(path)){
			cat("Type a file name describing this plot: ex.: toto","\n")
			nom<-scan(file="",what="character",quiet=TRUE,nline=1)
			cat("Do you want to save '",nom,
				".pdf' in the working directory: YES='Y' or NO='N'","\n",sep="")
			tmp<-scan(file="",what="character",quiet=TRUE,nline=1)
			if(tmp!="Y"){
				cat("Type the directory path you require ~/.../:","\n")
				tmp1<-scan(file="",what="character",quiet=TRUE,nline=1)
				path<-paste(tmp1,nom,".pdf",sep="")
			}else{
				path<-paste(getwd(),"/",nom,".pdf",sep="")
			}
		}
		pdf(file=path, onefile=TRUE, paper="special", width=width, 
		    height=height, family="Helvetica")
	} else {
		if(length(grep("darwin",R.Version()$os))!=1){
			windows(w=width,h=height)
		}else{
			quartz(w=width,h=height)
		}
	}	
	old.par <- par(no.readonly=TRUE) # all par settings that can be changed.
 	on.exit(par(old.par))
	
# Axis reversion filter
	xinv <- 1
	if(xax < 0) {
		xinv <- (-1)
		xax  <- -(xax)
	}
	yinv <- 1
	if(yax < 0) {
		yinv <- (-1)
		yax  <- -(yax)
	}
	k <- length(x$eig.values)
	if(colnames(x$U.sc1)[2] == "PCA_Axis1") k <- 1
	# cat("k =",k,"\n")
	if(xax > length(x$eig.values) | xax > k)
		stop("Their are not enough canonical axes; change the xax value")
	if(yax > length(x$eig.values))
		stop("Their are not enough canonical axes; change the yax value")

# Scaling 	
	if(scaling != 1) {
		x$F.sc1 <- x$F.sc2
		x$Z.sc1 <- x$Z.sc2
		x$U.sc1 <- x$U.sc2
		x$biplotScores1 <- x$biplotScores2
	}

	graph.type<-toupper(as.character(graph.type))
	if(graph.type=="notchosen" | graph.type!="F" & graph.type!="Z"){
		while(graph.type!="F" & graph.type!="Z") {
		cat("Choose a functional graph.type: type 'Z' for fitted site scores",
			" or 'F' for site scores","\n") 
		graph.type <- toupper(scan(file="", what="character", nlines=1,quiet=TRUE))
		}
	}
	
	if(graph.type != "Z") {
		x$Z.sc1 <- x$F.sc1
	}
	
	if(!is.null(binary)) {
		binvar <- matrix(NA,length(binary),ncol(x$Z.sc1))
		j=1
		for(i in binary) {
			tmp <- which(x$X.mat[,i]!=0)
			if(length(tmp) == 1) { 
			   binvar[j,] <- x$Z.sc1[tmp,] 
			   } else { 
			   binvar[j,] <- apply(x$Z.sc1[tmp,]*x$X.mat[tmp,i],2,mean)
			   }
			j=j+1
		}
		colnames(binvar)<-colnames(x$Z.sc1)
		rownames(binvar)<-colnames(x$X.mat)[binary]
		# rownames(binvar)<-c("Ra","Ri","Rs","Lv","Ln","Ll","Ls","Lf","Sc", 
		# 	"Iba","Vl","Vp","Vh","Vc","Tw","Tb","Hd","At","As","Pb","Li","Cl",
		#	"Ch","Cb")
		# rownames(binvar)<-c("a","b","c","d","e","f","g","h","i","j","k","l", 
		# 	"m","n","o","p","q","r","s","t","u","v","w","x")
	}
	# Text scaling
	mul.text <- mul.text+1

# Draw the triplot
	
	#Control marging size
	if(length(mar.perc)!=0) par(mai=c(height*mar.perc,width*mar.perc,
		height*mar.perc/2,width*mar.perc/2))

	#Interior margin of the plot (space between outer points and frame)
	if(is.null(xlim) & is.null(ylim)) {
	lf.Z=larger.frame(x$Z.sc1[,c(xax,yax)],percent=interior.mar.perc, xinv,yinv)
		xlim <- lf.Z[,1]
		ylim <- lf.Z[,2]
		}

    # Plot site points and labels
    if(plot.sites) { type="p" }
       else {
       type="n"
       label.sites=FALSE
       }
    if(colnames(x$U.sc1)[2] == "PCA_Axis1") {
    # cat("If TRUE, colnames(x$U.sc1)[2] =", colnames(x$U.sc1)[2],"\n")
	plot(xinv*x$Z.sc1[,xax], yinv*x$Z.sc1[,yax], xlim=xlim, ylim=ylim, asp=1,
	    cex=cex, cex.lab=cex.lab, cex.axis=cex.axis,
		xlab=paste("Canonical axis",xax), ylab=paste("Residual PCA axis 1"), 
		type=type)
	} else {
    # cat("If FALSE, colnames(x$U.sc1)[2] =", colnames(x$U.sc1)[2],"\n")
	plot(xinv*x$Z.sc1[,xax], yinv*x$Z.sc1[,yax], xlim=xlim, ylim=ylim, asp=1,
	    cex=cex, cex.lab=cex.lab, cex.axis=cex.axis,
		xlab=paste("Canonical axis",xax), ylab=paste("Canonical axis",yax), 
		type=type)
	}
	if(label.sites) {
		text(xinv*x$Z.sc1[,xax], yinv*x$Z.sc1[,yax], rownames(x$Z.sc1), 
		col=col.site, cex=cex, pos=pos.sites)
		}

	# Draw abscissa and ordinate at origin, point (0,0)
	abline(h = 0, lty = lty.axis, lwd=lwd)
    abline(v = 0, lty = lty.axis, lwd=lwd)
	
	mult.spe <- multi.factor(x$Z.sc1, x$U.sc1, xax=xax, yax=yax, xinv=xinv,
		yinv=yinv, percentage=mul.spe)
	mult.env <- multi.factor(x$Z.sc1, x$biplotScores1, xax=xax, yax=yax,
		xinv=xinv, yinv=yinv, percentage=mul.env)

	# Draw species (or other response variable) arrows
	if(plot.spe) {
	if(is.null(select.spe)){ vec <- 1:nrow(x$U.sc1) } else { vec <- select.spe }
	arrows(x0=0, y0=0, x1=xinv*x$U.sc1[vec,xax]*mult.spe, 
		y1=yinv*x$U.sc1[vec,yax]*mult.spe, col=col.spe, code=2, lty=lty.spe, 
		len=len, lwd=1)
	if(label.spe==TRUE) text(x=mul.text*mult.spe*xinv*x$U.sc1[vec,xax],
		y=mul.text*mult.spe*yinv*x$U.sc1[vec,yax], rownames(x$U.sc1[vec,]), 
		col=col.spe, cex=cex)
	}

    if(plot.env) {
	# Draw environmental variable points if there are binary variables	
	if(!is.null(binary)) {
		points(x=xinv*binvar[,xax],y=yinv*binvar[,yax], pch=10,col=col.env)
		if(label.env==TRUE)
			text(x=xinv*binvar[,xax],y=yinv*binvar[,yax], 
			label=rownames(binvar), col=col.env,pos=pos.bin,cex=cex)
 		remain<-seq(1, ncol(x$X.mat),by=1)[-binary]
	} else {
		remain<-seq(1, ncol(x$X.mat),by=1)	
	}
	if(length(remain)!=0) { # Draw environmental variable arrows
		arrows(x0=0,y0=0,xinv*x$biplotScores1[remain,xax]*mult.env,
			yinv*x$biplotScores1[remain,yax]*mult.env, lty=lty.env,col=col.env, 
			code=2, len=len,lwd=1)
		if(label.env==TRUE)
			text(x=mul.text*mult.env*xinv*x$biplotScores1[remain,xax], 
				y=mul.text*mult.env*yinv*x$biplotScores1[remain,yax], 
				rownames(x$biplotScores1)[remain], col = col.env,cex=cex)
    }	

	# Draw additional axes on top and right of graph for environmental variables
	len=9
	axis(3, at = seq(-mult.env,mult.env,length=len), 
		labels=round(seq(-1,1,length=len),digit=2), 
		col = col.env,cex.axis=cex.axis)
   	axis(4, at = seq(-mult.env,mult.env,length=len), 
   		labels=round(seq(-1,1,length=len),digit=2), 
   		col=col.env,cex.axis=cex.axis)

	if(!is.null(ell)) {
		ell<-as.factor(ell)
		name.ell=levels(as.factor(ell))
		for(i in name.ell) {
			numb<-which(ell==i)
			if(label.ell==TRUE) {
			ellipse(XY.mat=cbind(xinv*x$Z.sc1[numb,xax],yinv*x$Z.sc1[numb,yax]),
				alpha=alpha,ell.type=lty.ell,ell.lwd=lwd, p.axis=ell.axis,
				name.pos=pos.ell,col=col.ell,name=i)
			} else {
			ellipse(XY.mat=cbind(xinv*x$Z.sc1[numb,xax],yinv*x$Z.sc1[numb,yax]),
				alpha=alpha,ell.type=lty.ell,ell.lwd=lwd, p.axis=ell.axis,
				name.pos=pos.ell,col=col.ell,name=NULL)
			}
		}
	}
	}
	
	# Close the PDF plot file
	if(saveplot==TRUE) dev.off()
	invisible()
}

multi.factor<-
	function(plot.m, arrow.m,xax=1,yax=2,xinv=1,yinv=1,percentage=0.80)
{
# Scale to the same proportions two drawings, such as plot and arrow, to make 
#    sure all arrows will be seen
# Get the matrix column ranges
	aa<-xinv*range(plot.m[,xax])
	bb<-yinv*range(plot.m[,yax])
	cc<-xinv*range(arrow.m[,xax])
	dd<-yinv*range(arrow.m[,yax])
# Find the arrow requiring more correction if exceeding the plot size, or
# find the arrow closest to the plot border to extend to span of the arrow set
	dif<-list()
	if(aa[1]<0 & cc[1]<0 | aa[1]>0 & cc[1]>0)
		dif[1]<-(aa[1]/cc[1])#/aa[1]
	else
		dif[1]<-NA
	if(aa[2]<0 & cc[2]<0 | aa[2]>0 & cc[2]>0)
		dif[2]<-(aa[2]/cc[2])#/aa[2]
	else
		dif[2]<-NA
	if(bb[1]<0 & dd[1]<0 | bb[1]>0 & dd[1]>0)
		dif[3]<-(bb[1]/dd[1])#/bb[1]
	else
		dif[3]<-NA
	if(bb[2]<0 & dd[2]<0 | bb[2]>0 & dd[2]>0)
		dif[4]<-(bb[2]/dd[2])#/bb[2]
	else
		dif[4]<-NA
		
	sel<-which.min(dif)
	border1<-c(aa[1],aa[2],bb[1],bb[2])
	border2<-c(cc[1],cc[2],dd[1],dd[2])
	mult<-((border1[sel]*percentage)/border2[sel])#/percentage
	return(mult)
}

ellipse<-
	function(XY.mat, alpha=0.95, ell.res=5, ell.lwd=1, ell.type=3,
	p.axis=TRUE, p.centroid=FALSE, centroid.type=20, col="black", 
	name=NULL, name.pos=NULL, lwd=1) 
{
	#
	# Draw confidence ellipses around groups of points.
	# Reference: Biometry (Sokal & Rohlf 1995), Section 15.7
	#
	# Confidence regions: values between 0.01 and 0.999
	# ell.res: the default is to produce a point each 5 degrees
	# p.axis: print the ellipse major and minor axes
	# p.centroid: print the centroid
	# ell.type:  type of line used to draw the ellipse and its axis;
	#    see 'lty' in the help documentation file of 'par'
	# centroid.type: type of centroid used to be displayed;
	#    see 'pch' in the help documentation file of 'par'
	# col: the color used for the ellipses;
	#    see 'col' in the help documentation file of 'par'
	# names: names of the ellipse to be plot (near the centroid)
	# name.pos: the position relative to the centroid pos: values = '1','2','3'
    #    and '4', respectively, indicate positions below, to the left of, above,  
    #    or to the right of the specified coordinates.
    #
	if(!is.logical(c(p.axis,p.centroid))) 
	cat("Choose either 'FALSE' or 'TRUE' for both 'p.axis' and 'p.centroid':\n")
	step=5
	# Compute the 95% confidence ellipse 95% if alpha = 0.95
	n <- length(XY.mat[,1])
	# Variance
	Sx <- var(XY.mat[,1]);  Sy <- var(XY.mat[,2])
	# Mean
	Mx <- mean(XY.mat[,1]); My <- mean(XY.mat[,2])
	Eig <- eigen(cov(cbind(XY.mat)))$values
	
	# Covariance
	COVxy <- cov(XY.mat)
	Lx <- Eig[1]
	if(Lx<0) Lx <- 0 
	
	Ly=Eig[2]
	if(Ly<0) Ly <- 0
	# Slopes of major and minor axes
	b1 <- COVxy[1,2]/(Lx-Sx)
	if(is.infinite(b1) | is.na(b1)) b1 <- 0
	b2 <- -1/b1
	if(is.infinite(b2) | is.na(b1)) b2 <- 0
	
	if(n<=2) {
		cat("Not enough individuals in group: ",name,"\n")
	} else {
		C <- ((Lx*Ly)*((n-1)*2)/(n-2)*(qf(alpha,2,n-2)))
		
		locc <- sqrt(C/(Ly*(1+b1^2)))
		if(is.infinite(locc) | is.na(locc)) locc<-0
		#e
		a.y <- My+locc
		a.x <- Mx+b1*(a.y-My)
		#f
		b.y <- My-locc
		b.x <- Mx+b1*(b.y-My)
	
		locc <- sqrt((C/(Lx*(1+b2^2))))
		if(is.infinite(locc) | is.na(locc)) locc <- 0
		#g
		c.y <- My+locc
		c.x <- Mx+b2*(c.y-My)
		#h
		d.y <- My-locc
		d.x <- Mx+b2*(d.y-My)

		r.grand <- sqrt((a.y-My)^2+(a.x-Mx)^2)
		r.grand <- sqrt((b.y-My)^2+(b.x-Mx)^2)
		r.petit <- sqrt((c.y-My)^2+(c.x-Mx)^2)
		r.petit <- sqrt((d.y-My)^2+(d.x-Mx)^2)
		step  <- ell.res*pi/180
		angle <- seq(0,(2*pi),by=step)
		El<-matrix(0,length(angle),2)
		El[,1] <- r.grand*sin(angle)
		El[,2] <- r.petit*cos(angle)
		mat <- eigen(cov(XY.mat))$vectors
		El <- (El%*%t(mat))
		El[,1] <- El[,1]+Mx
		El[,2] <- El[,2]+My
		if(round(var(El[,1]),digit=10)==0 | round(var(El[,2]),digit=10)==0) {
		cat(name,"Ellipse not drawn because the group had its major and/or",
			"minor axis length near zero\n")
		} else {
			lines(El[,1],El[,2],lty=ell.type,col=col,lwd=ell.lwd)
		}
		if(p.centroid==TRUE) {
			points(x=Mx,y=My,col=col,pch = centroid.type)
		}
		if(p.axis==TRUE) {
			lines(x=c(a.x,b.x),y=c(a.y,b.y),lty=ell.type,col=col,lwd=lwd)
			lines(x=c(c.x,d.x),y=c(c.y,d.y),lty=ell.type,col=col,lwd=lwd)
		}
		if(length(names)!=0)
			text(x=Mx,y=My,label=name,col=col,pos=name.pos)
	}
}
