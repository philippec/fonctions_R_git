graph.rdaTest <- function(rdaTest, xax=1, yax=2,  scaling=1,
                plot.type="notchosen", binary=NULL, bin.pch=26,
            	width=6,height=6, xlim=NULL, ylim=NULL, 
                ell=NULL, alpha=0.95, ell.axis=FALSE, 
                
              	mul.spc=0.8, mul.env=0.90, mul.text=0.10,
                pos.ell=NULL, pos.bin=NULL, pos.site=1,
                label.env=TRUE, label.spc=TRUE, label.site=TRUE, label.ell=TRUE,
               	col.env="blue", col.spc="red", col.site="black",col.ell="black",
               	lty.env=1, lty.spc=1, lty.ell=1,  lty.axis=2,
                cex=1, cex.lab=1,cex.axis=1,
                lwd=1, len=0.1, 
                type="p",
                saveplot=FALSE, path=NULL,
                mai.perc=NULL)
# Version 1.0

# Function to produce triplots from the output list of function 'rdaTest'
#
# SŽbastien Durand, UniversitŽ de MontrŽal, April 2005
#
# rdaTest: name of the rdaTest output object. The most simple call to this function
#    is to provide only the 'rdaTest' output object name. Ex.: graph.rdaTest(toto)
#
# The following parameters can also be specified by the user:
#
# xax and yax: the matrix columns to be used to plot in the rda graph.
#    These value can be negative; if so, the plot axes will be reversed
# binary: a vector specifying which of the environmental variables are binary;
#    the corresponding data points will be averaged and represented by a symbol
#    instead of an arrow, as in Canoco
# scaling: allows the user to choose between the two scaling type 1 or 2
# plot.type: either "F" or "Z" -- "F"= site scores, "Z" = fitted site scores
# saveplot: if you want to save the plot to an .pdf file;
#    if set to FALSE, the plot will only be plotted in a window (quartz);
#    If set to TRUE, no graph will be plotted on screen, and a pdf file will created
# path: complete path to the file in which the plot will be saved; 
#    example: "~/Desktop/toto.pdf"
# width and height: sizes in inches of the drawn plot
# xlim and ylim: vectors describing the minimum and maximum plotted region 
# lwd: line width for the axis, arrow, and ellipses
# len: length of the arrow heads
#
# Parameters related to the confidence ellipses around groups of points:
# ell: draw confidence ellipses around groups of points defined by a vector 
#    defining the group assignments of the objects; see example below
# alpha: confidence region of the ellipses (e.g.: 0.90, 0.80, etc.) Default value: 0.95
# lty.ell: the drawing type; see 'lty' in the help documentation file of 'par'
# col.ell: the ellipses and groupment text color; see 'lty' in the help file of 'par'
# ell.axis: if TRUE, draw the major and minor axes of the ellipses
# pos.ell: offset the group names relative to the ellipse center:
#    1=bottom, 2=left, 3=top, 4=right
# pos.bin: offset the names of the binary variables:
#    1=bottom, 2=left, 3=top, 4=right
#
# lty.env: the environmental variables arrow drawing type
# lty.spc: the species variables arrow drawing type
# col.env: the environmental variables arrow and text color
# col.spc: the species variables arrow and text color
# lty.axis: the drawing type of axis, dotted, line, .... 
# text.cex: the font size
# mul.spc: multiplier based on the range of the plot (maximum for x and y); 
#    changes proportionally to the length of the species arrows; 
#    1 is the maximum relative length in the whole set of arrows
# mul.env: the same as Mul but for the explanatory variables
# mul.text: plots the text position at the end of each arrow plus a percentage; 
#    default is 10%
# type: by default, the sites are represented as empty circles in the plot; 
#    type = "n" will not display the sites 
# mai.perc:  allows to increase or decrease the margin size, using a percentage of the 
#	   total width and height of the graph; use a value between 0 and 1
# label...:  using logicals, define whether the labels will be plotted

# Use of the function in conjunction with 'rdaTest':
#
#
# First example: use all default values
# result <- rdaTest(Y,X)
# graph.rdaTest(result)
#
# Second example: assign objects to groups; change the default values to draw ellipses
# result <- rdaTest(Y,X)
# vec=c(1,1,1,2,3,2,3,2,3,2)
# graph.rdaTest(result, plot.type="F", ell=vec, ell.axis=TRUE)
#
# Third example: change the default values to draw symbols for binary variables
# result <- rdaTest(Y,X)
# graph.rdaTest(result, plot.type="F", binary=c(2:4) )
#
# Fourth example: reverse the orientation of the axes
# result <- rdaTest(Y,X)
# graph.rdaTest(result, xax=-1, yax=-2)

{
# Axis reversion filter
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
		}else{}
		pdf(file = path, onefile = TRUE,paper="special",width=width,height=height, family="Helvetica")
	}else{
		if(length(grep("darwin",R.Version()$os))!=1){
			windows(w=width,h=height)
		}else{
			quartz(w=width,h=height)
		}
	}	
	old.par <- par(no.readonly = TRUE) # all par settings which could be changed.
 	on.exit(par(old.par))
	
	
	xinv <- 1
	if(xax < 0) {
		xinv <- (-1)
		xax  <- -(xax)
	} else {}
	yinv <- 1
	if(yax < 0) {
		yinv <- (-1)
		yax  <- -(yax)
	} else {}
	if(xax > ncol(rdaTest$F))
		stop("Their are not enough canonical axes; change the xax value")
	if(yax > ncol(rdaTest$F))
		stop("Their are not enough canonical axes; change the yax value")

# cat("Please choose a file with the different partitions","\n") 
# g<-file.choose(new=FALSE)
# read.table(g)
# cat("Please the number of groups to be used","\n") 
# partition<-tolower(scan(file="",what="character",nlines=1,quiet=T))
# group.vec<-g[,partition]
		
# Scaling 	
	if(scaling != 1) {
		rdaTest$F<-rdaTest$FSc2
		rdaTest$Z<-rdaTest$ZSc2
		rdaTest$U<-rdaTest$USc2
		rdaTest$biplotScores1<-rdaTest$biplotScores2
	} else {}

	plot.type<-toupper(as.character(plot.type))
	if(plot.type=="notchosen" | plot.type!="F" & plot.type!="Z"){
		while(plot.type!="F" & plot.type!="Z"){
			cat("Choose a functionnal plot.type: type 'Z' for fitted site scores",
				" or 'F' for site scores","\n") 
			plot.type<-toupper(scan(file="",what="character",nlines=1,quiet=T))
		}
	} else {}
	
	if(plot.type != "Z") {
		rdaTest$Z<-rdaTest$F
	} else {}
	
	if(!is.null(binary)) {
		binvar<-matrix(NA,length(binary),ncol(rdaTest$Z))
		j=1
		for(i in binary) {
			tmp<-which(rdaTest$X.mat[,i]!=0)
			binvar[j,]<-apply(rdaTest$Z[tmp,]*rdaTest$X.mat[tmp,i],2,mean)
			j=j+1
		}
		colnames(binvar)<-colnames(rdaTest$Z)
		rownames(binvar)<-colnames(rdaTest$X.mat)[binary]
		# rownames(binvar)<-c("Ra","Ri","Rs","Lv","Ln","Ll","Ls","Lf","Sc","Iba","Vl",
		#    "Vp","Vh","Vc","Tw","Tb","Hd","At","As","Pb","Li","Cl","Ch","Cb")
		# rownames(binvar)<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n",
		#    "o","p","q","r","s","t","u","v","w","x")
	}
	# rownames(rdaTest$U)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13",
	#    "14","15","16","17","18","19","20","21")

	# Text adaptation
	mul.text=mul.text+1
	
	type<-tolower(as.character(type))
	if(type!="p" & type!="n") {
		while(type!="n" & type!="p") {
			cat("Please choose a functional character for 'type': use 'p' to display",
				" sites or 'n' to not display them","\n") 
			type<-tolower(scan(file="",what="character",nlines=1,quiet=T))
		}
	} else {}	

# Check and select plot.type and application
	
	#cc<-range(xinv*rdaTest$U[,xax]);dd<-range(yinv*rdaTest$U[,yax])
	# type="n"
	#Control the marging size
	if(length(mai.perc)!=0)
		par(mai=c(height*mai.perc,width*mai.perc,height*mai.perc/2,width*mai.perc/2))
		
	plot(xinv*rdaTest$Z[,xax],yinv*rdaTest$Z[,yax],xlim =xlim,ylim=ylim,asp=1,cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
		xlab=paste("Canonical axis",xax),ylab=paste("Canonical axis",yax),type=type)
	# Print zero axis
	abline(h = 0, lty = lty.axis, lwd=lwd)
    abline(v = 0, lty = lty.axis, lwd=lwd)
	if(label.site==TRUE)
		text(xinv*rdaTest$Z[,xax],yinv*rdaTest$Z[,yax],rownames(rdaTest$Z), col = col.site, cex = cex,pos=pos.site)
	
	mult.spc=multi.factor(rdaTest$Z,rdaTest$U,xax=xax,yax=yax,xinv=xinv,yinv=yinv,
			percentage=mul.spc)
	mult.env<-multi.factor(rdaTest$Z,rdaTest$biplotScores1,xax=xax,yax=yax,xinv=xinv,
		yinv=yinv,percentage=mul.env)
	# Drawing of species variable arrows
	arrows(x0=0,y0=0, x1=xinv*rdaTest$U[,xax]*mult.spc, y1=yinv*rdaTest$U[,yax]*mult.spc, 
		col=col.spc,code=2, lty=lty.spc, len=len,lwd=1)
	if(label.spc==TRUE)
		text(x=mul.text*mult.spc*xinv*rdaTest$U[,xax],y=mul.text*mult.spc*yinv*rdaTest$U[,yax], 
				rownames(rdaTest$U), col = col.spc, cex = cex)
	
	if(!is.null(binary)) {
	# Draw environmental variable points since they are binary variables	
		points(x=xinv*binvar[,xax],y=yinv*binvar[,yax], pch=10,col=col.env)
		if(label.env==TRUE)
			text(x=xinv*binvar[,xax],y=yinv*binvar[,yax], label=rownames(binvar), 
				col=col.env,pos=pos.bin,cex=cex)
 		remain<-seq(1, ncol(rdaTest$X.mat),by=1)[-binary]
	}else{
		remain<-seq(1, ncol(rdaTest$X.mat),by=1)	
	}
	if(length(remain)!=0) {
		arrows(x0=0,y0=0,xinv*rdaTest$biplotScores1[remain,xax]*mult.env,
			yinv*rdaTest$biplotScores1[remain,yax]*mult.env, lty=lty.env,col=col.env, 
			code=2, len=len,lwd=1)
		if(label.env==TRUE)
			text(x=mul.text*mult.env*xinv*rdaTest$biplotScores1[remain,xax], 
				y=mul.text*mult.env*yinv*rdaTest$biplotScores1[remain,yax], 
				rownames(rdaTest$biplotScores1)[remain], col = col.env,cex=cex)
    } else {}	
	# Draw environmental variable arrows	
	len=9
	axis(3, at = seq(-mult.env,mult.env,length=len), labels =  round(seq(-1,1,length=len),digit=2), col = col.env,cex.axis=cex.axis)
   	axis(4, at = seq(-mult.env,mult.env,length=len), labels =  round(seq(-1,1,length=len),digit=2), col = col.env,cex.axis=cex.axis)

	if(!is.null(ell)) {
		ell<-as.factor(ell)
		name.ell=levels(as.factor(ell))
		for(i in name.ell) {
			numb<-which(ell==i)
			if(label.ell==TRUE) {
				ellipse(XY.mat=cbind(xinv*rdaTest$Z[numb,xax],yinv*rdaTest$Z[numb,yax]),
					alpha=alpha,ell.type=lty.ell,ell.lwd=lwd, p.axis=ell.axis,
					name.pos=pos.ell,col=col.ell,name=i)
			} else {
				ellipse(XY.mat=cbind(xinv*rdaTest$Z[numb,xax],yinv*rdaTest$Z[numb,yax]),
				alpha=alpha,ell.type=lty.ell,ell.lwd=lwd, p.axis=ell.axis,
				name.pos=pos.ell,col=col.ell,name=NULL)
			}
		}
	} else {}
	# Close the drawn plot
	if(saveplot==TRUE)
		dev.off()
}

multi.factor<-function(plot.m, arrow.m,xax=1,yax=2,xinv=1,yinv=1,percentage=0.80){
# Scale to the same proportions two drawings, such as plot and arrow, to make 
#    sure all arrows will be seen
# Get the matrix column ranges
	aa<-xinv*range(plot.m[,xax])
	bb<-yinv*range(plot.m[,yax])
	cc<-xinv*range(arrow.m[,xax])
	dd<-yinv*range(arrow.m[,yax])
# Find the arrow requiring more correction if exceeding the plot size,
# or find the arrow closest to the plot border to extend to span of the arrow set
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
	
	#sel=3
	border1<-c(aa[1],aa[2],bb[1],bb[2])
	border2<-c(cc[1],cc[2],dd[1],dd[2])
#	mult=0
	mult<-((border1[sel]*percentage)/border2[sel])#/percentage
	#print(aa);print(bb);print(cc);print(dd);print(dif);print(sel);print(border1);
	#print(border2);print(mult);print("fin")
 	#lines(x=c(aa[1],aa[2],aa[2],aa[1],y=aa[1]),y=c(bb[1],bb[1],bb[2],bb[2],bb[1]))
	#lines(x=c(cc[1],cc[2],cc[2],cc[1],y=cc[1]),y=c(dd[1],dd[1],dd[2],dd[2],dd[1]))
	return(mult)
}
ellipse<-function(XY.mat, alpha=0.95, ell.res=5, ell.lwd=1, ell.type=3,
				p.axis=TRUE, p.centroid=FALSE, centroid.type=20,
				col="black", name=NULL, name.pos=NULL, lwd=1) {
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
	# name.pos: the position relative to the centroid pos: values of '1', '2', '3'
    #    and '4', respectively, indicate positions below, to the left of, above,  
    #    or to the right of the specified coordinates.
    #
	if(!is.logical(c(p.axis,p.centroid))){
		cat("Choose either 'FALSE' or 'TRUE' for both 'p.axis' and 'p.centroid':","\n")
	}
	step=5
	# Compute the 95% confidence ellipse 95% if alpha = 0.95
	n<-length(XY.mat[,1])
	# Variance
	Sx<-var(XY.mat[,1]);Sy<-var(XY.mat[,2])
	# Mean
	Mx<-mean(XY.mat[,1]);My<-mean(XY.mat[,2])
	Eig<-eigen(cov(cbind(XY.mat)))$values
	
	# Covariance
	COVxy<-cov(XY.mat)
	Lx=Eig[1]
	if(Lx<0) Lx=0 
	
	Ly=Eig[2]
	if(Ly<0) Ly=0
	# Slopes of major and minor axes
	b1<-COVxy[1,2]/(Lx-Sx)
	if(is.infinite(b1) | is.na(b1)) b1<-0
	b2<--1/b1
	if(is.infinite(b2) | is.na(b1)) b2<-0
	
	if(n<=2){
		cat("Not enough individuals in group : ",name,"\n")
	}else{
		C<-((Lx*Ly)*((n-1)*2)/(n-2)*(qf(alpha,2,n-2)))
		
		locc<-sqrt(C/(Ly*(1+b1^2)))
		if(is.infinite(locc) | is.na(locc)) locc<-0
		#e
		a.y<-My+locc
		a.x<-Mx+b1*(a.y-My)
		#f
		b.y<-My-locc
		b.x<-Mx+b1*(b.y-My)
	
		locc<-sqrt((C/(Lx*(1+b2^2))))
		if(is.infinite(locc) | is.na(locc)) locc<-0
		#g
		c.y<-My+locc
		c.x<-Mx+b2*(c.y-My)
		#h
		d.y<-My-locc
		d.x<-Mx+b2*(d.y-My)

		r.grand<-sqrt((a.y-My)^2+(a.x-Mx)^2)
		r.grand<-sqrt((b.y-My)^2+(b.x-Mx)^2)
		r.petit<-sqrt((c.y-My)^2+(c.x-Mx)^2)
		r.petit<-sqrt((d.y-My)^2+(d.x-Mx)^2)
		step=ell.res*pi/180
		angle<-seq(0,(2*pi),by=step)
		El<-matrix(0,length(angle),2)
		El[,1]=r.grand*sin(angle)
		El[,2]=r.petit*cos(angle)
		mat<-eigen(cov(XY.mat))$vectors
		El<-(El%*%t(mat))
		El[,1]<-El[,1]+Mx
		El[,2]<-El[,2]+My
		if(round(var(El[,1]),digit=10)==0 | round(var(El[,2]),digit=10)==0) {
			cat(name,"ellipse was not drawn because the group had its major and/or",
			    "minor axis with length near zero\n")
		}else{
			lines(El[,1],El[,2],lty=ell.type,col=col,lwd=ell.lwd)
		}
		if(p.centroid==TRUE){
			points(x=Mx,y=My,col=col,pch = centroid.type)
		}else{}
		if(p.axis==TRUE){
			lines(x=c(a.x,b.x),y=c(a.y,b.y),lty=ell.type,col=col,lwd=lwd)
			lines(x=c(c.x,d.x),y=c(c.y,d.y),lty=ell.type,col=col,lwd=lwd)
		}else{}
		if(length(names)!=0)
			text(x=Mx,y=My,label=name,col=col,pos=name.pos)
	}
}
