'map.hclust' <- 
function(res, coord, k1=2, k2=NULL, title=NULL, cex=2, pos=3, col="red", ...)
# Draw maps for constrained.clust output, from k1 (default: k1=2) to k2 groups.
# Maps are drawn to a pdf file.
{
if( is.integer(coord) ) coord = cbind(coord, rep(1,length(coord)))
if( (is.matrix(coord) | is.data.frame(coord)) & ncol(coord)==1 )
	coord = cbind(coord[,1], rep(1,nrow(coord)))
if(is.null(title)) {
	if(length(res$hcl)==0) {
		pdf(file="Unconstrained clustering maps.pdf",width=15,height=15,family="Times", pointsize=20)
		} else {
		pdf(file="Constrained clustering maps.pdf",width=15,height=15,family="Times", pointsize=20)
		}
	} else {
	pdf(file=title,width=15,height=15,family="Times", pointsize=20)
	}
n <- length(res$height) + 1
if(is.null(k2)) k2 <- trunc(n/2)
if(k2 > (n-1))  k2 <- trunc(n/2)
range.X <- max(coord[,1]) - min(coord[,1])
range.Y <- max(coord[,2]) - min(coord[,2])
outer <- TRUE
if(range.Y < 0.9*range.X) outer <- FALSE
#
# cat("n =",n,"  k1 =",k1,"  k2 =",k2,"\n")
for(k in k1:k2) {
	#
	plot(coord[,1], coord[,2], type="p", asp=1)
	text(coord[,1], coord[,2], cutree(res, k), cex=cex, pos=pos, col=col)
	#
if(length(res$hcl)==0) {
	title(paste("Unconstrained clustering map, k =",k), line=-2, outer=outer)
	} else {
	title(paste("Constrained clustering map, k =",k), line=-2, outer=outer)
	}
	#
	}
dev.off()
}
