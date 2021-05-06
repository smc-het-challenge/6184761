args <- commandArgs(T);

tb <- read.table(args[1],header=T,stringsAsFactors=F,sep="\t");
tb$Chr <- as.integer(sub("X","23",sub("Y","24",sub("MT*","25",sub("chr","",tb$chr)))));
tb$w <- tb$endpos-tb$startpos+1; tb$w<-tb$w/sum(tb$w); tb$cnA <- tb$nMaj1_A+tb$nMin1_A; tb$maf <- 1-tb$BAF;
tm <- read.table(pipe(paste("sed -e's/^#//g' ",args[2])),header=T,stringsAsFactors=F);
tm$Chr <- as.integer(sub("X","23",sub("Y","24",sub("MT*","25",sub("chr","",tm$CHROM)))));
tm$nRef <- as.integer(sub("^[^:]*:([0-9]*),.*","\\1",tm$tumor));
tm$nAlt <- as.integer(sub("^[^:]*:[0-9]*,([0-9]*):.*","\\1",tm$tumor));
tm$nDep <- tm$nRef+tm$nAlt;

batt.cell.lm <- function(b) {
	if(nrow(b)<1) return(NA);
	ix <- which(is.na(b$nMaj2_A));
	if(length(ix)==0) return(NA);
	return(max(min(1,lm(I(2*2**b$LogR[ix])~I((b$nMaj1_A+b$nMin1_A)[ix]),weights=(b$endpos-b$startpos+1)[ix])$coefficient[2])));
}

batt.cell.lrr <- function(b) {
	if(nrow(b)<1) return(NA);
	ix <- which(is.na(b$nMaj2_A) & b$pval > .999);
	if(length(ix)==0) return(NA);
	logr <- b$LogR[ix]; cx <- (b$nMaj1_A[ix]+b$nMin1_A[ix])/2-1; w <- b$endpos[ix]-b$startpos[ix]+1; w<-w/sum(w); mu<- .0; rho<-.5;
	for(j in 1:1000) {
		mu.n <- sum((logr-log2(1+cx*rho))*w);
		bi <- cx/(1+cx*rho);
		rho.n <- rho+sum((logr-mu.n-log2(1+cx*rho))*bi*w)/sum(bi**2*w);
		if(rho.n > 1) rho.n <- (1+rho)/2;
		if(rho.n < 0) rho.n <- rho/2;
		if(abs(mu.n/(1e-16*(mu==0)+mu)-1)+abs(rho.n/rho-1)<1e-8) break;
		rho <- rho.n; mu <- mu.n;
	}
	return(rho);
}

batt.cell.baf <- function(b) {
	if(nrow(b)<1) return(NA);
	ix <- which(is.na(b$nMaj2_A) & b$pval > .999 & (b$nMaj1_A > b$nMin1_A));
	if(length(ix)==0) return(NA);
	baf <- 1-b$BAF[ix]; nC <- (b$nMaj1_A[ix]+b$nMin1_A[ix]-2); nc <- (b$nMin1_A[ix]-1); w <- b$endpos[ix]-b$startpos[ix]+1; w<-w/sum(w); mu<- .0; rho<-.5;
	for(j in 1:10) {
		rho.n <- rho-sum((baf-(1+nc*rho)/(2+nC*rho))*(nC-2*nc)/(2+nC*rho)**2*w)/sum((nC-2*nc)**2/(2+nC*rho)**4*w);
		if(rho.n > 1) rho.n <- (1+rho)/2;
		if(rho.n < 0) rho.n <- rho/2;
		if(abs(rho.n/rho-1)<1e-8) break;
		rho <- rho.n;
	}
	return(rho);
}


batt.cns <- function(b) {
	ix <- which(is.na(b$nMaj2_A) & b$pval > .999);
	bcns <- data.frame(CN=sort(unique(b$cnA[ix])),nseg=0,w=0,mx=0,ma=0,b0=NA,b0w=NA,b1=NA,b1w=NA,b2=NA,b2w=NA,stringsAsFactors=F);
	for(j in 1:nrow(bcns)) {
		ix1 <- ix[b$cnA[ix]==bcns$CN[j]];
		bcns$nseg[j] <- length(ix1);
		bcns$w[j]  <- sum(b$w[ix1]);
		bcns$mx[j] <- sum(b$LogR[ix1]*b$w[ix1])/bcns$w[j];
		bcns$ma[j] <- sum(abs(b$LogR[ix1]-bcns$mx[j])*b$w[ix1])/bcns$w[j];
		k1 <- as.integer(floor((bcns$CN[j]+1)/2));
		for(i in 0:2) if(k1+i <= bcns$CN[j]) {
			ix2 <- ix1[b$nMaj1_A[ix1]==k1+i];
			if(length(ix2)==0) next;
			bcns[[sprintf("b%dw",i)]][j] <- aw <- sum(b$w[ix2]);
			bcns[[sprintf("b%d",i)]][j] <- sum(b$maf[ix2]*b$w[ix2])/aw;
		}
	}
	return(bcns);
}

batt.plcl <- function(b) {
	bcns <- batt.cns(b);
	mu <- NA; rho <- NA; rhot<-c(NA,NA,batt.cell.baf(b)); rvl <- NA;
	while(T) {
		if(length(ix2<-which(bcns$CN==2 & bcns$ma<.01))==1) {
			if(!is.na(bcns$b0[ix2]) && bcns$b0[ix2]>.5-.01 && bcns$b0w[ix2]>.05) {
				mu <- bcns$mx[ix2];
				if(length(ix4<-which(bcns$CN==4 & bcns$ma < .01))==1) {
					rhox <- min(1,2**(bcns$mx[ix4]-mu)-1);
					if((!is.na(bcns$b0[ix4]) && bcns$b0[ix4]>.5-.01 && bcns$b0w[ix4]>.05) ||
							(!is.na(bcns$b1[ix4]) && bcns$b1w[ix4]>.05 && abs(bcns$b1[ix4]-1/2/(1+rhox))<.01) ||
							(!is.na(bcns$b2[ix4]) && bcns$b1w[ix4]>.05 && abs(bcns$b2[ix4]-(1-rhox)/2/(1+rhox))<.01)) {
						rho <- rhox; rvl <- 1;
						break;
					}
				}
			}
			ix <- which(b$pval>.999 & is.na(b$nMaj2_A));
			if(length(unique(b$cnA[ix]))>1) {
				rhot[1] <- lm(I(2**(b$LogR[ix]-mu)-1)~I(b$cnA[ix]/2-1)-1,weights=b$w[ix])$coeff[1];
				logr <- b$LogR[ix]; cx <- b$cnA[ix]/2-1; w <- b$w[ix]; rhot[2]<-rhot[1];
				for(j in 1:1000) {
					bi <- cx/(1+cx*rhot[2]);
					rho.n <- rhot[2]+sum((logr-mu-log2(1+cx*rhot[2]))*bi*w)/sum(bi**2*w);
					if(rho.n > 1) rho.n <- (1+rhot[2])/2;
					if(rho.n < 0) rho.n <- rhot[2]/2;
					if(abs(rho.n/rhot[2]-1)<1e-8) break;
					rhot[2] <- rho.n;
				}
				rhos <- sort(rhot);
				if(rhos[3]-rhos[1]<.1) { rho<-rhos[2]; rvl<-2; break;  }
			}
		} 
		if(length(ix4<-which(bcns$CN==4 & bcns$ma<.01))==1) {
			ix <- which(b$pval>.999 & is.na(b$nMaj2_A));
			if(length(unique(b$cnA[ix]))>1 && (!is.na(bcns$b0[ix4]) && bcns$b0[ix4]>.5-.01 && bcns$b0w[ix4]>.05)) {
				mu4 <- bcns$mx[ix4];
				logr <- b$LogR[ix]; cx <- b$cnA[ix]/2-1; w <- b$w[ix]; rhot[2]<-rhot[1];
				rhot[1] <- max(0,min(1,1/lm(I(2**(logr-mu4)-cx)~I(1-cx)-1,weights=b$w[ix])$coeff[1]-1));
				logr <- b$LogR[ix]; cx <- b$cnA[ix]/2-1; w <- b$w[ix]; rhot[2]<-rhot[1];
				for(j in 1:1000) {
					bi <- cx/(1+cx*rhot[2])-1/(1+rhot[2]);
					rho.n <- rhot[2]+sum((logr-mu4+log2(1+rhot[2])-log2(1+cx*rhot[2]))*bi*w)/sum(bi**2*w);
					if(rho.n > 1) rho.n <- (1+rhot[2])/2;
					if(rho.n < 0) rho.n <- rhot[2]/2;
					if(abs(rho.n/rhot[2]-1)<1e-8) break;
					rhot[2] <- rho.n;
				}
				rhos <- sort(rhot);
				if(rhos[3]-rhos[1]<.1) { rho<-rhos[2]; rvl <- 3; mu <- mu4-log2(1+rho); break;  }
			}
		}
		rhot[1] <- batt.cell.lm(b);
		rhot[2] <- batt.cell.lrr(b);
		rhos <- sort(rhot);
		if(rhos[3]-rhos[1]<.1) { rho<-rhos[2]; rvl <- 4; break;  }
		rho<-rhot[3]; rvl<-5; 
		break;
	}
	if(is.na(mu) && !is.na(rho) && length(ix <- which(is.na(b$nMaj2_A) & b$pval > .999))>0) mu <- sum((b$LogR[ix]-log2(1+(b$cnA[ix]/2-1)*rho))*b$w[ix])/sum(b$w[ix]);
	return(c(mu,rho,rvl));
}


batt.cnv <- function(b,mu,rho,lrr.w=5) {
	if(missing(mu) || missing(rho)) {
		bpl <- batt.plcl(b);
		if(missing(mu))  mu <- bpl[1];
		if(missing(rho)) rho<- bpl[2];
		rvl <- bpl[3];
	} else {
		rvl <- NA;
	}
	bpl <- c(mu,rho,rvl);
	bcns <- batt.cns(b);
	bn <- bc <- pmax(0,round((2**(b$LogR-mu+1)-2)/rho+2));
	bm <- rep(NA,nrow(b));
	bv <- rep(Inf,nrow(b));
	for(k in (-1):1) {
		tb <- rep(Inf,nrow(b));
		tl <- rep(0,nrow(b)); tm <- rep(NA,nrow(b));
		ix <- which(bc<=-k); if(length(ix)>0) tl[ix]<-tl[ix]+(b$LogR[ix]-mu-3.5)**2;
		ix <- which(bc> -k); if(length(ix)>0) tl[ix]<-tl[ix]+(b$LogR[ix]-mu-log2(1+((bc[ix]+k)/2-1)*rho))**2;
		ix <- which(bc<=-k); if(length(ix)>0) { tb[ix] <- (b$maf[ix]-1/2)**2; tm[ix] <- 0;  }
		for(j in which(bc> -k)) {
			for(i in 0:(floor((bc[j]+k)/2))) {
				tc <- (b$maf[j]-(1+(i-1)*rho)/(2+(bc[j]+k-2)*rho))**2;
				if(tc < tb[j]) { tb[j] <- tc; tm[j] <- i;  };
			}
		}
		tv <- tb/lrr.w+tl;
		ix <- which(tv < bv); if(length(ix)>0) { bv[ix] <- tv[ix]; bn[ix] <- bc[ix]+k; bm[ix] <- tm[ix];  };
	}
	b$cn <- bn; b$cm <- bm; b$diff <- bv; b$prl <- mu+log2(1+(b$cn/2-1)*rho); b$prb <- (1+(b$cm-1)*rho)/(2+(b$cn-2)*rho);
	return(list(cn=bn,cm=bm,df=bv,b=b,plcl=bpl,cns=bcns));
}


batt.ncl <- function(b,cns) {
	if(missing(cns)) cns <- batt.cnv(b);
	dfmx <- sum(sqrt(cns$b$diff)*cns$b$w)/sum(cns$b$w);
	dfmd <- sum(abs(sqrt(cns$b$diff)-dfmx)*cns$b$w)/sum(cns$b$w);
	ix1 <- which(cns$b$LogR > cns$plcl[1]+log2(1+((cns$cn+.05)/2-1)*cns$plcl[2]));
	ix2 <- which(cns$b$LogR < cns$plcl[1]+log2(1+((cns$cn-.05)/2-1)*cns$plcl[2]));
	ix <- c(ix1,ix2);
	if(length(ix)==0) return(1);
	k <- 2;
	tx <- cns$b[ix,c("LogR","maf","w","cn","cm")];
	ix3 <- rep(c(1,-1),c(length(ix1),length(ix2)));
	for(j in 2:k) { 
		tx[[j*2+2]] <- tx$cn+(j-1)*ix3;
		tx[[j*2+3]] <- tx$cm; 
	 }
	prop <- rep(1/k,k);
	orgS <- sum((tx$LogR - cns$plcl[1] - log2(1+(tx[[4]]/2-1)*cns$plcl[2]))**2*tx$w);

	for(itw in 1:1000) {
		propo <- prop;
		oc <- prop[1]*tx[[4]]; for(j in 2:k) oc <- oc + prop[j]*tx[[j*2+2]];
		oc <- 1+(oc/2-1)*cns$plcl[2]; A <- tx$LogR - cns$plcl[1] - log2(oc);
		oldS <- sum(A**2*tx$w);
		B <- (tx[[k*2+2]]-tx[[4]])*cns$plcl[2]/2/oc;
		pdel <- sum(A*tx$w*B)/sum(tx$w*B**2); 
		for(ith in 1:1000) {
			propn <- prop + c(-sum(pdel),pdel);
			if(min(propn) <0) { pdel <- pdel/2; next;  }
			oc <- propn[1]*tx[[4]]; for(j in 2:k) oc <- oc + propn[j]*tx[[j*2+2]];
			newA <- cns$plcl[1]+log2(1+(oc/2-1)*cns$plcl[2])-tx$LogR;
			newS <- sum(newA**2*tx$w);
			if(newS < oldS) break;
			pdel <- pdel/2;
		}
		oc <- propn[1]*tx[[4]]; for(j in 2:k) oc <- oc + propn[j]*tx[[j*2+2]];
		tx$prl <- cns$plcl[1]+log2(1+(oc/2-1)*cns$plcl[2]);
		if(sum(abs(propn/(1e-16+propo)-1)) < 1e-8) break;
		prop <- propn;
	}
	return(sort(prop,decreasing=T));
}



cns <- batt.cnv(tb);
ncl <- batt.ncl(tb,cns);
tb$clust <- 1;
if(length(ncl)==2) {
	ix1 <- which(cns$b$LogR > cns$plcl[1]+log2(1+((cns$cn+.05)/2-1)*cns$plcl[2]));
	ix2 <- which(cns$b$LogR < cns$plcl[1]+log2(1+((cns$cn-.05)/2-1)*cns$plcl[2]));
	ix <- c(ix1,ix2);
	if(length(ix)>0) tb$clust[ix] <- 2;
}
tm$clust <- NA;
for(j in 1:nrow(tb)) if(length(ix<-which(tm$Chr==tb$Chr[j] & tm$POS >= tb$startpos[j] & tm$POS <= tb$endpos[j]))>0) tm$clust[ix] <- tb$clust[j];
if(length(ix<-which(is.na(tm$clust)))>0) {
	vals <- NA; 
	for(j in 1:length(ix)) {
		if(length(ix2<-which(tm$Chr==tm$Chr[ix[j]] & !is.na(tm$clust)))>0) {
			ixd <- abs(tm$POS[ix2]-tm$POS[ix[j]]);
			ix3 <- which(ixd==min(ixd))[1];
			tm$clust[ix[j]] <- tm$clust[ix3];
		}
		if(is.na(tm$clust[ix[j]])) tm$clust[ix[j]]<-1;
	}
}

write.table(cns$plcl[2],file="ans_1A.txt",col.names=F,row.names=F,sep="\t",quote=F);
write.table(length(ncl),file="ans_1B.txt",col.names=F,row.names=F,sep="\t",quote=F);
a1c <- data.frame(ncl=1:length(ncl),nssm=0,pcl=cns$plcl[2]*c(1,ncl)[1:length(ncl)],stringsAsFactors=F);
for(j in 1:length(ncl)) a1c$nssm[j] <- sum(tm$clust==j);
write.table(a1c,file="ans_1C.txt",col.names=F,row.names=F,sep="\t",quote=F);
write.table(tm$clust,file="ans_2A.txt",col.names=F,row.names=F,sep="\t",quote=F);
if(length(ncl)==1) {
	ccm <- matrix(1,nrow(tm),nrow(tm));
    write.table(ccm,file="ans_2B.txt",col.names=F,row.names=F,sep="\t",quote=F); rm(ccm);
	adm <- matrix(0,nrow(tm),nrow(tm));
	write.table(matrix(c(1,0),nrow=1,ncol=2),file="ans_3A.txt",col.names=F,row.names=F,sep="\t",quote=F);
} else {
	ccm <- diag(1,nrow(tm)); pr <- round(ncl[1],4);
	ix1 <- which(tm$clust==1); ccm[ix1,ix1] <- 1;
	ix2 <- which(tm$clust==2); ccm[ix2,ix2] <- 1;
	ccm[ix1,ix2] <- pr; ccm[ix2,ix1] <- pr;
    write.table(ccm,file="ans_2B.txt",col.names=F,row.names=F,sep="\t",quote=F); rm(ccm);
	adm <- matrix(0,nrow(tm),nrow(tm));
	adm[ix1,ix2] <- 1-pr;
	write.table(matrix(c(1,0,2,1),nrow=2,ncol=2,byrow=T),file="ans_3A.txt",col.names=F,row.names=F,sep="\t",quote=F);
}
write.table(adm,file="ans_3B.txt",col.names=F,row.names=F,sep="\t",quote=F);



