GetChromoStartEnd <- function(region){
  # "chr10:67,536,192-67,540,761" -> chr10, start, end as list
  jchromo <- strsplit(region, ":")[[1]][[1]]
  jstartend <- strsplit(region, ":")[[1]][[2]]
  jstart <- gsub(pattern = ",", replacement = "", strsplit(jstartend, "-")[[1]][[1]])
  jend <- gsub(pattern = ",", replacement = "", strsplit(jstartend, "-")[[1]][[2]])
  return(list(chromo = jchromo, start = jstart, end = jend))
}

MakeCoords <- function(x){
  # chr7  17058086  17058738 -> chr7:17058086-17058738
  xsplit <- strsplit(x, " ")[[1]][c(1, 3, 5)]
  return(paste0(xsplit[[1]], ":", xsplit[[2]], "-", xsplit[[3]]))
}

RenameGene<-function(GeneName){
  GeneInfo<-GenePos[GenePos$V16 == GeneName,]
  chr<-paste("chr",GeneInfo$V1,sep="")
  Pos <- paste0(GeneInfo$V4, " ", GeneInfo$V5)
  # Pos<-round(mean(c(GeneInfo$V4,GeneInfo$V5))/1000000,4)
  return(paste0(GeneName," ",chr, " ", Pos))
}

RenamePeak<-function(PeakID){
  PeakInfo<-PeakPos[PeakPos$V10==PeakID,]
  chr<-gsub("chr",'',PeakInfo$V1)
  chr<-paste("chr",chr,sep="")
  # Pos<-round(mean(c(PeakInfo$V4,PeakInfo$V5))/1000000,4)
  Pos <- paste0(PeakInfo$V4, " ", PeakInfo$V5)
  return(paste0("Peak ",chr, " ", Pos))
}

CalcDist <- function(GeneName, PeakID, units = "kb"){
  GeneInfo<-GenePos[GenePos$V16 == GeneName,]
  PeakInfo<-PeakPos[PeakPos$V10==PeakID,]

  # take min and max because Gene name can be assigned to multiple transcript. Take largest possible distance
  # length can be llarger than 1!
  GenePos.start <- GeneInfo$V4
  GenePos.end <- GeneInfo$V5

  PeakPos.start <- PeakInfo$V4
  PeakPos.end <- PeakInfo$V5
  PeakPos.center <- (PeakPos.start + PeakPos.end) / 2
  # calculate minimum distance relative to gene!
  left.or.right <- sign(GenePos.start - PeakPos.center)  # left is negative, right is positive
  dist.min.abs <- min(GenePos.start - PeakPos.center, GenePos.end - PeakPos.center)
  
  dist.min <- min(left.or.right * dist.min.abs)
  if (units == "kb"){
    dist.min <- round(dist.min / 1000)
  } else {
    print("Warning: Units only coded in kb, doing nothing...")
  }
  return(dist.min)
}

PlotTimeCourse<-function(Peak,Tr,.log=TRUE, .scale = TRUE, .center = TRUE, eps = 1, .ymin = -2, .ymax = 2, append.title = "", add.points = FALSE, add.rect = TRUE, TrExp.override = NULL, legend.pos = "topright"){
  # If TrExp.override then it should be a named vector. Names are times. Vector values are expression (in linear counts)
  #par(mfrow=c(2,2))


  if (!is.null(TrExp.override)){
    TrExp <- as.numeric(TrExp.override)
    names(TrExp) <- names(TrExp.override)
  }

  Time<-colnames(cpmvAtac)

  # THESE WERE REMOVED WHY???
  # idx<-which(Time %in% c(Time[grep('42S',Time)],Time[grep('70S',Time)],Time[grep('76S',Time)]))
  idx<-which(Time %in% c(Time[grep('70S',Time)],Time[grep('76S',Time)]))
  Time<-Time[-idx]

  PeakExp<-cpmvAtac[Peak,]
  PeakExp<-as.numeric(PeakExp[-idx])

  if (is.null(TrExp.override)){
    TrExp<-HCdGenes[Tr,]
    TrExp<-as.numeric(TrExp[-idx])
  }
  if (.log){
    PeakExp <- log2(PeakExp + eps)
  }

  
  Timey<-Time
  Timey[grep("T12N",Time)]<-"12"
  Timey[grep("T12S",Time)]<-"36"
  Timey[grep("T18N",Time)]<-"18"
  Timey[grep("T18S",Time)]<-"42"
  Timey[grep("T24S",Time)]<-"48"
  Timey[grep("T24N",Time)]<-"24"
  Timey[grep("T30S",Time)]<-"54"
  Timey[grep("T36S",Time)]<-"60"
  Timey[grep("T3N",Time)]<-"3"
  Timey[grep("T3S",Time)]<-"27"
  Timey[grep("T48S",Time)]<-"72"
  Timey[grep("T54S",Time)]<-"78"
  Timey[grep("T6S",Time)]<-"30"
  Timey[grep("T6N",Time)]<-"6"
  if (!is.null(TrExp.override)){
    # Set ZT0H and ZT0I to be ZT24: this affects only ATACseq
    Timey[grep("ZT0[B,C]",Time)]<-"0"
    Timey[grep("ZT0[H,I]",Time)]<-"24"
    # KEEP J70N1 and J70N2
    Timey[grep("J70N1",Time)]<-"0"
    Timey[grep("J70N2",Time)]<-"24"
    Timey[grep("J76N",Time)]<-"6"
    # Only ONE time point? I thought there were 3!
    Timey[grep("T42S",Time)]<-"66"
  } else {
    # Default put all 4 in T=0
    Timey[grep("ZT0",Time)]<-"0"
  }
 
  # This is not used ever 
  # Timep<-Time
  # Timep[grep("T12N",Time)]<-"12"
  # Timep[grep("T12S",Time)]<-"12"
  # Timep[grep("T18N",Time)]<-"18"
  # Timep[grep("T18S",Time)]<-"18"
  # Timep[grep("T24S",Time)]<-"24"
  # Timep[grep("T24N",Time)]<-"24"
  # Timep[grep("T30S",Time)]<-"30"
  # Timep[grep("T36S",Time)]<-"36"
  # Timep[grep("T3N",Time)]<-"3"
  # Timep[grep("T3S",Time)]<-"3"
  # Timep[grep("T48S",Time)]<-"48"
  # Timep[grep("T54S",Time)]<-"54"
  # Timep[grep("T6S",Time)]<-"6"
  # Timep[grep("T6N",Time)]<-"6"
  # Timep[grep("ZT0",Time)]<-"0"

  if (is.null(TrExp.override)){
    names(TrExp)<-Timey
  }
  names(PeakExp)<-Timey


  
  # Not used ever
  # Cond<-rep(8,length(Time))
  # Cond[grep("S",Time)]<-9

  # print("BEFORE")
  # print(PeakExp)
  # print(length(PeakExp))
  PeakExp<-as.vector(scale(PeakExp, center = .center, scale = .scale))
  names(PeakExp) <- Timey
  # print("AFTER")
  # print(PeakExp)
  names.before <- names(TrExp)
  TrExp<-as.vector(scale(TrExp, center = .center, scale = .scale))
  names(TrExp) <- names.before
  # names(TrExp) <- Timey
  
  TrM<-aggregate(TrExp,list(names(TrExp)),mean)
  PeakM<-aggregate(PeakExp,list(names(PeakExp)),mean)

  PeakM<-PeakM[order(as.numeric(PeakM$Group.1)),]
  TrM<-TrM[order(as.numeric(TrM$Group.1)),]

  # # scale before the mean for flexibility
  # PeakM$x<-scale(PeakM$x, center = .center, scale = .scale)
  # TrM$x<-scale(TrM$x, center = .center, scale = .scale)
  

  if (add.points){
    TrPts<-data.frame(Group.1 = as.numeric(names(TrExp)), x = TrExp)
    PeakPts<-data.frame(Group.1 = as.numeric(names(PeakExp)), x = PeakExp)
  }

  if (!add.points){
    minval<-min(c(TrM$x,PeakM$x))
    maxval<-max(c(TrM$x,PeakM$x))
  } else {
    minval<-min(c(TrPts$x,PeakPts$x))
    maxval<-max(c(TrPts$x,PeakPts$x))
  }
  
  # minmax<-abs(c(minval,maxval))

  if (maxval > .ymax){
    .ymax <- maxval
  } 
  if (minval < .ymin){
    .ymin <- minval
  }

  # minmax<-c(minval,maxval)
  # if (max(minmax) < .ymax ){
  #   minval <- .ymin
  #   maxval<- .ymax
  # } else {
  #   minval <- max(minmax)
  #   maxval <- max(minmax)
  # }
  
  par(mar=c(5,5,5,3))

  jdist <- CalcDist(Tr, Peak, units = "kb")
  jtitle <- paste0(RenameGene(Tr),"\n",RenamePeak(Peak),"\n", "Distance: ", jdist, " kb")
  jylab <- "RNA-seq or ATAC-seq levels"
  if (.scale){
  	jylab <- paste0(jylab, "\n(Scaled)")
  } else if (.scale == FALSE & .center == TRUE){
  	jylab <- paste0(jylab, "\n(centered log2 counts)")
  } else {
  	jylab <- paste0(jylab, "\n(log2 counts)")
  }
  jtitle <- paste0(jtitle, append.title)

  plot(PeakM$Group.1,PeakM$x,type='l', lwd = 5, col="blue",ylim=c(.ymin, .ymax),xlab="Time [h]",ylab=jylab, main=jtitle, xaxt="n")
  axis(1, at=seq(0,90,by=6),labels=seq(0,90,by=6), las=1)
  # rect(12,-maxval-.5,24,maxval+.5,col=rgb(0,0,0,.2))
  # rect(36,-maxval-.5,48,maxval+.5,col=rgb(0,0,0,.2))
  # rect(60,-maxval-.5,72,maxval+.5,col=rgb(0,0,0,.2))

  rect(12, .ymin - 0.2 * abs(.ymin), 24, .ymax + 0.2 * abs(.ymax), col=rgb(0,0,0,.2))
  rect(36, .ymin - 0.2 * abs(.ymin), 48, .ymax + 0.2 * abs(.ymax), col=rgb(0,0,0,.2))
  rect(60, .ymin - 0.2 * abs(.ymax), 72, .ymax + 0.2 * abs(.ymax), col=rgb(0,0,0,.2))
  # rect(12,.ymin-.5,24,.ymax+.5,col=rgb(0,0,0,.2))
  # rect(36,.ymin-.5,48,.ymax+.5,col=rgb(0,0,0,.2))
  # rect(60,.ymin-.5,72,.ymax+.5,col=rgb(0,0,0,.2))
  points(PeakM$Group.1,PeakM$x,pch=19,cex=.8,col="blue")
  
  lines(c(c(0,3,6,12,18),c(0,3,6,12,18)+24,c(0,3,6,12,18)+48,c(0,3,6,12,18)+72),rep(TrM$x[1:5],4),col=rgb(1,0,0,0.5),lty=2)
  points(c(c(0,3,6,12,18),c(0,3,6,12,18)+24,c(0,3,6,12,18)+48,c(0,3,6,12,18)+72),rep(TrM$x[1:5],4),col=rgb(1,0,0,0.5),pch=19,cex=.5)
  
  lines(c(c(0,3,6,12,18),c(0,3,6,12,18)+24,c(0,3,6,12,18)+48,c(0,3,6,12,18)+72),rep(PeakM$x[1:5],4),col=rgb(0,0,1,0.5),lty=2)
  points(c(c(0,3,6,12,18),c(0,3,6,12,18)+24,c(0,3,6,12,18)+48,c(0,3,6,12,18)+72),rep(PeakM$x[1:5],4),col=rgb(0,0,1,0.5),pch=19,cex=.5)
  
  lines(TrM$Group.1,TrM$x,col="red", lwd = 5)
  points(TrM$Group.1,TrM$x,pch=19,cex=.8,col="red")

  if (add.points){
    points(TrPts$Group.1, TrPts$x, col = "red")
    points(PeakPts$Group.1, PeakPts$x, col = "blue")
  }

  abline(v=c(0,24,48,72),col=rgb(0,0,0,0.5))
  if (add.rect){
    rect(24,.ymin,30,.ymin * 1.1, col=rgb(1,0,0,0.7))
  }

  legend(legend.pos, legend=c("RNAseq", "ATACseq"),
	        col=c("red", "blue"), cex=0.8, lty = 1)
  
}
