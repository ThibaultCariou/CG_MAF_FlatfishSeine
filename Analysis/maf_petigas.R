maf.f <- function(data, xy, hmin, hmax) {
  # Rscript to compute MAFs on a space-time data set
  # Pierre Petitgas (Ifremer), Didier Renard and Nicolas Desassis (Mines-ParisTech), July 2018
  # Requires library(ade4) for principle components analysis (PCA)
  # data: Matrix (nx,nt) on which to apply PCA
  # nx (lines) are the spatial locations
  # nt (columns) are the time locations
  # xy: Matrix of data coordinates [,1]=x, [,2]=y
  # hmin, hmax: Minimum, maximum lag distance at which to compute MAFs
  # Returns:
  # eig= Eigenvalues
  # x= MAFs (zero mean and unit variance)
  # c1= Normed scores of data vectors on MAFs
  ### 1. Original data are transformed into normalized PCs (y)
  prc1 = dudi.pca(data,center=T,scale=T,scannf=F,nf=dim(data)[2])
  y = as.matrix(prc1$l1) # Row normed scores
  c1.1 = as.matrix(prc1$c1) # Normed scores
  npc = dim(y)[2]
  ### 2. Increments of normalized PCs (datdif)
  matdis=sqrt(outer(xy[,1],xy[,1],"-")^2+outer(xy[,2],xy[,2],"-")^2) #Euclidean distance
  datdif=matrix(NA,nrow=sum(!(matdis<hmin | matdis>hmax)),ncol=npc)
  for (i in 1:npc) {
    dif.y = outer(y[,i],y[,i],"-")
    dif.y[matdis<hmin | matdis>hmax]=NA
    datdif[,i] = dif.y[! is.na(dif.y)]
  }
  datdif=datdif[apply(!is.na(datdif),1,all),]
  ### 3. PCA of increments, MAFs
  prc2 = dudi.pca(datdif,center=F,scale=F,scannf=F,nf=npc)
  c1.2 = as.matrix(prc2$c1)
  eig.2 = as.matrix(prc2$eig[npc:1])
  x = y %*% c1.2
  x = x[,npc:1]
  maf.c1 = c1.1 %*% c1.2
  maf.c1 = maf.c1[,npc:1]
  dimnames(x) = list(NULL,paste0("MAF",1:npc))
  dimnames(maf.c1) = list(names(data),paste0("MAF",1:npc))
  #
  return(list(eig=eig.2, x=x, c1=maf.c1))
}