# from argument delete from the first one until it meets --args
commandArgs()[-(1:which(commandArgs()=="--args"))]->raw.args 
titles <- sub("=.*","",raw.args)
data <- sub(".*=","",raw.args)
names(data) <- titles

p_size = as.numeric(data[2])
print(titles)
print(p_size)

print(  commandArgs()[-(1:which(commandArgs()=="--args"))] )
norm='null'

for(i in 1:length(data))
       assign(names(data)[i],data[i])

# default Histone value:
#H = try(get("Histone"))  #get a value of a variable with a name "Histone"
#if (inherits(H,"try-error"))
Histone = read.table("name.cell")
Title = as.matrix(Histone[,1])
#print(Title)

# Quantile values are used to determine how to filter by each mark; you can specify the % filtering that each one does (q=1 is no filtering)
#QuantileValues = as.numeric(c(qH3K4me1,qH3K4me2,qH3K4me3,qH3K9ac,qK3K9me1,qH3K27ac,qH3K27me3,qH3K36me3,qH4K20me1,qCTCF) )
QuantileValues= rep(1,length(Title))
outputTitle = ''

gp =2 

#dYMAX = Histone[,2] #rep(10,10)
#type = tail(unlist(strsplit(infile,"\\.")),1) # get extension
#tfname = sub(paste(".*/(.*)\\.",type,"$",sep=""),"\\1",infile)
fname = tail(unlist(strsplit(infile,"/")),1)
celltype = unlist(strsplit(fname,"\\."))[1]
part = unlist(strsplit(fname,"\\."))[2]
#print(c(celltype,part))

data_org = read.table(infile,comment.char=">",fill=TRUE,na.strings="<")
data_p = na.exclude(data_org)

if (( dim(data_p)[1] %% p_size)!=0 ){
	print (dim(data_p)[1] %% p_size)
	print ("size does not match")
	print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	div = dim(data_p)[1] %/% p_size
	#print(dim(data_p))
	data_p = data_p[1:(p_size*div),]
	#print(dim(data_p))
}

if (dim(data_p)[2]!=length(Title)) { stop("Dimension does not match") }
numMarks = length(Title) 
imap = seq(1:numMarks)

# First get max values for each transcript
maxes = matrix(0,nc=numMarks,nr=nrow(data_p)/p_size)
dataMax = numeric(numMarks)
keep = matrix(TRUE,nc=ncol(maxes),nr=nrow(maxes))

for (i in 1:length(imap)) {	
	dat2 = data_p[,i]
	#print(dim(data_p))
	mat2 =matrix(dat2,nc=p_size,byrow=T)
	maxes[,i] = apply(mat2,1,max)		# Max across all positions for each gene
	dataMax[i] = max(maxes[,i])+1		# Record the max over all genes of this mark (for plotting)
}	
for (i in 1:length(imap)) {	
	keep[,i]  = maxes[,i] >= quantile(maxes[,i],1-QuantileValues[i])	# Keep those genes where the max is above the min cutoff
	# Record these vars in output name
	if (QuantileValues[i] != 1)
		outputTitle= paste(outputTitle,"_q",i,"_",QuantileValues[i],sep="")
}
print(outputTitle)
if (norm != 'none')
	outputTitle= paste(outputTitle,"_norm_",norm,sep="")

keep2draw = apply(keep,1,all) # apply 'all' function to keep only rows that satisfy all criteria
len = length(which(keep2draw))

#pngname = paste("Results/",fname,".png",sep="")
pngname = paste(infile,".png",sep="")
png(pngname,width=1200, height=800)
par(mfcol=c(2,7),oma=c(1,1,4,1))

for (it in 1:length(imap)) {
	i = imap[it]
	
	print(c(i,Title[i]))
	if (i<0) {  
		frame(); frame();
	}
	else {
	
		dat2 = data_p[,i]
		mat2 =matrix(dat2,nc=p_size,byrow=T)
		mat2 = mat2[keep2draw,]
		plotmax = dataMax[i]
		#if(norm=='rowMax'){
			t_max = max(maxes[,it])
			t_max = quantile(mat2,seq(0,1,by=0.001))[1000]# max(maxes[,it])
			
			#mat3 = t(apply(cbind(maxes[,it],mat2),1,function(x)x[-1]/x[1]))			
			mat3 = t(apply(cbind(maxes[,it],mat2),1,function(x)x[-1]/t_max))
			#T =  mat3 >1	# Keep those genes where the max is above the min cutoff
			#mat3[which(T)]=1
			mat3[mat3>0.3]=0.3
			plotmax = 0.3
		#}

		plotmax = 0.98 
		t_max=15
		mat3 = t(apply(cbind(maxes[,it],mat2),1,function(x)x[-1]/t_max))
		mat3[mat3>plotmax]=plotmax
		

		par(mar=c(0,gp+1,gp+5,0))		

		#mat3=mat2
		# Heatmap		
		image(t(mat3),col=rgb(seq(0,1,.01),0,0),zlim=c(0,plotmax),main=Title[i],axes=FALSE,cex.main=gp*1.9)
		box()
				
		## Average plot		 	
		ll_sum = colMeans(mat2,na.rm=T)
		
		tmpmax =15 
		#if (it>4) tmpmax=70
		par(mar=c(gp+4, gp+1,2, 0))					
		plot(ll_sum,ylim=c(0,tmpmax),type='l',xlim=c(1,ncol(mat2)),ylab="count",xlab="",xaxt="n",axes=FALSE)
		box()
		axis(2,at=c(0,tmpmax),cex.axis=gp)
		axis(1,at=c(2,(p_size-1)/2,p_size-2),cex.axis=gp,
		label=c(paste("     ",-(p_size)/100,"Kbp",sep=""),"0", paste((p_size)/100,"Kbp",sep="") ))
	} 
}
par(mar=c(0,gp+2,gp+1,0))
y<-seq(0,1,0.01)
x<-seq(0,10,1)

plot(0,0,xlim=c(0,10), ylim=c(0,1),axes=FALSE)
for (k in 1:1) {
	for (l in 1:100) {points(x[k],y[l],pch=7,col=rgb(y[l],0,0) )}
}
 axis(2, at=c(0,0.5,1))
       
dev.off()
