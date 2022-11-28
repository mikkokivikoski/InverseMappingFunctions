##
## A toy example of inferring the probability mass function of bivalent crossover counts from gametic crossovers, and 
## utilizing those information to predict the recombination frequency for a marker pair with the p0(k) function.
##  
##

# Step 1: Estimate the likelihoods for the bivalent crossovers.
source("../Functions/emAlgorithmForCoProbabilities.r") 
crossover.data=read.table("gametic_crossovers.csv",header=T,stringsAsFactors=F,sep="\t")

maximum.co.count=max(crossover.data$Crossovers)
MLEstimates=data.frame()
for (chr in unique(crossover.data$Chromosome)){
    gametic.crossovers=concatenate.data(subset(crossover.data, Chromosome == chr)$Crossovers) 
    out.vector=main(DATA=gametic.crossovers, ITER=1000, r=100,ITER.CONF=100, ALPHA=0.025, maximum.co.count=maximum.co.count)
    out.data.frame=data.frame(c(chr,length(gametic.crossovers)-1,out.vector))
    colnames(out.data.frame)=c("Chromosome","MaximumCOcount",paste0("n",0:maximum.co.count), paste0("MLEp",0:((2*maximum.co.count)-1)), paste0("MLEpRestricted",0:((2*maximum.co.count)-1)),"BootstrapPvalue","LowerBound","UpperBound")
    MLEstimates=rbind(MLEstimates,out.data.frame)
}
print(MLEstimates)

# Step 2: Predict the likelihoods for the bivalent crossovers, utilizing the estimates derived above.
source("../Functions/p0kFunction.r")
marker.data=read.table("markers.csv",header=T,stringsAsFactors=F,sep="\t")
for(i in 1:nrow(marker.data)){
	chr.tmp=marker.data[i,"Chromosome"]
	d.tmp=marker.data[i,"MapLength"]
	mi.tmp=marker.data[i,"MarkerA"]
	mj.tmp=marker.data[i,"MarkerB"]
	p_0.tmp=as.numeric(MLEstimates[MLEstimates$Chromosome==chr.tmp,"MLEp0"])
	pk.tmp=as.numeric(MLEstimates[MLEstimates$Chromosome==chr.tmp,paste0("MLEp",1:3)])
	n.tmp=as.numeric(MLEstimates[MLEstimates$Chromosome==chr.tmp,"MaximumCOcount"])
	p0=p0k(p_0=p_0.tmp, pk=pk.tmp, n=n.tmp, d=d.tmp, mi=mi.tmp, mj=mj.tmp)
	marker.data[i,"p0"]=p0
	marker.data[i,"RecombinationFrequency_p0k"]=0.5*(1-p0)
	marker.data[i,"RecombinationFrequency_Haldane"]=0.5*(1-exp(-2*0.01*(mj.tmp-mi.tmp)))
	marker.data[i,"RecombinationFrequency_Kosambi"]=0.5*tanh(2*0.01*(mj.tmp-mi.tmp))
	marker.data[i,"RecombinationFrequency_linear"]=min(c(0.5,0.01*(mj.tmp-mi.tmp)))
}
print(marker.data)
