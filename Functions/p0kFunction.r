# R implementation of the p_0(k) function described in Kivikoski et al. (2022)

# INPUT PARAMETERS
# p_0: likelihood for no crossovers in meiosis (Estimated from data or 0 if obligate crossover is assumed)
# pk: a vector of likelihoods for 1...n crossovers in the bivalent (Estimated from data)
# n: the highest possible number of bivalent crossovers (Estimated from data)
# d: map length of the chromosome (cM)
# mi: map position of marker 1 (cM) 
# mj: map position of marker 2 (cM)

# OUTPUT
# p0.out: Probability of no bivalent crossovers between the markers

p0k=function(p_0,pk,n,d,mi,mj){
		p0.out=0	
		for(k in 1:n){#iterate over the crossover numbers
		    p=pk[k]
	 	    region.length=d/k
	 	    ri=ceiling(mi/region.length)
	 	    rj=ceiling(mj/region.length)
	 	    bi=ri*region.length
	 	    bj=rj*region.length
	 	    	 	    
	 	    if(rj==ri){
	 	      p0=1-((mj-mi)/region.length)
	 	    } else if (rj-ri == 1){
	 	      p0=(1-((bi-mi)/region.length))*((bj-mj)/region.length)
	 	    } else if (rj-ri>1){
	 	      p0=0
	 	    } else {
	 	      print(paste("COULDNT SOLVE",mi,mj,k,p0)) #Trouble shoot
	 	    }
			if(p0<0){print(paste("GOT NEGATIVE p0",mi,mj,k,p0))} #Trouble shoot
			if(mj>d){print(paste("mj>d", mj,d))} #Trouble shoot
			
			p0.out=p0.out+(p*p0) #add the result 
		}
		p0.out=p0.out+p_0 #Add the likelihood of no crossovers
		if(mj-mi==0){p0.out=1}
		return(p0.out)
}
