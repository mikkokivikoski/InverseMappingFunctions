
##INPUT DATA:

 recombinationSitesHelsinki = #Ppungitius_Helsinki_linkagemap.txt.gz
 ninespinedCentromeres = #Centromere positions according to P. pungitius reference genome ver. 7 (Kivikoski et al. 2021)
 v7Index = #Index (.fai) of the P. pungitius reference genome ver. 7 (Kivikoski et al. 2021)
 emAlgorithmFunctions = #emAlgortihmForCoProbabilities.r
  
 p0kFunction = #p0kFunction.r

#STEP 1: Concatenate the recombinations

           IN = read.table(recombinationSitesHelsinki,sep="\t",stringsAsFactors=F,skip=3,header=F)
           header = read.table(recombinationSitesHelsinki,sep="\t",stringsAsFactors=F,nrows=3,header=F)
           table.out=data.frame()
           raw.offspring = data.frame()

           lgs=paste0("LG",1:21)
       for(LG in lgs){
            d=subset(IN,V1==LG)
                        for(i in 6:ncol(d)){
                                MALE = unlist(strsplit(header[1,i],split="f"))[2]
                                FEMALE = paste0(unlist(strsplit(header[1,i],"-"))[1],"-",unlist(strsplit(header[1,i],"-"))[2])
                                FAMILY= header[1,i]
                                OFFSPRING = header[2,i]
                                SEX = header[3,i]

                                tmp = data.frame(d$V1,d$V2,d$V4,d$V5,unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[1]))),unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[2]))))
                                colnames(tmp) = c("CHR","SITE","PATERNALMAP","MATERNALMAP","PATERNAL","MATERNAL")

                                sites.paternal = c()
                                precise.sites.paternal = c()
                                precise.sites.paternal.CM = c()
                                sites.maternal = c()
                                precise.sites.maternal = c()
                                precise.sites.maternal.CM = c()
                                intervals.paternal = c()
                                intervals.maternal = c()

                                first.paternal.non.zero=min(which(tmp$PATERNAL!=0))
                                first.maternal.non.zero=min(which(tmp$MATERNAL!=0))

                                previous.site.MATERNAL = tmp[first.maternal.non.zero,"SITE"]
                                previous.site.MATERNAL.CM = tmp[first.maternal.non.zero,"MATERNALMAP"]
                                previous.site.PATERNAL = tmp[first.paternal.non.zero,"SITE"]
                                previous.site.PATERNAL.CM = tmp[first.paternal.non.zero,"PATERNALMAP"]
                                previous.gt.PATERNAL = tmp[first.paternal.non.zero,"PATERNAL"]
                                previous.gt.MATERNAL = tmp[first.maternal.non.zero,"MATERNAL"]
                                haplotypes.paternal=c(previous.gt.PATERNAL)
                                haplotypes.maternal=c(previous.gt.MATERNAL)

                                for (ROW in 1:nrow(tmp)){
                                  paternal.gt = tmp[ROW,"PATERNAL"]
                                  maternal.gt = tmp[ROW,"MATERNAL"]
                                  if(ROW>=first.paternal.non.zero){
                                          #PATERNAL                         
                                          if (paternal.gt == previous.gt.PATERNAL){
                                                  previous.site.PATERNAL = tmp[ROW,"SITE"]
                                                  previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]

                                          } else if (paternal.gt != previous.gt.PATERNAL & paternal.gt != 0){
                                             precise.sites.paternal = c(precise.sites.paternal, floor(mean(c(previous.site.PATERNAL,tmp[ROW,"SITE"]))))
                                             precise.sites.paternal.CM = c(precise.sites.paternal.CM, mean(c(previous.site.PATERNAL.CM,tmp[ROW,"PATERNALMAP"])))
                                     previous.gt.PATERNAL = paternal.gt
                                                 haplotypes.paternal=c(haplotypes.paternal,previous.gt.PATERNAL)
                                                 previous.site.PATERNAL = tmp[ROW,"SITE"]
                                                 previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]
                                          }
                                  }
                              if(ROW>=first.maternal.non.zero){
                                      #MATERNAL   
                                          if (maternal.gt == previous.gt.MATERNAL){
                                                  previous.site.MATERNAL = tmp[ROW,"SITE"]
                                                  previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]

                                          } else if (maternal.gt != previous.gt.MATERNAL & maternal.gt != 0){
                                             precise.sites.maternal = c(precise.sites.maternal, floor(mean(c(previous.site.MATERNAL,tmp[ROW,"SITE"]))))
                                             precise.sites.maternal.CM = c(precise.sites.maternal.CM, mean(c(previous.site.MATERNAL.CM,tmp[ROW,"MATERNALMAP"])))
                                     previous.gt.MATERNAL = maternal.gt
                                     haplotypes.maternal=c(haplotypes.maternal,previous.gt.MATERNAL)
                                                 previous.site.MATERNAL = tmp[ROW,"SITE"]
                                                 previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]
                                          }
                                  }
                    } #FOR INDIVIDUAL ENDS

                    raw.offspring.tmp = data.frame(OFFSPRING=OFFSPRING, MALE=MALE, FEMALE=FEMALE,CHR=LG,
                    MATERNALCOUNT=length(precise.sites.maternal),PATERNALCOUNT=length(precise.sites.paternal),MATERNALSITES=paste(precise.sites.maternal,collapse=","),
                            SEX=SEX,PATERNALSITES=paste(precise.sites.paternal,collapse=","), PATERNALSITESCM=paste(precise.sites.paternal.CM,collapse=","), MATERNALSITESCM=paste(precise.sites.maternal.CM,collapse=","),
                            PATERNALHAPLOTYPES=paste(haplotypes.paternal,collapse=","),MATERNALHAPLOTYPES=paste(haplotypes.maternal,collapse=","))
                            raw.offspring=rbind(raw.offspring,raw.offspring.tmp)

               }#FOR DATA ENDS
       } #FOR LG ENDS

  max.co.count=max(c(raw.offspring$MATERNALCOUNT,raw.offspring$PATERNALCOUNT))
  for(i in 1:max.co.count){
          raw.offspring[[paste0("MATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITES),split=","),function(x) x[i])
          raw.offspring[[paste0("MATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITESCM),split=","),function(x) x[i])
          raw.offspring[[paste0("PATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITES),split=","),function(x) x[i])
          raw.offspring[[paste0("PATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITESCM),split=","),function(x) x[i])
  }
  table.out=raw.offspring
  
### STEP 1 ends

#### STEP 2: Run the EM algorithm
	        table1 = table.out #From Step 1
	    
	        source(emAlgorithmFunctions)
		    table.out=data.frame()
		    OUT.ARRAY=list()
		     
		    lgs = paste0("LG",seq(1,21))
		    lgs=lgs[1:21]
		    table1=subset(table1, !(OFFSPRING %in% c("Tank1-344", "Tank2-134", "Tank2-31","Tank2-452"))) #Remove the 4 single-offspring families, number of offspring from 938 -> 934.
		    
	  	    maximum.co.count=max(c(table1$MATERNALCOUNT,table1$PATERNALCOUNT))
		    OUT.PATERNAL=data.frame(CHR=lgs)
		    OUT.MATERNAL=data.frame(CHR=lgs)
		    for(i in 0:maximum.co.count){
		        OUT.PATERNAL[[paste0("n",i)]]=100
		        OUT.MATERNAL[[paste0("n",i)]]=100
		    }
		    for(i in 0:((2*maximum.co.count-1))){
		        OUT.PATERNAL[[paste0("MLEp",i)]]=100
		        OUT.MATERNAL[[paste0("MLEp",i)]]=100
		    }
		    
		    for(i in 0:((2*maximum.co.count-1))){
		        OUT.PATERNAL[[paste0("MLEpRestricted",i)]]=100
		        OUT.MATERNAL[[paste0("MLEpRestricted",i)]]=100
		    }
		    
		    BootstrapPvalue=100
		    LowerBound=100
		    UpperBound=100
		    OUT.PATERNAL=cbind(OUT.PATERNAL,BootstrapPvalue,LowerBound,UpperBound)
		    OUT.MATERNAL=cbind(OUT.MATERNAL,BootstrapPvalue,LowerBound,UpperBound)
		    
			ITER=100
			ITER.CONF=100
	 		r=100
	 		ALPHA=0.025		    
	   	    PERFAMILY=FALSE
	   	    NORMAL.LENGTH=T
	        
	        if(NORMAL.LENGTH){
	         ITER=5000
	         ITER.CONF=5000
			 r=1000
	 		}
	
		    main=function(DATA,ITER,r,ALPHA){
		        M=sum(DATA)
		        if(length(DATA)<maximum.co.count+1){
	 		        DATA.print=c(DATA,rep(0,(maximum.co.count+1)-length(DATA)))
	 		    } else{
	 		        DATA.print=DATA
				}
	 		    p.non.restricted = MLEforP(DATA,ITER,FALSE,FALSE,TRUE,prepare.prior(DATA,TRUE))
	 		    q.non.restricted = generate.Qp(DATA,p.non.restricted)
	 		    Generated.data.sets = generateDataFrom.Qp(M,r,q.non.restricted,DATA)
	 		    All.MLEs = MLE.From.Genrated.Data(Generated.data.sets,ITER.CONF)
	            conf.intervals = Find.Conf.Intervals(All.MLEs,ALPHA)
	 		    
	 		    		    
	 		    p.restricted = MLEforP(DATA,ITER,FALSE,TRUE,FALSE,3)
	 		    q.restricted = generate.Qp(DATA,p.restricted)    
	            generated.restricted.datasets = generateDataFrom.Qp(M,r,q.restricted,DATA)
	            
	            L.non.restricted=calculate.likelihood.Function(q.non.restricted,DATA,TRUE)
	            L.restricted=calculate.likelihood.Function(q.restricted,DATA,TRUE)
	            delta.N = 2*L.non.restricted - 2*L.restricted
	            delta.N.simulated = calculate.delta.n(delta.N,generated.restricted.datasets,q.restricted,q.non.restricted,TRUE)
	            P.val=calculate.p.value(delta.N,delta.N.simulated)
	            
	 		    
	 		    if(length(p.non.restricted)<(2*maximum.co.count)){
	 		        OUT.VECTOR=c(DATA.print,round(p.non.restricted,8),rep(0,(2*maximum.co.count)-length(p.non.restricted)),round(p.restricted,8),rep(0,(2*maximum.co.count)-length(p.restricted)),P.val,conf.intervals)
	 		    } else {
			        OUT.VECTOR=c(DATA.print,round(p.non.restricted[1:(2*maximum.co.count)],8),round(p.restricted[1:(2*maximum.co.count)],8),P.val,conf.intervals)
			    }
	 		    
	 		    
	           return(OUT.VECTOR)
		       	    
		    }
		    
		    for(LG in lgs){
		        print(LG)
		        perOffspring=subset(table1,CHR==LG)
	 		    
			    if(TRUE){
		 		    maternal.data=concatenate.data(perOffspring$MATERNALCOUNT)
		 		    paternal.data=concatenate.data(perOffspring$PATERNALCOUNT)
		 		    OUT.PATERNAL[OUT.PATERNAL$CHR==LG,2:(1+1+maximum.co.count+2*maximum.co.count+2*maximum.co.count+3)]=main(paternal.data,ITER,r,ALPHA)
		 		    OUT.MATERNAL[OUT.MATERNAL$CHR==LG,2:(1+1+maximum.co.count+2*maximum.co.count+2*maximum.co.count+3)]=main(maternal.data,ITER,r,ALPHA)
		 		    
		 	    }
		    }
	 		
	 		OUT.ARRAY[["PATERNAL"]]=OUT.PATERNAL
	 		OUT.ARRAY[["MATERNAL"]]=OUT.MATERNAL
			
			colnames(OUT.PATERNAL)=paste0("PATERNAL",colnames(OUT.PATERNAL))
			colnames(OUT.MATERNAL)=paste0("MATERNAL",colnames(OUT.MATERNAL))


# STEP 2 ends

# STEP 3 Predict the recombination frequencies
    

        library(ggplot2)
        table1=table.out #From Step 1
		table2=metadata #Chromosome lengths
		IN=OUT.ARRAY #From Step 2
        source(p0kFunction)
		
	    table.out=data.frame()
		table1=subset(table1, !(OFFSPRING %in% c("Tank1-344", "Tank2-134", "Tank2-31","Tank2-452"))) #Remove the 4 single-offspring families, number of offspring from 938 -> 934.
		####
		haldane=function(x){0.5*(1-exp((-1)*0.02*x))}
		kosambi=function(x){(exp(0.04*x)-1)/(2*(exp(0.04*x)+1))}
		linear=function(x){x/100}
		linear2=function(x){min(c(0.5,x/100))}
		x.values=seq(0,200,0.01)
		mapping.functions=data.frame(x=x.values,y=c(haldane(x.values),kosambi(x.values),linear(x.values)),mapping.function=rep(c("Haldane","Kosambi","Linear"),each=length(x.values)))
		crossover.events=table1
		metadata=table2
		
		sum.of.squares=data.frame()
        rbar.out=data.frame()

		for(lg in 1:21){
		cm.length.maternal=metadata[lg,"CHRLENGTHCMMATERNAL"]
		t = read.table(gzfile(paste0("chr", lg, ".rf.gz")), header=F) #Empirical recombination frequencies and map distances based on the linkage maps
		t=subset(t,V1>V2)
		
       ## t=t[sort(sample(1:nrow(t),5000,replace=F)),] #for quicker execution
         
        rbar.loci.male=nrow(subset(t, V4=="male"))
		rbar.loci.female=nrow(subset(t, V4=="female"))
		
		t$V4=factor(t$V4,levels=c("female","male"),labels=c("Maternal","Paternal"))
		t$recomb=t$V3/100
		t$CHR=paste0("LG",lg)
		
		
		p=ggplot(t,aes(x=V1-V2,y=recomb)) + 
			geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Maternal"),size=0.5,shape=21,stroke=0) +
			geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Paternal"),size=0.5,shape=21,stroke=0) + 
			geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
			ylim(0,0.60)+xlim(0,cm.length.maternal)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
			scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3","firebrick")) + 
			scale_fill_manual(values=c("tan1","steelblue1")) + 
			labs(fill="Empirical data", colour="Mapping function",title=paste0("LG",lg))  +    
			guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
			theme(legend.position="right",aspect.ratio=1)#+facet_wrap(~CONTEXT,ncol=6)
			ggsave(paste0("lg",lg,"empirical.png"),device="png",width=100,height=60,units="mm",dpi=100)
			p=p+facet_wrap(~V4,ncol=2)
            ggsave(paste0("lg",lg,"empirical_2panels.png"),device="png",width=100,height=60,units="mm",dpi=100)
			
			p=ggplot(t,aes(x=V1-V2,y=recomb)) + 
			geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Maternal"),size=0.5,shape=21,stroke=0) +
			geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
			ylim(0,0.60)+xlim(0,cm.length.maternal)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
			scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3","firebrick")) + 
			scale_fill_manual(values=c("tan1","steelblue1")) + 
			labs(fill="Empirical data", colour="Mapping function",title=paste0("LG",lg))  +    
			guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
			theme(legend.position="right",aspect.ratio=1)#+facet_wrap(~CONTEXT,ncol=6)
			ggsave(paste0("lg",lg,"empirical_female.png"),device="png",width=100,height=60,units="mm",dpi=300)
		
			p=ggplot(t,aes(x=V1-V2,y=recomb)) + 
			geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Paternal"),size=0.5,shape=21,stroke=0) +
			geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
			ylim(0,0.60)+xlim(0,cm.length.maternal)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
			scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3","firebrick")) + 
			scale_fill_manual(values=c("steelblue1")) + 
			labs(fill="Empirical data", colour="Mapping function",title=paste0("LG",lg))  +    
			guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
			theme(legend.position="right",aspect.ratio=1)#+facet_wrap(~CONTEXT,ncol=6)
			ggsave(paste0("lg",lg,"empirical_male.png"),device="png",width=100,height=60,units="mm",dpi=300)
		
			
			plot.table1.rbar=t
            plot.table1.rbar[plot.table1.rbar$recomb >0.5, "recomb"]=1-plot.table1.rbar[plot.table1.rbar$recomb >0.5, "recomb"]
          
            rbar.tmp=data.frame(Chr=lg,Sex=c("Maternal","Paternal"), Nrbarloci=c(rbar.loci.female,rbar.loci.male), SumRij=as.numeric(aggregate(recomb ~ V4, sum, data=plot.table1.rbar)$recomb),Source="Empirical")
           
			
	    Kosambi=kosambi(t$V1-t$V2)
	    Haldane=haldane(t$V1-t$V2)
	    Linear=sapply(t$V1-t$V2,function(x) linear2(x))
	    
	    P0K=c()
	    plot.p0k=data.frame()

	       for(sex in c("Paternal","Maternal")){
		       n=max(which(metadata[metadata$CHR==paste0("LG",lg),paste0(casefold(sex,upper=T),"MLEp",1:9)]>0))
		       pk=as.numeric(metadata[metadata$CHR==paste0("LG",lg),paste0(casefold(sex,upper=T),"MLEp",1:9)])
		       p_0=as.numeric(metadata[metadata$CHR==paste0("LG",lg),paste0(casefold(sex,upper=T),"MLEp",0)])
               tmp.meta=IN[[casefold(sex,upper=T)]]
		       if(as.numeric(tmp.meta[tmp.meta$CHR==paste0("LG",lg),"BootstrapPvalue"])>=0.05){
		               pk=as.numeric(tmp.meta[tmp.meta$CHR==paste0("LG",lg),grepl("MLEpRestricted",colnames(tmp.meta))])
		               p_0=pk[1]
		               pk=pk[-1]
			   }
		       
		       tmp=subset(t,V4==sex)
		       length.cm=max(tmp$V1,tmp$V2,as.numeric(metadata[metadata$CHR==paste0("LG",lg),paste0("CHRLENGTHCM",casefold(sex,upper=T))]))
		            
			   p0k.prediction=as.numeric(apply(tmp,1,function(x) p0k(pk=pk,p_0=p_0,n=n,d=length.cm,mi=as.numeric(x[2]),mj=as.numeric(x[1]))))
			   #Estimate with two crossover classes
			   max2=floor(length.cm/50)
			   p_max=(length.cm-(max2*50))/length.cm

		       if(sex=="Paternal"){COLOR=rgb(0,0,1,0.05)}else{COLOR=rgb(1,0.65,0.31,0.05)}
		       plot.p0k=rbind(plot.p0k,data.frame(x=(tmp$V1-tmp$V2),y=(0.5*(1-p0k.prediction)),SEX=sex))		       
		       
	      P0K=c(P0K,(0.5*(1-p0k.prediction)))
	     }
	     
	     
        plot.p0k$SEX=factor(plot.p0k$SEX, levels=c("Maternal","Paternal")) 
			
		p=ggplot(plot.p0k,aes(x=x,y=y)) + 
		geom_point(aes(fill=SEX),alpha=0.05,data=subset(plot.p0k,SEX=="Maternal"),size=0.5,shape=21,stroke=0) +
		#geom_point(aes(fill=SEX),alpha=0.05,data=subset(plot.p0k,SEX=="Paternal"),size=0.5,shape=21,stroke=0) + 
		geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		scale_fill_manual(values=c("tan1","steelblue1")) + 
		labs(fill="Prediction", colour="Mapping function",title=paste0("LG",lg))  +    
		guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
		theme(legend.position="right",aspect.ratio=1)
		ggsave(paste0("lg",lg,"prediction_female.png"),device="png",width=100,height=60,units="mm",dpi=300)
		
		p=ggplot(plot.p0k,aes(x=x,y=y)) + 
		#geom_point(aes(fill=SEX),alpha=0.05,data=subset(plot.p0k,SEX=="Maternal"),size=0.5,shape=21,stroke=0) +
		geom_point(aes(fill=SEX),alpha=0.05,data=subset(plot.p0k,SEX=="Paternal"),size=0.5,shape=21,stroke=0) + 
		geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		scale_fill_manual(values=c("steelblue1")) + 
		labs(fill="Prediction", colour="Mapping function",title=paste0("LG",lg))  +    
		guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
		theme(legend.position="right",aspect.ratio=1)
		ggsave(paste0("lg",lg,"prediction_male.png"),device="png",width=100,height=60,units="mm",dpi=300)
		
	    t=rbind(t,t,t,t)
	    t$predicted=c(Kosambi,Haldane,Linear,P0K)
	    t$mappingFunction=rep(c("Kosambi","Haldane","Linear","P0K"),each=length(Kosambi))
	    t$mappingFunction=factor(t$mappingFunction,labels=c("Haldane","Kosambi","Linear","p0(k)"))
	    print(head(t))
	   
	    t$loglikelihood=as.numeric(apply(t,1,function(x) dbinom(x=round(as.numeric(x[5])*934), size=934, prob=as.numeric(x[7]),log=T )))
        t$loglikelihood=log(as.numeric(apply(t,1,function(x) dbinom(x=round(as.numeric(x[5])*934), size=934, prob=as.numeric(x[7]))))+10^(-9))#add small number to avoid log(0)
                  
        rbar.tmp.function=aggregate(predicted ~ V4 + mappingFunction, sum, data=t)
        rbar.tmp.function=data.frame(Chr=lg,Sex=rbar.tmp.function$V4,Nrbarloci=c(rbar.loci.female,rbar.loci.male), SumRij=as.numeric(rbar.tmp.function$predicted),Source=rbar.tmp.function$mappingFunction)
        
        rbar.out=rbind(rbar.out,rbar.tmp,rbar.tmp.function)
        
        t$sq.error=(t$recomb-t$predicted)^2
	    t$error=(t$recomb-t$predicted)
	    t$error.absolute=abs(t$error)
       
	    p1=ggplot(subset(t,mappingFunction!="p0(k)*"),aes(x=V1-V2,y=recomb-predicted)) + 
	    geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Maternal" & mappingFunction!="p0(k)*"),size=0.5,shape=21,stroke=0) +
	    geom_point(aes(fill=V4),alpha=0.05,data=subset(t,V4=="Paternal" & mappingFunction!="p0(k)*"),size=0.5,shape=21,stroke=0) + 
	    geom_hline(yintercept=0,linetype=2) +
	    xlim(0,170)+ xlab("Map distance (cM)")+ylab("Error")+
	    scale_fill_manual(values=c("tan1","steelblue1")) + 	
	    labs(fill="Sex",title=paste0("LG",lg))  +    
	    guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
	    theme(legend.position="none",text=element_text(size=5))+
	    facet_grid(V4~mappingFunction,scales="free_y")
	    ggsave(paste0("lg",lg,"residual_scatter.png"),device="png",width=100,height=60,units="mm",dpi=250)
		
		

		
	    t2=t
	    tmp.sum.of.squares=aggregate(t2$sq.error,by=list(t2$V4,t2$mappingFunction,t2$CHR),sum)
	    colnames(tmp.sum.of.squares)=c("V4","mappingFunction","chromosome","sq.sum")
	    tmp.sum.of.squares$sq.mean=aggregate(t2$sq.error,by=list(t2$V4,t2$mappingFunction,t2$CHR),mean)$x
	    tmp.sum.of.squares$sq.median=aggregate(t2$sq.error,by=list(t2$V4,t2$mappingFunction,t2$CHR),median)$x
	    tmp.sum.of.squares$sq.sd=aggregate(t2$sq.error,by=list(t2$V4,t2$mappingFunction,t2$CHR),sd)$x
	    tmp.sum.of.squares$logLikelihood=aggregate(t2$loglikelihood,by=list(t2$V4,t2$mappingFunction,t2$CHR),sum)$x
	    tmp.sum.of.squares$TotalError=aggregate(error ~ V4 + mappingFunction + CHR,sum,data=t2)$error
	    tmp.sum.of.squares$TotalErrorAbsolute=aggregate(error.absolute ~ V4 + mappingFunction + CHR,sum,data=t2)$error.absolute
            tmp.sum.of.squares$TotalErrorLow95=aggregate(error ~ V4 + mappingFunction + CHR, function(x) quantile(x,c(0.025,0.975))[1], data=t2)$error
	    tmp.sum.of.squares$TotalErrorUpper95=aggregate(error ~ V4 + mappingFunction + CHR, function(x) quantile(x,c(0.025,0.975))[2], data=t2)$error
	    tmp.sum.of.squares$TotalErrorSd=aggregate(error ~ V4 + mappingFunction + CHR, function(x) sd(x), data=t2)$error
	   
		    
		n.markers=nrow(t)/4 #4 for the 4 inverse mapping functions 
		tmp.sum.of.squares$NumberOfMarkers=n.markers
	    
		sum.of.squares=rbind(sum.of.squares,tmp.sum.of.squares)
		
	 }
	 table.out=sum.of.squares
	 
     OUT=list()
	 OUT[["rbar"]]=rbar.out
	 
#### STEP 3 ends


### STEP 4 Collect the results
	 library(ggplot2)

     table1=table.out #from STEP 3
     table.out=data.frame()
	 ####
	 OUT=list()
	 sum.of.squares=table1
	 table.out=sum.of.squares
	 sum.of.squares.total=aggregate(sq.sum ~ V4 + mappingFunction, sum, data=sum.of.squares)

	 sum.of.squares$chromosome=factor(sum.of.squares$chromosome,levels=paste0("LG",1:21)) 
	 
	 
	 best=aggregate(logLikelihood ~ chromosome+V4,data=sum.of.squares,max)
	 best=subset(sum.of.squares,logLikelihood %in% best$logLikelihood)
	 best$LABEL="*"
	 best[best$mappingFunction=="Linear","LABEL"]="L"
	 best[best$mappingFunction=="Kosambi","LABEL"]="K"
	 best[best$mappingFunction=="Haldane","LABEL"]="H"
	 best[best$mappingFunction=="p0(k)","LABEL"]="p0(k)"
	 best[best$mappingFunction=="p0(k)*","LABEL"]="p0(k)*"
	 sum.of.squares$mappingFunction=factor(sum.of.squares$mappingFunction,labels=c("Haldane", "Kosambi", "Linear", "p0(k)", "p0(k)*"))
	 ggplot(sum.of.squares,aes(x=chromosome,y=logLikelihood,fill=mappingFunction))+geom_jitter(shape=23,stroke=0.1,alpha=0.75,size=1,width=0.05,height = 0)+
	 scale_color_manual(values=c("#F0E442","#009E73","orchid","tan2","firebrick")) + 	
	 scale_fill_manual(values=c("#F0E442","#009E73","orchid","tan2","firebrick")) +
	 ylab("log-ikelihood of data")+xlab("Chromosome")+labs(fill="Function") + 
	 geom_text(data=best,aes(x=chromosome,y=0,label=LABEL),size=0.7,inherit.aes=F)+
	 expand_limits(y=0) +
	 facet_wrap(~V4,ncol=2,scales="free_y")+
	 guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))+
	 theme(text=element_text(size=7),axis.text.x=element_text(angle=60, hjust = 1),panel.grid.minor=element_blank(),legend.position="bottom")
	 ggsave("SumOfLikelihoodsOfData.png",device="png",width=110,height=80,units="mm",dpi=1000)
	 ggsave("SumOfLikelihoodsOfData.pdf",device="pdf",width=160,height=60,units="mm")
	
	likelihood_total=aggregate(sum.of.squares$logLikelihood,by=list(sum.of.squares$V4,sum.of.squares$mappingFunction),sum)
	  colnames(likelihood_total)=c("Sex","mappingFunction","Sum")
	  
	  
	  
	  
	  plot.table.error=aggregate(cbind(TotalError, NumberOfMarkers) ~ V4 + mappingFunction, sum, data=table1)
	  plot.table.error$MeanError=plot.table.error$TotalError/plot.table.error$NumberOfMarkers
	  plot.table.error$MeanAbsoluteError=aggregate(TotalErrorAbsolute~ V4 + mappingFunction, sum, data=table1)$TotalErrorAbsolute/plot.table.error$NumberOfMarkers
	  plot.table.error$TotalErrorSd=aggregate(TotalErrorSd ~ V4 + mappingFunction, mean, data=table1)$TotalErrorSd
	  
	  ggplot(plot.table.error,aes(x=mappingFunction,y=MeanError,ymin=MeanError-TotalErrorSd,ymax=MeanError+TotalErrorSd,fill=mappingFunction))+
	  geom_errorbar(width=0.2) + geom_point(shape=23) +
	  scale_color_manual(values=c("#F0E442","#009E73","orchid","tan2","firebrick")) + 	
	  scale_fill_manual(values=c("#F0E442","#009E73","orchid","tan2","firebrick")) +
	  ylab("Error")+xlab("Mapping function")+labs(fill="Mapping function")  + 
	  facet_wrap(~V4,ncol=1,scales="free_y")+
	  theme(legend.position="bottom",text=element_text(size=6))
	  ggsave("TotalError.png",device="png",width=110,height=80,units="mm",dpi=1000)
	 
	 
	  
	  OUT[["SumOfSquares"]]=sum.of.squares.total
	  OUT[["likelihood"]]= likelihood_total
	  OUT[["Error"]]=plot.table.error
	  

#STEP 4 ends

#STEP 5: Calculate rbar according to (Veller et al. 2019)
      table1=metadata #Chromosome lengths
      IN=OUT #From step 3
      table.out=data.frame()
	  #############
	
		rbar=IN[["rbar"]]
	    rbar=subset(rbar, Chr !=12) #Remove the sex chromosome
	    rbar$IntraContrib=rbar$SumRij/rbar$Nrbarloci #This gives the intra-chromosomal component of one chromosome
	    
	    rbar$ChrLength=rep(table1$LENGTH[-12], each=12)
	    rbar$IntraContrib2=rbar$IntraContrib*((rbar$ChrLength/sum(as.numeric(table1$LENGTH[-12])))^2) #This gives the weighted component
	    print(aggregate(IntraContrib ~ Sex + Source, sum,data=rbar))
	    print(aggregate(IntraContrib2 ~ Sex + Source, sum,data=rbar))
	    inter.chrom=data.frame(Species="Ppun",InterChromComponent=0.5*(1-sum((table1$LENGTH[-12]/sum(table1$LENGTH[-12]))^2)))
	    OUT=list()
	    OUT[["rbar"]]=rbar
	    OUT[["interChrom"]]=inter.chrom
	    table.out=aggregate(IntraContrib ~ Sex + Source, sum,data=rbar)
	    table.out$IntraContrib2=aggregate(IntraContrib2 ~ Sex + Source, sum,data=rbar)$IntraContrib2
	    total.rbar=aggregate(IntraContrib2 ~ Source, sum,data=table.out)
	    total.rbar$varIBD=0.125*(1-total.rbar$IntraContrib2-2*inter.chrom$InterChromComponent)  
	    OUT[["ibdVarRbar"]]=total.rbar
	    




