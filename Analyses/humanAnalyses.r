

  chromosomeLengths = read.table("")#table of bp lengths of human autosomes
  humanCO = read.table("")# Data S4 of Halldorsson et al. (2019)
  paternalMap = read.table("")# Data S1 of Halldorsson et al. 2019
  maternalMap = read.table("")# Data S2 of Halldorsson et al. 2019
  p0kFunction="p0kFunction.r" #Source for the p_0(k) function
  emAlgorithmFunctions="emAlgortihmForCoProbabilities.r"

# STEP 1: Process crossovers

    table1=humanCO
	probands.m=unique(subset(table1,ptype=="M")$pid)
	probands.p=unique(subset(table1,ptype=="P")$pid)
	probands=intersect(probands.m,probands.p)
	probands.for.co.count=probands
	probands.for.rec.count=probands
	
	out=data.frame()
	table1=subset(table1,pid %in% probands)
	table1$count=1
	##ADD A PSEUDO CROSSOVER TO LIST ALL PROBANDS IN ALL CHROMOSOMES
	##THE BASE PAIR LOCATION OF THE PSEUDO CROSSOVER IS 10^9 TO ENSURE IT CAN BE OMITTED.
	table1=rbind(table1,data.frame(Chr=rep(unique(table1$Chr),each=2*length(probands)),lbnd=NA,ubnd=NA,pid=rep(probands,each=2),ptype=c("P","M"),medsnp=10**9,complex="NA",count=0))
	###
	
	maternal=subset(table1,ptype=="M")
	paternal=subset(table1,ptype=="P")
	
	maternal.co=aggregate(maternal$count,by=list(maternal$Chr,maternal$pid),sum)
	maternal.sites=aggregate(maternal$medsnp,by=list(maternal$Chr,maternal$pid),function(x) paste(x,collapse=","))

	paternal.co=aggregate(paternal$count,by=list(paternal$Chr,paternal$pid),sum)
	paternal.sites=aggregate(paternal$medsnp,by=list(paternal$Chr,paternal$pid),function(x) paste(x,collapse=","))

	table.out=data.frame(PROBAND=maternal.co$Group.2,CHR=maternal.co$Group.1,MATERNALCOUNT=maternal.co$x,
	PATERNALCOUNT=paternal.co$x,MATERNALSITES=maternal.sites$x,PATERNALSITES=paternal.sites$x)

	maternal.co.out=subset(maternal.co, Group.2 %in% probands.for.co.count)
	paternal.co.out=subset(paternal.co, Group.2 %in% probands.for.co.count)

	maternal.co.out.sites=subset(maternal.sites, Group.2 %in% probands.for.co.count)
	paternal.co.out.sites=subset(paternal.sites, Group.2 %in% probands.for.co.count)
	
	maternal.rec.out=subset(maternal.co, Group.2 %in% probands.for.rec.count)
	paternal.rec.out=subset(paternal.co, Group.2 %in% probands.for.rec.count)
	maternal.rec.out.sites=subset(maternal.sites, Group.2 %in% probands.for.rec.count)
	paternal.rec.out.sites=subset(paternal.sites, Group.2 %in% probands.for.rec.count)
		
	table.out=data.frame(PROBAND=maternal.co.out$Group.2,CHR=maternal.co.out$Group.1,MATERNALCOUNT=maternal.co.out$x,
	
	PATERNALCOUNT=paternal.co.out$x,MATERNALSITES=maternal.co.out.sites$x,PATERNALSITES=paternal.co.out.sites$x)
    OUT=list()
    OUT[["SITES"]]=data.frame(PROBAND=maternal.rec.out$Group.2,CHR=maternal.rec.out$Group.1,MATERNALCOUNT=maternal.rec.out$x,
	PATERNALCOUNT=paternal.rec.out$x,MATERNALSITES=maternal.rec.out.sites$x,PATERNALSITES=paternal.rec.out.sites$x)
  
  

# STEP 1 ends

# STEP 2: Run the EM algorithm

		    source(emAlgorithmFunctions)
		    table.out=data.frame()
		    table1=read.table("",header=T,stringsAsFactors=F,sep="\t") #table.out from Step 1 
		    chrs = paste0("chr",seq(1,22))
		    chrs=chrs[1:22]
	  	    maximum.co.count=max(c(table1$MATERNALCOUNT,table1$PATERNALCOUNT))
	  	    print(maximum.co.count)
		    OUT.PATERNAL=data.frame(CHR=chrs)
		    OUT.MATERNAL=data.frame(CHR=chrs)
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
		    
			ITER=10
			ITER.CONF=10
	 		r=10
	 		ALPHA=0.025		    
	   	    PERFAMILY=FALSE
	   	    NORMAL.LENGTH=F #Set TRUE for final results
	        
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
		    
		    for(chr in chrs){
		        print(chr)
		        perOffspring=subset(table1,CHR==chrs)
	 		    
			    if(TRUE){
		 		    maternal.data=concatenate.data(perOffspring$MATERNALCOUNT)
		 		    paternal.data=concatenate.data(perOffspring$PATERNALCOUNT)
		 		    OUT.PATERNAL[OUT.PATERNAL$CHR==chr,2:(1+1+maximum.co.count+2*maximum.co.count+2*maximum.co.count+3)]=main(paternal.data,ITER,r,ALPHA)
		 		    OUT.MATERNAL[OUT.MATERNAL$CHR==chr,2:(1+1+maximum.co.count+2*maximum.co.count+2*maximum.co.count+3)]=main(maternal.data,ITER,r,ALPHA)
		 		    
		 	    }
		    }
	 		

# STEP 2 ends

# STEP 3: Plot gamtic and bivalent crossover  counts
 
 	library(ggplot2)
	#################
	COLORS <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    tmp=OUT.MATERNAL #From STEP 2
	
	no.probands=sum(as.numeric(tmp[1,grepl("n[0-9]",colnames(tmp))]))
	plot.table=data.frame(CHR=tmp$CHR,n=as.numeric(unlist(tmp[paste0("n",0:17)])),CLASS=rep(0:17,each=22),RESULT="Observed")
	plot.table=rbind(plot.table,data.frame(CHR=tmp$CHR,n=no.probands*as.numeric(unlist(tmp[paste0("MLEp",0:33)])),CLASS=rep(0:33,each=22),RESULT="Inferred"))
	plot.table$Sex="Maternal"
	
    tmp=OUT.MATERNAL #From STEP 3	
    plot.table=rbind(plot.table,data.frame(CHR=tmp$CHR,n=as.numeric(unlist(tmp[paste0("n",0:17)])),CLASS=rep(0:17,each=22),RESULT="Observed",Sex="Paternal"))
	plot.table=rbind(plot.table,data.frame(CHR=tmp$CHR,n=no.probands*as.numeric(unlist(tmp[paste0("MLEp",0:33)])),CLASS=rep(0:33,each=22),RESULT="Inferred",Sex="Paternal"))
	plot.table$CLASS.2=plot.table$CLASS
	plot.table[plot.table$CLASS.2>17,"CLASS.2"]=-1
	plot.table[plot.table$CLASS.2>5,"CLASS.2"]=-2
	plot.table$CLASS.2=as.character(plot.table$CLASS.2)
	plot.table[plot.table$CLASS.2=="-2","CLASS.2"]="6-17"
	plot.table[plot.table$CLASS.2=="-1","CLASS.2"]="18-33"
	plot.table$CLASS.2=factor(plot.table$CLASS.2,levels=c("0","1","2","3","4","5","6-17","18-33"))
	
	plot.table$CHR=factor(plot.table$CHR,levels=paste0("chr",c(21,22,20:1)))
	plot.table$CHR2=factor(sapply(strsplit(as.character(plot.table$CHR),split="chr"),function (x) x[2]),levels=c(21,22,20:1))
	plot.table$Species="Human"
	
	ggplot(plot.table, aes(x=CHR2,fill=as.factor(CLASS.2)))+geom_bar(aes(weight=n/no.probands))+scale_x_discrete("Chromosome")+
	   	      scale_fill_manual(values=COLORS)+	  
   	      labs(fill = "Crossover count") + ylab(label="Proportion") + scale_y_continuous(breaks=c(0,0.5,1)) +
   		  facet_grid(RESULT~Sex) + guides(fill = guide_legend(nrow = 1)) +
          theme(axis.text.x=element_text(angle=60, hjust = 1,size=8),legend.position="bottom")   		

# STEP 3 ends

# STEP 4: Get the empirical recombination frequencies


	table1=table.out #From STEP 1 
	table2=paternalMap
	table3=maternalMap
	
   recombination.out=data.frame()

	for (chr in paste0("chr",1:22)){
	
       paternal.map=subset(table2,Chr==chr)
       set.seed(2021)
       markers.to.use=sort(sample(1:nrow(paternal.map),round(0.015*nrow(paternal.map)),replace=F)) # SAMPLE 1.5% of the markers
         
         paternal.map=paternal.map[markers.to.use,]
         maternal.map=subset(table3,Chr==chr)
         maternal.map=maternal.map[markers.to.use,]
     
         markers_bp=paternal.map$End 
         markers_cm=paternal.map$cM 

	 	 paternal.markers=data.frame(markerA_bp=rep(markers_bp,each=length(markers_bp)),
		 markerA_cm=rep(markers_cm,each=length(markers_cm)),
		 markerB_bp=markers_bp,
		 markerB_cm=markers_cm,recombination=-1
		 )		
		 paternal.markers=subset(paternal.markers,markerA_bp<markerB_bp)
		 
		 
		 markers_bp=maternal.map$End #c(maternal.map$Begin,max(maternal.map$End))
         markers_cm=maternal.map$cM #c(0,maternal.map$cM)

		maternal.markers=data.frame(markerA_bp=rep(markers_bp,each=length(markers_bp)),
		markerA_cm=rep(markers_cm,each=length(markers_cm)),
		markerB_bp=markers_bp,
		markerB_cm=markers_cm,recombination=-1
		)		
		
		maternal.markers=subset(maternal.markers,markerA_bp<markerB_bp)
	    
	    #Calculate empirical recombination frequency (paternal)
	    rec=as.numeric(apply(paternal.markers,1,function(x) prop.table(table(factor(as.numeric(apply(data,1,function(y) abs(length(which(as.numeric(unlist(strsplit(as.character(y[6]),split=",")))<as.numeric(x[1])))%%2 - length(which(as.numeric(unlist(strsplit(as.character(y[6]),split=",")))<as.numeric(x[3])))%%2))),levels=0:1),useNA="no"))[2]))
	    paternal.markers$recombination=rec	
		paternal.markers$sex="Paternal"
		table.out=paternal.markers
        paternal.markers$Chr=chr
	    recombination.out=rbind(recombination.out,paternal.markers)
        
        #Calculate empirical recombination frequency (maternal)
        rec=as.numeric(apply(maternal.markers,1,function(x) prop.table(table(factor(as.numeric(apply(data,1,function(y) abs(length(which(as.numeric(unlist(strsplit(as.character(y[5]),split=",")))<as.numeric(x[1])))%%2 - length(which(as.numeric(unlist(strsplit(as.character(y[5]),split=",")))<as.numeric(x[3])))%%2))),levels=0:1),useNA="no"))[2]))
        maternal.markers$recombination=rec	
		maternal.markers$sex="Maternal"
		maternal.markers$Chr=chr
		recombination.out=rbind(recombination.out,maternal.markers)

	}	
	recombination.out #output table with empirical recombination frequencies

#STEP 4 ends

#STEP 5: Predict recombination frequencies and compare to empirical frequencies

   
        table.out=data.frame()
        
        OUT=list()
		OUT[["sumOfSquares"]]=data.frame()
		OUT[["rbar"]]=data.frame()
		
	    library(ggplot2)	
		source(p0kFunction)
        haldane=function(x){0.5*(1-exp((-1)*0.02*x))}
		kosambi=function(x){(exp(0.04*x)-1)/(2*(exp(0.04*x)+1))}
		linear2=function(x){min(c(0.5,x/100))}
		x.values=seq(0,350,0.01)
		mapping.functions=data.frame(x=x.values,y=c(haldane(x.values),kosambi(x.values),linear(x.values)),mapping.function=rep(c("Haldane","Kosambi","Linear"),each=length(x.values)))
	  
	   
	    for(chr in paste0("chr",1:22)){
		    table1=subset(recombination.out,Chr==chr)#Output from STEP 4
	   	    rbar.loci.male=nrow(subset(table1, sex=="Paternal"))
	        rbar.loci.female=nrow(subset(table1, sex=="Maternal"))
	        print(paste0(chr, " For rbar, the number of pairs: ", nrow(table1)))
	        n.probands=no.probands #See Step 1
	        
	        t=table1
	        Kosambi=kosambi(t$markerB_cm-t$markerA_cm)
		    Haldane=haldane(t$markerB_cm-t$markerA_cm)
		    Linear=sapply(t$markerB_cm-t$markerA_cm,function(x) linear2(x))
		    P0K=c()
		    plot.p0k=data.frame()   
		    
		    for(sex_tmp in c("Paternal","Maternal")){
		           
		           if(sex_tmp=="Maternal"){metadata=co.likelihoods.maternal}
		           if(sex_tmp=="Paternal"){metadata=co.likelihoods.paternal}
		           
		           pk=as.numeric(metadata[metadata$CHR==chr,grepl("MLEp[0-9]",colnames(metadata))])
		           n=max(which(pk>0))-1
		           p_0=pk[1]
		           pk=pk[-1]
		           if(as.numeric(metadata[metadata$CHR==chr,"BootstrapPvalue"])>=0.05){
		               pk=as.numeric(metadata[metadata$CHR==chr,grepl("MLEpRestricted",colnames(metadata))])
		               p_0=pk[1]
		               pk=pk[-1]
			       }
		           tmp=subset(t,sex==sex_tmp)
			       if(sex_tmp=="Maternal"){length.cm=max(tmp$markerA_cm,tmp$markerB_cm)+0.001}
			       if(sex_tmp=="Paternal"){length.cm=max(tmp$markerA_cm,tmp$markerB_cm)+0.001}
			       print(paste(chr,sex_tmp,length.cm))
			       tmp$markerA_cm=round(tmp$markerA_cm,5) #Round to avoid numerical problems			       
			       tmp$markerB_cm=round(tmp$markerB_cm,5) #Round to avoid numerical problems		       
			       p0k.prediction=as.numeric(apply(tmp,1,function(x) p0k(pk=pk,p_0=p_0,n=n,d=length.cm,mi=as.numeric(x[2]),mj=as.numeric(x[4]))))
		           P0K=c(P0K,(0.5*(1-p0k.prediction)))
		     }
		     
		     plot.table1=rbind(subset(t,sex=="Maternal"),subset(t,sex=="Paternal"),make.row.names=F)
      
		p1=ggplot(plot.table1,aes(x=(markerB_cm-markerA_cm), y=recombination))+
            geom_point(aes(fill=sex),alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		    scale_fill_manual(values=c("tan1","steelblue1")) + 
		    labs(fill="Empirical data", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"empirical.png"),device="png",width=100,height=60,units="mm")
            
            p1=p1+facet_wrap(~sex,ncol=2)
            ggsave(paste0(chr,"empirical_2panels.png"),device="png",width=100,height=60,units="mm")
            
            p1=ggplot(subset(plot.table1,sex=="Maternal"),aes(x=(markerB_cm-markerA_cm), y=recombination))+
            geom_point(fill="tan1",alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		    labs(fill="Empirical data", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"empirical_female.png"),device="png",width=100,height=60,units="mm")

            p1=ggplot(subset(plot.table1,sex=="Paternal"),aes(x=(markerB_cm-markerA_cm), y=recombination))+
            geom_point(fill="steelblue1",alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		    labs(fill="Empirical data", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"empirical_male.png"),device="png",width=100,height=60,units="mm")
            
		    
		    
		    plot.table1.rbar=plot.table1
            plot.table1.rbar[plot.table1.rbar$recombination >0.5, "recombination"]=1-plot.table1.rbar[plot.table1.rbar$recombination >0.5, "recombination"]
        
            rbar.tmp=data.frame(Chr=chr,Sex=c("Maternal","Paternal"), Nrbarloci=c(rbar.loci.female,rbar.loci.male), SumRij=as.numeric(aggregate(recombination ~ sex, sum, data=plot.table1.rbar)$recombination),Source="Empirical")
           
           
            ##############
            plot.table=rbind(t,t,t,t) #Four for the four inverse mapping functions
            
		    plot.table$prediction=c(Haldane,Kosambi,Linear,P0K)
            plot.table$MappingFunction=rep(c("Haldane","Kosambi","Linear","P0K"),each=nrow(t))
            plot.table$MappingFunction=factor(plot.table$MappingFunction,labels=c("Haldane","Kosambi","Linear","p0(k)"))
         	rbar.tmp.function=aggregate(prediction ~ sex + MappingFunction, sum, data=plot.table)
            rbar.tmp.function=data.frame(Chr=chr,Sex=rbar.tmp.function$sex,Nrbarloci=c(rbar.loci.female,rbar.loci.male), SumRij=as.numeric(rbar.tmp.function$prediction),Source=rbar.tmp.function$MappingFunction)
            
		    
		    plot.table$loglikelihood=as.numeric(apply(plot.table,1,function(x) dbinom(x=round(as.numeric(x[5])*n.probands), size=n.probands, prob=as.numeric(x[8]),log=T )))
            plot.table$loglikelihood=log(as.numeric(apply(plot.table,1,function(x) dbinom(x=round(as.numeric(x[5])*n.probands), size=n.probands, prob=as.numeric(x[8]))))+10^(-9))#add small number to avoid log(0)
            
            plot.table$sq.error=(plot.table$recombination-plot.table$prediction)^2
            plot.table$error=(plot.table$recombination-plot.table$prediction)
            plot.table$error.absolute=abs(plot.table$error)
         	plot.table$p0=1-2*plot.table$prediction
         	plot.table$binom.test.pvalue=as.numeric(apply(plot.table,1,function(x) binom.test(x=round(as.numeric(x[5])*n.probands), n=n.probands, p = max(0,as.numeric(x[8])))$p.value))
 		    
		    
            plot.table=rbind(subset(plot.table,sex=="Maternal"),subset(plot.table,sex=="Paternal"),make.row.names=F)
            p1=ggplot(subset(plot.table,MappingFunction=="p0(k)"),aes(x=(markerB_cm-markerA_cm), y=prediction))+
            geom_point(aes(fill=sex),alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		    scale_fill_manual(values=c("tan1","steelblue1")) + 
		    labs(fill="Prediction", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"prediction.png"),device="png",width=100,height=60,units="mm")
            p1=p1+facet_wrap(~sex,ncol=2)
            ggsave(paste0(chr,"prediction_2panels.png"),device="png",width=100,height=60,units="mm")

            p1=ggplot(subset(plot.table,MappingFunction=="p0(k)" & sex=="Maternal"),aes(x=(markerB_cm-markerA_cm), y=prediction))+
            geom_point(aes(fill=sex),alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3")) + 
		    scale_fill_manual(values=c("tan1","steelblue1")) + 
		    labs(fill="Prediction", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"prediction_female.png"),device="png",width=100,height=60,units="mm")

            p1=ggplot(subset(plot.table,MappingFunction=="p0(k)" & sex=="Paternal"),aes(x=(markerB_cm-markerA_cm), y=prediction))+
            geom_point(aes(fill=sex),alpha=0.05,size=0.5,shape=21,stroke=0)+           
            geom_line(aes(x=x,y=y,col=mapping.function),data=mapping.functions,inherit.aes=F,size=0.5) +
		    ylim(0,0.60)+xlim(0,length.cm)+ xlab("Map distance (cM)")+ylab("Recombination frequency")+
		    scale_color_manual(values=c("#F0E442","#009E73","orchid","tan3","firebrick")) + 
		    scale_fill_manual(values=c("steelblue1")) + 
		    labs(fill="Prediction", colour="Mapping function",title=chr)  +  
		    theme(aspect.ratio=1,panel.grid.minor=element_blank())+  
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))
            ggsave(paste0(chr,"prediction_male.png"),device="png",width=100,height=60,units="mm")

            
            plot.table=rbind(subset(plot.table,sex=="Maternal"),subset(plot.table,sex=="Paternal"),make.row.names=F)
           
            ggplot(subset(plot.table, MappingFunction!="p0(k)*"),aes(x=(markerB_cm-markerA_cm), y=(recombination-prediction)))+
            geom_point(aes(col=sex,fill=sex),alpha=0.05,size=0.5,shape=21,stroke=0)+        
            geom_hline(yintercept=0,linetype=2)+   
            xlab("Map distance (cM)")+ylab("Error")+
		    scale_color_manual(values=c("tan1","steelblue1")) + 
		    scale_fill_manual(values=c("tan1","steelblue1")) + 
		    labs(fill="Sex",title=chr)  +    
		    guides(fill = guide_legend(override.aes = list(alpha=1,size=2))) +
		    theme(legend.position="none",text=element_text(size=5))+
		    facet_grid(sex~MappingFunction,scales="free_y")
		    ggsave(paste0(chr,"residual_scatter.png"),device="png",width=100,height=60,units="mm",dpi=250)
		    
		    
            sq=aggregate(plot.table$sq.error,by=list(plot.table$Chr,plot.table$sex,plot.table$MappingFunction),sum)
            sq$likelihood=aggregate(plot.table$likelihood,by=list(plot.table$Chr,plot.table$sex,plot.table$MappingFunction),prod)$x
			sq$sumOfLogLikelihood=aggregate(plot.table$loglikelihood,by=list(plot.table$Chr,plot.table$sex,plot.table$MappingFunction),sum)$x
			sq$TotalError=aggregate(plot.table$error,by=list(plot.table$Chr,plot.table$sex,plot.table$MappingFunction),sum)$x
            sq$TotalAbsoluteError=aggregate(plot.table$error.absolute,by=list(plot.table$Chr,plot.table$sex,plot.table$MappingFunction),sum)$x
			sq$TotalErrorLow95=aggregate(error ~ Chr + sex + MappingFunction, function(x) quantile(x,c(0.025,0.975))[1], data=plot.table)$error
			sq$TotalErrorUpper95=aggregate(error ~ Chr + sex + MappingFunction, function(x) quantile(x,c(0.025,0.975))[2], data=plot.table)$error
			sq$TotalErrorSd=aggregate(error ~ Chr + sex + MappingFunction, function(x) sd(x), data=plot.table)$error
	    
	   
			colnames(sq)[1:4]=c("chromosome","sex","MappingFunction","sq.sum")
			sum.of.squares=rbind(sum.of.squares,sq)
			n.markers=nrow(table1) 
			sum.of.squares$NumberOfMarkers=n.markers
			table.out=rbind(sum.of.squares,table.out)
  
			sq.total=aggregate(sum.of.squares$sq.sum,by=list(sum.of.squares$sex,sum.of.squares$MappingFunction),sum)
			colnames(sq.total)=c("Sex","MappingFunction","SumOfSquares")
			
			OUT[["sumOfSquares"]]=rbind(OUT[["sumOfSquares"]], sq.total)
			OUT[["rbar"]]=rbind(OUT[["rbar"]], rbar.tmp,rbar.tmp.function)
			
		}
# STEP 5 ends			

# STEP 6: Collect the results

      library(ggplot2)
	  plot.table=table.out #From Step 2
      
      plot.table$chromosome = factor(plot.table$chromosome,levels=paste0("chr",1:22))      
      best=aggregate(sumOfLogLikelihood ~ chromosome+sex,data=plot.table,max)
	  best=subset(plot.table,sumOfLogLikelihood %in% best$sumOfLogLikelihood)
	  best$LABEL="*"
	  best[best$MappingFunction=="Linear","LABEL"]="L"
	  best[best$MappingFunction=="Kosambi","LABEL"]="K"
	  best[best$MappingFunction=="Haldane","LABEL"]="H"
	  best[best$MappingFunction=="p0(k)","LABEL"]="p0(k)"
	  plot.table$MappingFunction=factor(plot.table$MappingFunction,labels=c("Haldane", "Kosambi", "Linear", "p0(k)"))
      
      ggplot(plot.table,aes(x=chromosome,y=sumOfLogLikelihood,fill=MappingFunction))+
      geom_jitter(shape=23,stroke=0.1,alpha=0.75,size=1,width=0.05,height = 0)+
      scale_color_manual(values=c("#F0E442","#009E73","orchid","tan2")) + 	
	  scale_fill_manual(values=c("#F0E442","#009E73","orchid","tan2")) +
	  ylab("log-likelihood of data")+xlab("Chromosome")+labs(fill="Function") + 
	  geom_text(data=best,aes(x=chromosome,y=0,label=LABEL),size=0.7,inherit.aes=F)+
	  expand_limits(y=0) +
	  facet_wrap(~sex,ncol=2,scales="free_y")+
	  guides(fill = guide_legend(override.aes = list(alpha=1,size=2)))+
	  theme(text=element_text(size=7),axis.text.x=element_text(angle=60, hjust = 1),panel.grid.minor=element_blank(),legend.position="bottom")
	  
	  plot.table.sum=aggregate(plot.table$sumOfLogLikelihood,by=list(plot.table$MappingFunction,plot.table$sex),sum)
      colnames(plot.table.sum)=c("MappingFunction","sex","sumOfLogLikelihood")
      ggplot(plot.table.sum,aes(x=MappingFunction,y=sumOfLogLikelihood,fill=MappingFunction))+geom_point(shape=23,alpha=0.8)+
	  scale_color_manual(values=c("#F0E442","#009E73","orchid","tan2")) + 	
	  scale_fill_manual(values=c("#F0E442","#009E73","orchid","tan2")) +
	  ylab("Likelihood of data")+xlab("Chromosome")+labs(fill="Mapping function")  + 
	  facet_wrap(~sex,ncol=1,scales="free_y")+
	  theme(legend.position="bottom",text=element_text(size=6))
 	  
	  total.marker.number=sum(unique(plot.table$NumberOfMarkers))
	  
	  plot.table.error=aggregate(cbind(TotalError, NumberOfMarkers) ~ sex + MappingFunction, sum, data=plot.table)
	  plot.table.error$MeanError=plot.table.error$TotalError/plot.table.error$NumberOfMarkers
      plot.table.error$TotalAbsoluteError=aggregate(TotalAbsoluteError ~ sex + MappingFunction, sum, data=plot.table)$TotalAbsoluteError
      plot.table.error$MeanAbsoluteError=plot.table.error$TotalAbsoluteError/plot.table.error$NumberOfMarkers
	  
	  plot.table$WeightedSd=plot.table$TotalErrorSd*(plot.table$NumberOfMarkers/total.marker.number)
	  
	  plot.table.error$TotalErrorSd=aggregate(TotalErrorSd ~ sex + MappingFunction, mean, data=plot.table)$TotalErrorSd
	  plot.table.error$TotalErrorSdWeighted=aggregate(WeightedSd ~ sex + MappingFunction, sum, data=plot.table)$WeightedSd
	  
	  ggplot(plot.table.error,aes(x=MappingFunction,y=MeanError,ymin=MeanError-TotalErrorSd,ymax=MeanError+TotalErrorSd,fill=MappingFunction))+
	  geom_errorbar(width=0.2) + geom_point(shape=23) +
	  scale_color_manual(values=c("#F0E442","#009E73","orchid","tan2")) + 	
	  scale_fill_manual(values=c("#F0E442","#009E73","orchid","tan2")) +
	  ylab("Error")+xlab("Mapping function")+labs(fill="Mapping function")  + 
	  facet_wrap(~sex,ncol=1,scales="free_y")+
	  theme(legend.position="bottom",text=element_text(size=6))
 	   
	  table.out=plot.table.sum
	  
	  OUT=list()
	  OUT[["likelihoods"]]=plot.table
	  sum.of.squares.total=aggregate(sq.sum ~ sex + MappingFunction, sum,data=plot.table)
	  OUT[["sumOfSquares"]]=sum.of.squares.total
	  OUT[["absoluteError"]]=plot.table.error

# STEP 6 ends

# STEP 7: calculate rbar

        table1=chromosomeLengths
        table.out=data.frame()
	    library(ggplot2)
		#############
		table1$LENGTH=as.numeric(table1$LENGTH)
		rbar=OUT[["rbar"]] #From Step 2
	    
	    rbar$IntraContrib=rbar$SumRij/rbar$Nrbarloci
	    rbar$IntraContrib2=rbar$IntraContrib*((rbar$ChrLength/sum(as.numeric(table1$LENGTH)))^2)
	    inter.chrom=data.frame(Species="Hsap",InterChromComponent=0.5*(1-sum((table1$LENGTH/sum(table1$LENGTH))^2)))
	    
	    OUT=list()
	    OUT[["rbar"]]=rbar
	    OUT[["interChrom"]]=inter.chrom
	    
	    table.out=aggregate(IntraContrib ~ Sex + Source, sum,data=rbar)
	    table.out$IntraContrib2=aggregate(IntraContrib2 ~ Sex + Source, sum,data=rbar)$IntraContrib2
	  
	    total.rbar=aggregate(IntraContrib2 ~ Source, sum,data=table.out)
	    total.rbar$varIBD=0.125*(1-total.rbar$IntraContrib2-2*inter.chrom$InterChromComponent)  
	    OUT[["ibdVarRbar"]]=total.rbar
	   
# STEP 7 ends
	   
