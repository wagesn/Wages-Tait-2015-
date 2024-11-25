###install required R packages
library(binom)
library(nnet)
library(dfcrm)


###Load the function 'bpocrm' 
bpocrm<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb,n.ar){
 
# if a single ordering is inputed as a vector, convert it to a matrix
	if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel));
	
	nord.tox = nrow(p.skel);
	mprior.tox = rep(1/nord.tox, nord.tox);  # prior for each toxicity ordering

# if a single ordering is inputed as a vector, convert it to a matrix
	if(is.vector(q.skel)) q.skel=t(as.matrix(q.skel));
	
	nord.eff = nrow(q.skel);
	mprior.eff = rep(1/nord.eff, nord.eff);  # prior for each efficacy ordering

post.tox<-function(a,p,y,n){
	s2=1.34
	lik=1
	for(j in 1:length(p)){
		pj=p[j]**exp(a)
		lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
		}
	return(lik*exp(-0.5*a*a/s2));
    }

#the posterior mean of ptox
    posttoxf <- function(a, p, y, n, j) { p[j]^(exp(a))*post.tox(a, p, y, n); }
	
post.eff<-function(b,q,z,n){
	s2=1.34
	lik=1
	for(j in 1:length(q)){
		qj=q[j]**exp(b)
		lik=lik*qj^z[j]*(1-qj)^(n[j]-z[j]);
		}
	return(lik*exp(-0.5*b*b/s2));
    }

#the posterior mean of peff
    postefff <- function(b, q, z, n, j) {q[j]^(exp(b))*post.eff(b, q, z, n); }



### run a trial 	
    ncomb = ncol(p.skel);   #number of combos
    y=rep(0,ncomb);  #number of toxicity/responses at each dose level
    z=rep(0,ncomb);   #number of efficacy at each dose level
    n=rep(0,ncomb);  #number of treated patients at each dose level
    comb.curr = start.comb;  # current dose level	 
    ptox.hat = numeric(ncomb); # estimate of toxicity prob
    peff.hat = numeric(ncomb); # estimate of efficacy prob
    comb.select=rep(0,ncomb); # a vector of indicators for dose selection
    stop=0; #indicate if trial stops early

 for(i in 1:ncohort)
    {
	# generate data for a new cohort of patients
	
		y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,p0[comb.curr]);
		z[comb.curr] = z[comb.curr] + rbinom(1,cohortsize,q0[comb.curr]);
		n[comb.curr] = n[comb.curr] + cohortsize;

		
		marginal.tox = rep(0, nord.tox);
		for(k in 1:nord.tox)
		{
			marginal.tox[k] = integrate(post.tox,lower=-Inf,upper=Inf, p=p.skel[k,], y=y,n=n)$value;
		}
		
		postprob.tox = (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox);

		marginal.eff = rep(0, nord.eff);
		for(k in 1:nord.eff)
		{
			marginal.eff[k] = integrate(post.eff,lower=-Inf,upper=Inf, q=q.skel[k,], z=z, n=n)$value;
		}
		
		postprob.eff = (marginal.eff*mprior.eff)/sum(marginal.eff*mprior.eff);
		
		# toxicity model selection, identify the model with the highest posterior prob
			if(nord.tox>1){ 
				mtox.sel = which.is.max(postprob.tox); 
			} else{
				mtox.sel = 1;
			}

		# efficacy model selection, identify the model with the highest posterior prob
			if(nord.eff>1){
				meff.sel = which.is.max(postprob.eff); 
			} else{
				meff.sel = 1;
			}

		# calculate posterior mean of toxicity probability at each combo
			for(j in 1:ncomb){
				ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, p.skel[mtox.sel,], y, n,j)$value/marginal.tox[mtox.sel]; 
			}
		

		# calculate posterior mean of efficacy probability at each combo
			for(j in 1:ncomb){
				peff.hat[j] = integrate(postefff,lower=-Inf,upper=Inf, q.skel[meff.sel,], z, n,j)$value/marginal.eff[meff.sel]; 
			}
		
		
		aset=which(ptox.hat<=tul)
		if(length(aset)==0){aset=which.min(ptox.hat)}
		peff.hat.aset=rep(0,ncomb)
		peff.hat.aset[aset]=peff.hat[aset]
		ar.prob=peff.hat.aset/sum(peff.hat.aset)

		if(length(aset)==1){
			comb.best=aset
		} else {
			ifelse(sum(n)<=n.ar,comb.best<-sample(1:ncomb,1,prob=ar.prob),comb.best<-which.max(peff.hat.aset))
 		}

	try=length(n[n>0])
	if(try<ncomb){
	if(comb.best==comb.curr || comb.best<comb.curr){
		comb.curr<-comb.best
	} else{
		comb.curr<-comb.curr+1
		}
	} else{
		comb.curr<-comb.best
	}


		##########stopping rules
		safety=binom.confint(y[1],n[1],conf.level=0.95,methods="exact")$lower
		if(safety>tul){
			stop=1
			break
			}

		if(sum(n) > n.ar){
		futility=binom.confint(z[comb.curr],n[comb.curr],conf.level=0.95,methods="exact")$upper
		if(futility<ell){
			stop=2
			break
			}
   		   } 
		}
	if(stop==0){
		comb.select[comb.curr]=comb.select[comb.curr]+1;
		}
	return(list(comb.select=comb.select,tox.data=y,eff.data=z,pt.allocation=n,stop=stop))
}
##########'bpocrm' end here

###Load the function 'bpocrm.sim' 
bpocrm.sim<-function(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial){
	ncomb=length(p0)
	
	comb.select<-y<-z<-n<-matrix(nrow=ntrial,ncol=ncomb)
	nstop=0
	
	for(i in 1:ntrial){
		result<-bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb,n.ar)
		comb.select[i,]=result$comb.select
		y[i,]=result$tox.data
		z[i,]=result$eff.data
		n[i,]=result$pt.allocation
		nstop=nstop+result$stop
	}
	cat("True tox probability: ", round(p0,3), sep="    ",  "\n");
	cat("True eff probability: ", round(q0,3), sep="    ",  "\n");
	cat("selection percentage: ", formatC(colMeans(comb.select)*100, digits=1, format="f"), sep="    ",  "\n");
	cat("number of toxicities:    ", formatC(colMeans(y), digits=1, format="f"), sep="    ",   "\n");
	cat("number of responses:    ", formatC(colMeans(z), digits=1, format="f"), sep="    ",   "\n");
	cat("number of patients treated:     ", formatC(colMeans(n), digits=1, format="f"), sep="    ",   "\n");
	cat("percentage of stop:    ", nstop/ntrial*100, "\n");
}
##########'bpocrm.sim' end here




####################################
#
#	Input
#
#	p0 = true toxicity probabilities 
#	q0 = true efficacy probabilities
#	p.skel = toxicity skeleton values
# 	q.skel = efficacy skeleton values
#	tul = toxicity upper limit allowed
# 	ell = efficacy lower limit allowed
#	cohortsize = size of each cohort inclusion
# 	nchort = number of cohorts in a sinlge trial
#	start.comb = starting combination
#	n.ar = size of adaptive randomization phase
#
#	ntrial = number of simulated trials
# 
#	
#
####################################


######################################
#
#        Generate results in Table 2
#
######################################


#####Specify the total number of combinations
d<-6

#####Specify the number of possible toxicity orderings
s<-1   

###Specify a set of toxicity skeleton values
p.skel<-c(0.01,0.08,0.15,0.22,0.29,0.36)

#####Specify the number of possible toxicity orderings
g<-11   #efficacy

###Specifiy the possible efficacy orderings of the drug combinations
q.skel<-matrix(nrow=g,ncol=d)
q.skel[1,]<-c(0.60,0.50,0.40,0.30,0.20,0.10)  
q.skel[2,]<-c(0.50,0.60,0.50,0.40,0.30,0.20)  
q.skel[3,]<-c(0.40,0.50,0.60,0.50,0.40,0.30)  
q.skel[4,]<-c(0.30,0.40,0.50,0.60,0.50,0.40)  
q.skel[5,]<-c(0.20,0.30,0.40,0.50,0.60,0.50)  
q.skel[6,]<-c(0.10,0.20,0.30,0.40,0.50,0.60)  
q.skel[7,]<-c(0.20,0.30,0.40,0.50,0.60,0.60)  
q.skel[8,]<-c(0.30,0.40,0.50,0.60,0.60,0.60)  
q.skel[9,]<-c(0.40,0.50,0.60,0.60,0.60,0.60)  
q.skel[10,]<-c(0.50,0.60,0.60,0.60,0.60,0.60)  
q.skel[11,]<-c(rep(0.60,6))  
 

p0<-c(0.05,0.10,0.20,0.28,0.50,0.50)   ##true toxicity probability scenario
q0<-c(0.05,0.23,0.47,0.70,0.70,0.70)
#c(0.05,0.13,0.25,0.38,0.50,0.63)   ##true efficacy probability scenario
tul<-0.33 ##toxicity upper limit 
ell<-0.05 ##efficacy lower limit
cohortsize=1 ##cohort size for each inclusion
ncohort=64   ##number of cohorts
start.comb=1 ##starting combination
n.ar=16      ##size of AR phase
ntrial=250   ##number of simulated trials 
set.seed(580)  ##random seed

##simulate a single trial
#bpocrm(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,start.comb,n.ar)

##simulate many trials
bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial)


######################################
#
#        Generate results in Table 3
#
######################################

#####Specify the total number of doses
d<-5

#####Specify the number of possible toxicity orderings
s<-1   

###Specify a set of toxicity skeleton values
p.skel<-c(0.01,0.08,0.15,0.22,0.29)

#####Specify the number of possible toxicity orderings
g<-9   #efficacy

###Specifiy the possible efficacy orderings of the doses
q.skel<-matrix(nrow=g,ncol=d)
 q.skel[1,]<-c(0.60,0.70,0.60,0.50,0.40)
 q.skel[2,]<-c(0.70,0.60,0.50,0.40,0.30)
 q.skel[3,]<-c(0.50,0.60,0.70,0.60,0.50)
 q.skel[4,]<-c(0.40,0.50,0.60,0.70,0.60)
 q.skel[5,]<-c(0.30,0.40,0.50,0.60,0.70)
 q.skel[6,]<-c(0.70,0.70,0.70,0.70,0.70) 
 q.skel[7,]<-c(0.60,0.70,0.70,0.70,0.70)
 q.skel[8,]<-c(0.50,0.60,0.70,0.70,0.70) 
 q.skel[9,]<-c(0.40,0.50,0.60,0.70,0.70)

p1<-c(0.01,0.05,0.10,0.15,0.20)
q1<-c(0.30,0.50,0.60,0.40,0.25)

p2<-c(0.02,0.06,0.12,0.30,0.40)
q2<-c(0.38,0.50,0.40,0.30,0.25)

p3<-c(0.03,0.09,0.16,0.28,0.42)
q3<-c(0.25,0.35,0.48,0.65,0.52)

p4<-c(0.02,0.05,0.07,0.09,0.11)
q4<-c(0.68,0.56,0.49,0.40,0.33)

tul<-0.33 ##toxicity upper limit 
ell<-0.20 ##efficacy lower limit
cohortsize=1 ##cohort size for each inclusion
ncohort=48   ##number of cohorts
start.comb=1 ##starting dose
n.ar=36      ##size of AR phase
ntrial=1000   ##number of simulated trials 
set.seed(580)  ##random seed


p0<-p4
q0<-q4
##simulate many trials
bpocrm.sim(p0,q0,p.skel,q.skel,tul,ell,cohortsize,ncohort,ntrial)


