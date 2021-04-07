library(UsingR)
library(wmwpow) #calculate wilcoxon power
library(stats)


###Problem 1
x<-rnorm(40,1,2)
y<-rnorm(40,2,2)


a<-wilcox.test(x,y,alternative="two.sided")
b<-t.test(x,y,alternative="two.sided")
a$p.value
b$p.value

#Let's do this a bunch of times
f<-function(x=1,y=2){
  x1<-rnorm(40,x,2)
  y1<-rnorm(40,y,2)
  a1<-wilcox.test(x1,y1,alternative="two.sided")
  b1<-t.test(x1,y1,alternative="two.sided")
  return(c(a1$p.value,b1$p.value))
}

p_list<-c()
for(i in 1:100){
  p_list<-c(p_list,f())
}

#p_list      #It's a bunch of p-values


odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0
#create p-value lists for each tests; divide up big list into 2 lists
wp<-p_list[odd(1:length(p_list))]  #wilcox p-values
tp<-p_list[even(1:length(p_list))] #t-test p-values

#Proportion of p-values that are significant for each test
sum(wp<0.05)/100
sum(tp<0.05)/100


#Function power calculation for each of the tests
wmwpowd(40,40,alpha=0.05,distn="norm(1,2)",distm = "norm(2,2)")  #The Wilcoxon power is 0.575
power.t.test(40,1,2,0.05)   #0.598
#Formal power calculation for each of the tests

#Heck yah, the t-test has a power of 0.59 and wilcox 0.57 in our function and theirs!!!
#If they're not *perfectly* the same, it's because of randomness




###Problem 4.2
fms <- read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt", 
                  header=T, sep="\t")
NamesAkt1Snps <- names(fms)[substr(names(fms),1,5)=="actn3" | names(fms)=="HDL_C"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
#FMSgeno$HDL_C

#Looking at the frequency of each genotype
table(FMSgeno$actn3_r577x)
table(FMSgeno$actn3_rs540874)
table(FMSgeno$actn3_rs1815739)
table(FMSgeno$actn3_1671064)

#I defined variant allele as the more common homozygous genotype
FMSgeno$actn3_r577x<-ifelse(FMSgeno$actn3_r577x=="CC",0,1)
FMSgeno$actn3_rs540874<-ifelse(FMSgeno$actn3_rs540874=="GG",0,1)
FMSgeno$actn3_rs1815739<-ifelse(FMSgeno$actn3_rs1815739=="CC",0,1)
FMSgeno$actn3_1671064<-ifelse(FMSgeno$actn3_1671064=="AA",0,1)

(t1<-t.test(FMSgeno$HDL_C[FMSgeno$actn3_r577x==0],FMSgeno$HDL_C[FMSgeno$actn3_r577x==1]))
(t2<-t.test(FMSgeno$HDL_C[FMSgeno$actn3_rs540874==0],FMSgeno$HDL_C[FMSgeno$actn3_rs540874==1]))
(t3<-t.test(FMSgeno$HDL_C[FMSgeno$actn3_rs1815739==0],FMSgeno$HDL_C[FMSgeno$actn3_rs1815739==1]))
(t4<-t.test(FMSgeno$HDL_C[FMSgeno$actn3_1671064==0],FMSgeno$HDL_C[FMSgeno$actn3_1671064==1]))

p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value),method="bonferroni") #Bonferroni
p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value),method="BH") #Benjamini Hochberg

#Bonferroni results in *much* higher p-values than unadjusted. One of them is 1.00
#BH adjustments cause the p-values to slightly less significant than unadjusted (BH p-value:~0.2)
#Bonferroni corrects much more than BH -- perhaps overcorrecting
#BH seems to be a better choice...
#We're adjusting for 4 tests!





###Problem 4.3
fms <- read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt", 
                  header=T, sep="\t")
Res <- names(fms)[names(fms)=="resistin_c180g" | names(fms)=="CHOL" | names(fms)=="Gender"]
FMSgeno <- fms[,is.element(names(fms),Res)]
#FMSgeno$CHOL
summary(lm(CHOL~resistin_c180g,data=FMSgeno))      #The difference in cholesterol between CC & GG is significant.
summary(lm(CHOL~resistin_c180g+Gender,data=FMSgeno))    #stratify by gender this way
TukeyHSD(aov(data=FMSgeno,CHOL~resistin_c180g))    #The difference in cholesterol between any of the genotypes is not significant
TukeyHSD(aov(data=FMSgeno,CHOL~resistin_c180g+Gender))  #stratify by gender this way

###Stratify by gender with a subset
###Only the bros in this section
FMSgenoM<-FMSgeno[FMSgeno$Gender=="Male",]
#Same analyses but only for Males
summary(lm(CHOL~resistin_c180g,data=FMSgenoM)) #The difference in cholesterol between CC & GG is significant.
TukeyHSD(aov(data=FMSgenoM,CHOL~resistin_c180g)) #The difference in cholesterol between any of the genotypes is not significant


###Only the ladies in this section
FMSgenoF<-FMSgeno[FMSgeno$Gender=="Female",]
#Same analyses but only for Males
summary(lm(CHOL~resistin_c180g,data=FMSgenoF)) #The difference in cholesterol between CC & GG is significant.
TukeyHSD(aov(data=FMSgenoF,CHOL~resistin_c180g)) #The difference in cholesterol between any of the genotypes is not significant


#For both genders individually there is no significant association between total Cholesterol and the resistin_c180g SNP
#When both genders are together, GG and CC were significantly different, but not after multiple test correction w/Tukey
#In other words, the GG genotype is associated with a significantly lower total cholesterol than the CC genotype *but not after correction*







###Problem 4.6
#Free Step-down resampling -- FAMUSS
fms <- read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt", 
                  header=T, sep="\t")
NamesAkt1Snps <- names(fms)[substr(names(fms),1,5)=="actn3" | names(fms)=="NDRM.CH"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]

table(FMSgeno$actn3_r577x) #Major is C
table(FMSgeno$actn3_rs540874)  #Major is G
table(FMSgeno$actn3_rs1815739)  #Major is C
table(FMSgeno$actn3_1671064)  #Major is A

Actn3Bin <- data.frame(FMSgeno$actn3_r577x!="CC", FMSgeno$actn3_rs540874!="GG",
                       FMSgeno$actn3_rs1815739!="CC", 
                       FMSgeno$actn3_1671064!="AA") #Minor alleles are T,A,T,and G respectively

quantile(FMSgeno$NDRM.CH,c(0.25,0.75),na.rm=TRUE)  #number > 66.7 are top quartile
Actn3Bin$NDRM.CH<-FMSgeno$NDRM.CH>66.7

colnames(Actn3Bin)

#How do linear model w/ only binary variables?
linmod<-lm(data=Actn3Bin,NDRM.CH~.)
(Mod<-summary(linmod))


TestStatObs <- Mod$coefficients[-1,3]      #pick out the good coefficients
Tobs <- as.vector(sort(abs(TestStatObs)))  #sort by strength
MissDat <- apply(is.na(Actn3Bin),1,any)    #grab locations of NAs
Actn3BinC <- Actn3Bin[!MissDat,]           #Fix the NAs
sum(is.na(Actn3BinC$NDRM.CH))              #Check that it worked
Ord <- order(abs(TestStatObs))             #Order by strength
M <- 1000                                  #Number of reps
NSnps <- 4 #Number of SNPs
Nobs <- sum(!MissDat)
TestStatResamp <- matrix(nrow=M,ncol=NSnps)
for (i in 1:M){ 
  Ynew <- sample(Mod$residuals,size=Nobs,replace=T)
  ModResamp <- summary(lm(Ynew~.,data=Actn3BinC))
  TestStatResamp[i,] <- 
    abs(ModResamp$coefficients[-1,3])[Ord]
}
Qmat <- t(apply(TestStatResamp, 1, cummax))
Padj <- apply(t(matrix(rep(Tobs,M), NSnps)) < Qmat, 2, mean)

#Qmat
Padj
#Sweet! We know from this technique that none of our p-values are significant before or after resampling






###Problem 4.7
hgdp_akt1 <- read.delim("http://www.biostat.umn.edu/~cavanr/HGDP_AKT1.txt", header=T, sep="\t",stringsAsFactors = TRUE)
levels(hgdp_akt1$Gender)
table(hgdp_akt1$AKT1.C0756A) #AA is minor
table(hgdp_akt1$AKT1.C6024T) #TT is minor
table(hgdp_akt1$AKT1.G2347T) #TT is minor
table(hgdp_akt1$AKT1.G2375A) #AA is minor


#Do them all together!!!!!!~!!!!~!~!~~!@#$%^&_
#The trait we choose is completely irrelevant to the effective number of tests
#mEff only takes into account genetics data related to linkage disequilibrium
AKT1_Bin<-data.frame(hgdp_akt1$AKT1.C0756A!="AA",hgdp_akt1$AKT1.C6024T!="TT",
                      hgdp_akt1$AKT1.G2347T!="TT",hgdp_akt1$AKT1.G2375A!="AA")
MissDat <- apply(is.na(AKT1_Bin),1,any) #Find all rows where any column has value NA
AKT1_BinThere <- AKT1_Bin[!MissDat,] #remove missing values, matrix of TRUE/FALSE


corAKT1 <- cor(AKT1_BinThere)
eig <- eigen(corAKT1)$values
mEff <- 1+(4-1)*(1-var(eig)/4)
mEff   #Total effective tests!

#The effective number of tests is 2.077








###Problem 6
#2000 columns and 100 rows
#rbinom(p=0.5)
#mu=SNP[,1]*10+SNP[,2]*9 + SNP[,3]*8 + ...
#y = rnorm(100,mu,sd=1)

longlist<-rbinom(n=2000*100,size=1,p=0.5) #a bunch of 0s and 1s
length(longlist)

SNPmat<-matrix(longlist,nrow=100) #Divide into our dataframe
dim(SNPmat) #Yup it worked!
head(SNPmat[,1:10]) #Yup it worked!!!

SNPdat<-as.data.frame(SNPmat)
head(SNPdat[,1:10])


mu=c(rep(0,100))
for(i in 1:10){   #For loop creating mu based on first 10 markers
  mu=mu+(10-i+1)*SNPdat[,i]
}
mu
SNPdat$Outcome<-rnorm(100,mu,1)  #add some randomness so it's not perfect

summary(lm(data=SNPdat,Outcome~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10)) #Check that Outcome IS correlated
#Outcome is indeed correlated
summary(lm(data=SNPdat,Outcome~V4)) #not significant for 1 variable

###NOTE###
#Linear regression is significant, but t-test is not bc conditial probability
###NOTE###

#Doing one column
t.test(SNPdat$Outcome[SNPdat$V4==1],SNPdat$Outcome[SNPdat$V4==0]) 

#DOING ALL THE VARIABLES
listy<-c()
for(i in 1:(dim(SNPdat)[2]-1)){
  listy<-c(listy,t.test(SNPdat$Outcome[SNPdat[,i]==1],SNPdat$Outcome[SNPdat[,i]==0])$p.value)
}
head(listy,25) #Yup! Seems to work!

hist(listy)



if (!requireNamespace("BiocManager", quietly =TRUE)) +
  install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)
qs<-qvalue(listy)
qs$qvalues[1:25] #it seems to have worked
hist(qs$qvalues)

#Here's the local FDR
qs$lfdr    

#Here's the manual calculation for FDR
#We know that Markers 1:10 are the only relevant variables, the others are false positives
sum(qs$qvalues[11:length(qs$qvalues)]<0.05)
#Our FDR is 0 because there are 0 false positives and 3 predicted positives







###Problem 7
set.seed(1)
normies<-rnorm(100,0,1)
weirdos<-rnorm(50,3,0.5)
head(normies); head(weirdos)

people<-c(normies,weirdos)
hist(people) #wow, incredible; mixed sample distribution


B <-1000
sim <- rep(NA,1000)
for (i in 1:B){
  Samps <- sample(people, replace=T)
  sim[i]=sd(Samps)
  }
sim
hist(sim) #bootstrap distribution
quantile(sim,c(0.025,0.975))

#The confidence interval is 1.419706 & 1.649527 (on MY computer)








