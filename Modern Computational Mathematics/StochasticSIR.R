######Start Tests
numAgents=10000
currState=c(numAgents,10,0)
currState


numSteps <- 100

S=c()
I=c()
R=c()

for(k in 1:10){
  numAgents=10000
  currState=c(numAgents,10,0)
  S=c()
  I=c()
  R=c()
for(m in 1:50){
  S=c(S,currState[1])
  I=c(I,currState[2])
  R=c(R,currState[3])
  M=matrix(c(1-0.5*currState[2]/numAgents,0.5*currState[2]/numAgents,0,
             0,0.95,0.05,
             0,0,1),nrow=3);M
  nextState1 <- sample(1:3, currState[1], rep=TRUE, prob=M[,1])
  nextState2 <- sample(1:3, currState[2], rep=TRUE, prob=M[,2])
  nextState3 <- sample(1:3, currState[3], rep=TRUE, prob=M[,3])
  
  fromState1 <- as.vector(tabulate(nextState1,nbins=3))
  fromState2 <- as.vector(tabulate(nextState2,nbins=3))
  fromState3 <- as.vector(tabulate(nextState3,nbins=3))
  
  
  currState <- fromState1 + fromState2 + fromState3
  
}
currState



if(k==1){
plot(S~seq(1,length(S)),type="l",col="red",main="SIR Plotted Stochastically",ylab="Number of People",xlab="Days Passed")
lines(I~seq(1,length(I)),col="blue")
lines(R~seq(1,length(R)),col="magenta")
}else{
  lines(S~seq(1,length(S)),col="red")
  lines(I~seq(1,length(I)),col="blue")
  lines(R~seq(1,length(R)),col="magenta")
}
}
legend("topright",legend=c("Susceptible","Infected","Recovered"),col=c("red","blue","magenta"),type="lty")



####Deterministic Method
library(deSolve)
sir=function(t,y,parms){ 
  beta=parms[1] 
  gamma=1/13 
  dS=-beta*y[1]*y[2] 
  dI=beta*y[1]*y[2]-gamma*y[2]
  dR=gamma*y[2] # equation for dR/dt
  dy=c(dS,dI,dR)
  list(c(dS,dI,dR))
}

yini=c(160000,10,0) # define S(0), I(0) and R(0)
time=seq(0,60,length=100) # define time steps
parms=c(beta=5e-6,gamma=1/3) # define parameters
soln=ode(y=yini,func=sir,times=time,parms=parms) # solve ode, assign output to "soln"

plot(soln[,1],soln[,2],col="red",type="l")
lines(soln[,1],soln[,3],col="blue")
lines(soln[,1],soln[,4],col="magenta")
legend("topright",legend=c("Susceptible","Infected","Recovered"),col=c("red","blue","magenta"),lty=1)


niamey=read.csv("~/Downloads/niamey.csv")





  
  
  
  
######Start Tests MSEIRS
numAgents=100
currState=c(50,numAgents,50,50,0)
currState


numSteps <- 100

M<-c()
S=c()
E=c()
I=c()
R=c()


numAgents=8000
numAgents
for(k in 1:10){
  currState=c(numAgents-7600,7500,50,50,0,0)
  M<-c()
  S=c()
  E=c()
  I=c()
  R=c()
  totAgents=numAgents

  
  for(m in 1:300){
    M=c(M,currState[1])
    S=c(S,currState[2])
    E=c(E,currState[3])
    I=c(I,currState[4])
    R=c(R,currState[5])
    
    B=9.6468
    del=0.1059
    mew=0.0001024138
    beta=0.1014887
    eps=0.148299
    gamma=0.051054
    
    Mat=matrix(c(1-del-mew,del,0,0,0,mew,0,1-beta*currState[4]/totAgents-mew,beta*currState[4]/totAgents,
               0,0,mew,0,0,1-eps-mew,eps,0,mew,0,0,0,1-gamma-mew,gamma,mew,0,0,0,0,1-mew,mew,0,0,0,0,0,1),nrow=6);Mat
    nextState1 <- sample(1:6, currState[1], rep=TRUE, prob=Mat[,1])
    nextState2 <- sample(1:6, currState[2], rep=TRUE, prob=Mat[,2])
    nextState3 <- sample(1:6, currState[3], rep=TRUE, prob=Mat[,3])
    nextState4 <- sample(1:6, currState[4], rep=TRUE, prob=Mat[,4])
    nextState5 <- sample(1:6, currState[5], rep=TRUE, prob=Mat[,5])
    nextState6 <- sample(1:6, currState[6], rep=TRUE, prob=Mat[,6])
    
    
    fromState1 <- as.vector(tabulate(nextState1,nbins=6))
    fromState2 <- as.vector(tabulate(nextState2,nbins=6))
    fromState3 <- as.vector(tabulate(nextState3,nbins=6))
    fromState4 <- as.vector(tabulate(nextState4,nbins=6))
    fromState5 <- as.vector(tabulate(nextState5,nbins=6))
    fromState6 <- as.vector(tabulate(nextState6,nbins=6))

    
    currState <- fromState1 + fromState2 + fromState3 + fromState4 + fromState5 + fromState6
    currState[1]<-currState[1]+B
    totAgents<-numAgents -currState[6] +B
  }
  currState

  if(k==1){
    plot(S~seq(1,300,by=300/length(S)),ylim=c(0,8000),type="l",col="orange",xlab="Days Passed",ylab="Number of People",main="Number of People per Compartment in MSEIR Model")
    lines(M~seq(1,300,by=300/length(M)),col="red")
    lines(E~seq(1,300,by=300/length(E)),col="yellow")
    lines(I~seq(1,300,by=300/length(I)),col="green")
    lines(R~seq(1,300,by=300/length(R)),col="blue")
  }else{
    lines(S~seq(1,300,by=300/length(S)),type="l",col="orange")
    lines(M~seq(1,300,by=300/length(M)),col="red")
    lines(E~seq(1,300,by=300/length(E)),col="yellow")
    lines(I~seq(1,300,by=300/length(I)),col="green")
    lines(R~seq(1,300,by=300/length(R)),col="blue")
  }
}
legend("topright",legend=c("Maternal Immunity","Susceptible","Exposed", "Infected","Recovered"),col=c("red","orange",
                              "yellow","green","blue"),lty=1)
points(niamey$cases~niamey$days)







#MSIR model

numAgents=8000
numAgents
for(k in 1:10){
  currState=c(numAgents-7600,7500,50,50,0)
  M<-c()
  S=c()
  I=c()
  R=c()
  totAgents=numAgents

  
  for(m in 1:300){
    M=c(M,currState[1])
    S=c(S,currState[2])
    I=c(I,currState[3])
    R=c(R,currState[4])
    
    # B=rnorm(1,30,sd=1)
    # del=rnorm(1,0.1,sd=0.0025)
    # mew=rnorm(1,0.0001,sd=3*0.00001/10)
    # beta=rnorm(1,0.1,sd=0.025)
    # eps=rnorm(1,0.15,sd=0.0025)
    # gamma=rnorm(1,0.05,sd=0.0025)
    
    B=20
    del=0.1
    mew=0.0001
    beta=0.1
    gamma=0.05
    
    Mat=matrix(c(1-del-mew,del,0,0,mew,0,1-beta*currState[3]/totAgents-mew,beta*currState[3]/totAgents,0,mew,
                 0,0,1-gamma-mew,gamma,mew,0,0,0,1-mew,mew,0,0,0,0,1),nrow = 5);Mat
    nextState1 <- sample(1:5, currState[1], rep=TRUE, prob=Mat[,1])
    nextState2 <- sample(1:5, currState[2], rep=TRUE, prob=Mat[,2])
    nextState3 <- sample(1:5, currState[3], rep=TRUE, prob=Mat[,3])
    nextState4 <- sample(1:5, currState[4], rep=TRUE, prob=Mat[,4])
    nextState5 <- sample(1:5, currState[5], rep=TRUE, prob=Mat[,5])

    
    fromState1 <- as.vector(tabulate(nextState1,nbins=5))
    fromState2 <- as.vector(tabulate(nextState2,nbins=5))
    fromState3 <- as.vector(tabulate(nextState3,nbins=5))
    fromState4 <- as.vector(tabulate(nextState4,nbins=5))
    fromState5 <- as.vector(tabulate(nextState5,nbins=5))

    
    
    currState <- fromState1 + fromState2 + fromState3 + fromState4 + fromState5
    currState[1]<-currState[1]+B
    totAgents<-numAgents -currState[5] +B
  }
  currState
  
  if(k==1){
    plot(S~seq(1,300,by=300/length(S)),ylim=c(0,8000),type="l",col="orange",xlab="Days Passed",ylab="Numper of People",main="MSIR Plotted Stochastically")
    lines(M~seq(1,300,by=300/length(M)),col="red")
    lines(I~seq(1,300,by=300/length(I)),col="green")
    lines(R~seq(1,300,by=300/length(R)),col="blue")
  }else{
    lines(S~seq(1,300,by=300/length(S)),type="l",col="orange")
    lines(M~seq(1,300,by=300/length(M)),col="red")
    lines(I~seq(1,300,by=300/length(I)),col="green")
    lines(R~seq(1,300,by=300/length(R)),col="blue")
  }
}
legend("topright",legend=c("Maternal Immunity","Susceptible", "Infected","Recovered"),col=c("red","orange","green","blue"),lty=1)
points(niamey$cases~niamey$days)


