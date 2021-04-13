install.packages("ggfortify")
library(ggfortify)
install.packages("ggrepel")
library(ggrepel)


PrePCAProt<-read.csv("~/Downloads/pre_PCA.csv")

head(PrePCAProt[,1:10])
dim(PrePCAProt)
rownames(PrePCAProt)<-PrePCAProt$Go.Terms
head(PrePCAProt[,1:10])
PrePCAProtNew<-PrePCAProt[,2:dim(PrePCAProt)[2]]
head(PrePCAProtNew[,1:10])

PrePCAsamp<-PrePCAProtNew[,sample(dim(PrePCAProtNew)[2],100)]
# dim(PrePCAsamp)
# head(PrePCAsamp)
# colnames(PrePCAsamp)



tPrePCAsamp<-as.data.frame(t(PrePCAsamp))
head(tPrePCAsamp[,1:10])
colnames(tPrePCAsamp)
pcasamp<-prcomp(tPrePCAsamp)
pcasamp$rotation
dim(pcasamp$x)
dim(pcasamp$rotation)

pcasamp$rotation<-pcasamp$rotation[1:100,]
pcasamp$center<-pcasamp$center[1:100]
pcasamp$rotation
autoplot(pcasamp)
biplot(pcasamp)







tPrePCAProtNew<-t(PrePCAProtNew)
head(tPrePCAProtNew[,1:10])
tPCAProt<-prcomp(as.data.frame(t(PrePCAProtNew)))

# dim(tPCAProt$rotation)
# dim(tPCAProt$x)
# head(tPCAProt$rotation[,1:10])
# head(tPCAProt$x[,1:10])
# head(tPCAProt$center)


ROT<-as.data.frame(tPCAProt$rotation)

#Test the variables, see upper limit of PC1
ROTorder<-ROT[order(-ROT$PC1),]
head(ROTorder[,1:5],20)

#Good division is at PC>0.1
ROT$PC1<-abs(ROT$PC1)
ROTfilter<-ROT[ROT$PC1>0.1,]
dim(ROTfilter) #42 variables >0.1
ROTfilter<-ROTfilter[1:42,] #However many you want

SigVar<-rownames(ROTfilter) #extract the names of the biggest players
tPCAProtFilter<- tPCAProt
tPCAProtFilter$rotation<-tPCAProt$rotation[rownames(tPCAProt$rotation) %in% SigVar,] #only the top 42
dim(tPCAProtFilter$rotation) #Check that it worked!


write.csv(tPCAProtFilter$rotation[,1:10],file="~/Downloads/GOTermsPCA.csv")

          
?write.csv          
Trigrams<-as.data.frame(tPCAProtFilter$rotation)
Trigrams<-Trigrams[,1:10]
Trigrams$Zero<-rep(0,10)



tscandy<-autoplot(tPCAProt)
tscandy + ggtitle("GO Term Prediction Score PCA") +
  theme(plot.title=element_text(size=13, vjust=3)) +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_segment(data=Trigrams,aes(x=Zero,y=Zero,xend=PC1/5,yend=PC2/5,colour="magenta"),arrow=arrow(length=unit(0.5,"cm")),size=1.25) +
  geom_text_repel(data=Trigrams,aes(x=PC1/5,y=PC2/5,label=rownames(Trigrams)),colour="black",direction="y",hjust = -0.5,size=5,fontface="bold")+
  geom_text_repel(data=Trigrams,aes(x=PC1/5,y=PC2/5,label=rownames(Trigrams)),colour="cyan",direction="y",hjust = -0.5,size=5) + 
  theme(legend.position = "none")







z<-read.csv("~/Downloads/pre_PCAblast.csv")
head(z[,1:10])

rownames(z)<-z$Go.Terms
head(z[,1:10])
zNew<-z[,2:dim(z)[2]]
head(zNew[,1:10])

tz<-t(zNew)
head(tz[,1:10])
tPCAz<-prcomp(as.data.frame(tz))

head(tPCAz$rotation[,1:10])
head(tPCAz$x[,1:10])
tpanda<-autoplot(tPCAz)
tpanda + ggtitle("GO Term Prediction Score PCA") +
  theme(plot.title=element_text(size=13, vjust=3)) +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  theme(plot.title = element_text(hjust = 0.5))

