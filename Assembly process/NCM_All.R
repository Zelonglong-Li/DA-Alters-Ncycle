library(Hmisc)
library(minpack.lm)
library(stats4)
library(stringr)



spp<-read.table('AllOTU.tsv',head=T,stringsAsFactors=F,row.names=1,sep = "\t")

CK <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="CK"))
DA <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="DA"))
AOM <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="A"))
AD <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="AD"))

NCM <- function(grp){
  
        spp<-t(grp)
        
        N <- mean(apply(spp, 1, sum))
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
        spp.bi <- 1*(spp>0)
        freq <- apply(spp.bi, 2, mean)
        freq <- freq[freq != 0]
        C <- merge(p, freq, by=0)
        C <- C[order(C[,2]),]
        C <- as.data.frame(C)
        C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
        p <- C.0[,2]
        freq <- C.0[,3]
        names(p) <- C.0[,1]
        names(freq) <- C.0[,1]
        d = 1/N
        m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
        m.fit  
        m.ci <- confint(m.fit, 'm', level=0.95)
        freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
        pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
        Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
        Rsqr  
        Rsqr <- round(Rsqr,3)
        Nm <- round(coef(m.fit)*N,3)
        m <- round(coef(m.fit),3)
        
       
        
       
        bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
        bacnlsALL$Group <- deparse1(substitute(grp))
        bacnlsALL$Rsqr <- Rsqr
        bacnlsALL$Nm <- Nm
        bacnlsALL$m <- m
        bacnlsALL$U_L <-  "In prediction"
        bacnlsALL$U_L[bacnlsALL$freq <= bacnlsALL$Lower] <- "Below prediction"
        bacnlsALL$U_L[bacnlsALL$freq >= bacnlsALL$Upper] <- "Above prediction"
        return(bacnlsALL)
}

NCM_CK<- NCM(CK)  
NCM_DA<- NCM(DA)
NCM_AOM <- NCM(AOM)
NCM_AD<- NCM(AD)

Resultall <- rbind(NCM_CK,NCM_DA,NCM_AOM,NCM_AD<- NCM(AD))
Resultall$Group <- factor(Resultall$Group,levels = c('CK','DA','AD','AOM'))

Resultall$R_N <- paste("R^2=",Resultall$Rsqr,"     ","m=",Resultall$m)

p  <- ggplot(data = Resultall,aes(x=log10(p),y=freq,color=U_L)) +
    geom_point(position = "identity",size = .01)+
    scale_color_manual(values = c('#A52A2A',"#29A6A6","black"))+
    labs(x='Mean Relative Abundance (log10)',y='Frequency of Occurance')+
    # annotate(geom = "text",x=-Inf,y=Inf,vjust=1.5,hjust= -0.1,label= paste("R^2=",Resultall$Rsqr,"\n","Nm=",Resultall$Nm))+
    geom_line(aes(y=freq.pred),color="blue",lwd=1,alpha=0.5)+
    geom_line(aes(y=Lower),color="blue",lwd=0.5,lty=2)+
    geom_line(aes(y=Upper),color="blue",lwd=0.5,lty=2)+
    facet_wrap(~Group+R_N,nrow = 2,ncol = 2,scales = "free_x")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0,0),
          legend.position = c(0.25,0.05),
          legend.spacing.y = unit(0,"line"),
          legend.text = element_text(size = unit(7,"pt")),
          legend.margin = margin(rep(5,4),unit = "pt"),
          legend.key.size = unit(10,"pt"))
p  
ggsave("NCM.svg",p,width = 5,height = 5,units = "in")

