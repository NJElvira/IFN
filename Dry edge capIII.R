###08/07/2024 Nur J elvira
###Dry edge capIII
###data from script ifn2merge3

library(gridExtra)
library(lme4); library(glmmTMB); library(dotwhisker)
library(broom); library(broom.mixed)
library(r2glmm)
library(MuMIn)
library(merDeriv)
library(ks)
library(sjPlot)
library(sjmisc)

library(plyr); library(dplyr);library(jtools)

library(proxy)

###PINUS HALEPENSIS####
ph_inter <-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/EUF_Phalepensis_PCA12Bio.csv")
hpi_inter <- Hscv (x=as.matrix(ph_inter[, c(13,14)]),  binned=T, pilot="samse")# estimating kernel bandwidth with crossvalidation, it can be also selected manually
niche_volume_inter <- kde(x=ph_inter[, c(13,14)], gridsize =  2500, H=hpi_inter)# weighted by specie' abundance #
niche_rast_inter <- raster(niche_volume_inter)# do that in the complex way (niche_volume$estimate, xmx, etc) give a especular plot
niche_rast_inter@data@values[which(niche_rast_inter@data@values<niche_volume_inter$cont[95])] <- 0# here select the desired percentile
niche_inter_spdf <- as(niche_rast_inter, "SpatialPixelsDataFrame")
niche_inter_df <- data.frame(niche_inter_spdf@coords, niche_inter_spdf@data)
niche_inter_df[niche_inter_df$layer< contourLevels(niche_volume_inter, prob=0.05), "layer"] <- 0

cograv <- COGravity(x=niche_inter_df[, 1], y=niche_inter_df[, 2], 
                    wt=niche_inter_df[, 3]^40)# it could also work with lower exponent
cograv# option a centre of mass
centroid <- data.frame(matrix(nrow=1, ncol=2))
centroid[, 1] <- cograv[[1]]
centroid[, 2] <- cograv[[3]]

## coordenades centroide
centroid
xcph <- 0.3054232
ycph <- 0.9686951 

niche_rast1 <- niche_rast_inter
niche_rast1@data@values[which(is.na(niche_rast1@data@values))] <- 0
niche_rast1@data@values[niche_rast1@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
a <- rasterToContour(niche_rast1)
b <-  as(a, "SpatialPointsDataFrame")
b <- data.frame(b@coords)
names(b) <- c("pca1", "pca2")

limitph<-b



#datbase con los ejes
DPhalpreprova <- ddply(DistPhalpre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                      coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                      edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                      xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                      dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                      drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                      drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                      ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                      sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire),
                      pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                      pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                      pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                      pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPhalpostprova <- ddply(DistPhalpost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire),
                       pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                       pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

## selecció//definició del wet edge (queno s'ha fet servir) i el dry edge
wet_points_limit <- rbind(subset(limitph,pca1>-5)%>%subset(.,pca2<(-1)),subset(limitph,pca1<(-2.5))%>%subset(.,pca2>-4))
## !!! dry limit només agafant la projecció de bio04 i la de bio6 ~~ aproximadament a 0.51 de l'eix 1 del PCA
dry_points_limit5 <- rbind(subset(limitph,pca1>0.5)%>%subset(.,pca2<3),subset(limitph,pca1>2.5)%>%subset(.,pca2>-4))

wpl <- wet_points_limit[,1:2]
dpl5 <- dry_points_limit5[,1:2]
colnames(wpl) <- c('pc1l','pc2l')
colnames(dpl5) <- c('pc1l','pc2l')

##distance to the wet and dry limits of all observations pre fire
DistPhalpre$wetdist <- apply(DistPhalpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistPhalpre$drydist5 <- apply(DistPhalpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistPhalpre$fedge <- apply(DistPhalpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitph[,1:2],method='Euclidean')))
DistPhalpre$fcore <- apply(DistPhalpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcph,ycph),ncol=2)),method='Euclidean')))

##distance to the wet and dry limits of all observations post fire
DistPhalpost$wetdist <- apply(DistPhalpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistPhalpost$drydist5 <- apply(DistPhalpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistPhalpost$fedge <- apply(DistPhalpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitph[,1:2],method='Euclidean')))
DistPhalpost$fcore <- apply(DistPhalpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcph,ycph),ncol=2)),method='Euclidean')))

####sACAR LAS OCURRENCIAS FUERA DE LIMITE Y PONDERAR LAS DISTANCIAS BORDE SECO
###PREFIRE
inoutPh <- in.kde(x=DistPhalpre[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutPh)
DistPhalpre <- cbind(DistPhalpre, inoutPh)

DistPhalpre$DistWDE
indice_false <- which(inoutPh == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistPhalpre[k,'DistWDE'] <- DistPhalpre[k,"min.dist.cen"] / (DistPhalpre[k,"min.dist.cen"] - DistPhalpre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutPh == TRUE)
for(k in indice_true){
  
  DistPhalpre[k,'DistWDE'] <- DistPhalpre[k,"min.dist.cen"] / (DistPhalpre[k,"min.dist.cen"] + DistPhalpre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}
###POSTFIRE
inoutPh <- in.kde(x=DistPhalpost[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutPh)
DistPhalpost <- cbind(DistPhalpost, inoutPh)

DistPhalpost$DistWDE
indice_false <- which(inoutPh == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistPhalpost[k,'DistWDE'] <- DistPhalpost[k,"min.dist.cen"] / (DistPhalpost[k,"min.dist.cen"] - DistPhalpost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutPh == TRUE)
for(k in indice_true){
  
  DistPhalpost[k,'DistWDE'] <- DistPhalpost[k,"min.dist.cen"] / (DistPhalpost[k,"min.dist.cen"] + DistPhalpost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}


##base datos pre y severidad
DPhalpreprova <- ddply(DistPhalpre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
                       #pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       #pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       #pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                       #pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPhalpreprova$mindpl5 <- apply(DPhalpreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DPhalpreprova$minwpl <- apply(DPhalpreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

sev_ph<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/DataFirePhal_Severity.csv")

DPhalpreprova <-dplyr::full_join(DPhalpreprova,sev_ph, by="IDParcdom")%>% unique()

DPhalpreprova <- cbind(DPhalpreprova,scale(DPhalpreprova[,c(2:21,29:30,42)]))
colnames(DPhalpreprova)[43:ncol(DPhalpreprova)] <- paste0(colnames(DPhalpreprova)[43:ncol(DPhalpreprova)],'_sc')


###base datos post


DPhalpostprova <- ddply(DistPhalpost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
#pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
#pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPhalpostprova$mindpl5 <- apply(DPhalpostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DPhalpostprova$minwpl <- apply(DPhalpostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

DPhalpostprova <- cbind(DPhalpostprova,scale(DPhalpostprova[,c(2:21,29:30)]))
colnames(DPhalpostprova)[31:ncol(DPhalpostprova)] <- paste0(colnames(DPhalpostprova)[31:ncol(DPhalpostprova)],'_sc')

###union pre y post

DPhalprova <-dplyr::full_join(DPhalpreprova,DPhalpostprova, by="IDParcdom")%>% unique()

DPhalprova$intercambio.x[DPhalprova$intercambio.x != "0"] <- "1"
DPhalprova$intercambio.x <-as.numeric(DPhalprova$intercambio.x)




write.csv2(DPhalprova,"C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Phalmodel.csv")



mphp1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x+wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y++ wet_sc.y+wetm_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DPhalprova)
mphp2<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x+wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y++ wet_sc.y+wetm_sc.y+ diffire3.x +severity_sc+(1|ID_fire.x) ,DPhalprova)
mphp3<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x+wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y++ wet_sc.y+wetm_sc.y+ diffire3.x +severity_sc+(1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DPhalprova)

summ(mphp3)
pph1 <- dwplot(mphp1) + ggtitle(paste('progresion PH     ','R2m =',round(unlist(r.squaredGLMM(mphp1))[1],2),'R2c =',round(unlist(r.squaredGLMM(mphp1))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph2 <- dwplot(mphp2) + ggtitle(paste('sucesion PH   ','R2m =',round(unlist(r.squaredGLMM(mphp2))[1],2),'R2c =',round(unlist(r.squaredGLMM(mphp2))[2],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph3 <- dwplot(mphp3) + ggtitle(paste('intercambio PH   ','R2m =',round(unlist(r.squaredGLMM(mphp3))[1],2),'R2c =',round(unlist(r.squaredGLMM(mphp3))[2],2))) + scale_color_manual(values=c("orchid")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/dryprovaph.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pph1,pph2,pph3)
dev.off()





###PINUS NIGRA####

ph_inter <-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/EUF_Pnigra_PCA12Bio.csv")
hpi_inter <- Hscv (x=as.matrix(ph_inter[, c(13,14)]),  binned=T, pilot="samse")# estimating kernel bandwidth with crossvalidation, it can be also selected manually
niche_volume_inter <- kde(x=ph_inter[, c(13,14)], gridsize =  2500, H=hpi_inter)# weighted by specie' abundance #

niche_rast_inter <- raster(niche_volume_inter)# do that in the complex way (niche_volume$estimate, xmx, etc) give a especular plot
niche_rast_inter@data@values[which(niche_rast_inter@data@values<niche_volume_inter$cont[95])] <- 0# here select the desired percentile
niche_inter_spdf <- as(niche_rast_inter, "SpatialPixelsDataFrame")
niche_inter_df <- data.frame(niche_inter_spdf@coords, niche_inter_spdf@data)
niche_inter_df[niche_inter_df$layer< contourLevels(niche_volume_inter, prob=0.05), "layer"] <- 0

cograv <- COGravity(x=niche_inter_df[, 1], y=niche_inter_df[, 2], 
                    wt=niche_inter_df[, 3]^40)# it could also work with lower exponent
cograv# option a centre of mass
centroid <- data.frame(matrix(nrow=1, ncol=2))
centroid[, 1] <- cograv[[1]]
centroid[, 2] <- cograv[[3]]
## coordenades centroide
centroid
xcpn <- 1.236544
ycpn <- -0.1479887

  niche_rast1 <- niche_rast_inter
niche_rast1@data@values[which(is.na(niche_rast1@data@values))] <- 0
niche_rast1@data@values[niche_rast1@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
a <- rasterToContour(niche_rast1)
b <-  as(a, "SpatialPointsDataFrame")
b <- data.frame(b@coords)
names(b) <- c("pca1", "pca2")

limitpn <-b

#datbase con los ejes
DPnigpreprova <- ddply(DistPnigpre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire),
                       pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                       pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPnigpostprova <- ddply(DistPnigpost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                        coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                        edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                        ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                        sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire),
                        pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                        pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                        pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                        pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

## selecció//definició del wet edge (queno s'ha fet servir) i el dry edge
wet_points_limit <- rbind(subset(limitpn,pca1<(-4.5))%>%subset(.,pca2<3.2),subset(limitpn,pca1>-4.9)%>%subset(.,pca2>2))
## !!! dry limit només agafant la projecció de bio04 i la de bio6 ~~ aproximadament a 0.51 de l'eix 1 del PCA
dry_points_limit5 <- rbind(subset(limitpn,pca1>1)%>%subset(.,pca2<4.2),subset(limitpn,pca1>2)%>%subset(.,pca2>-3.5))








wpl <- wet_points_limit[,1:2]
dpl5 <- dry_points_limit5[,1:2]
colnames(wpl) <- c('pc1l','pc2l')
colnames(dpl5) <- c('pc1l','pc2l')

##distance to the wet and dry limits of all observations pre fire
DistPnigpre$wetdist <- apply(DistPnigpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistPnigpre$drydist5 <- apply(DistPnigpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistPnigpre$fedge <- apply(DistPnigpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitpn[,1:2],method='Euclidean')))
DistPnigpre$fcore <- apply(DistPnigpre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcpn,ycpn),ncol=2)),method='Euclidean')))

##distance to the wet and dry limits of all observations post fire
DistPnigpost$wetdist <- apply(DistPnigpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistPnigpost$drydist5 <- apply(DistPnigpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistPnigpost$fedge <- apply(DistPnigpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitpn[,1:2],method='Euclidean')))
DistPnigpost$fcore <- apply(DistPnigpost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcpn,ycpn),ncol=2)),method='Euclidean')))

####sACAR LAS OCURRENCIAS FUERA DE LIMITE Y PONDERAR LAS DISTANCIAS BORDE SECO
###PREFIRE
inoutPn <- in.kde(x=DistPnigpre[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutPn)
DistPnigpre <- cbind(DistPnigpre, inoutPn)

DistPnigpre$DistWDE
indice_false <- which(inoutPn == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistPnigpre[k,'DistWDE'] <- DistPnigpre[k,"min.dist.cen"] / (DistPnigpre[k,"min.dist.cen"] - DistPnigpre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutPn == TRUE)
for(k in indice_true){
  
  DistPnigpre[k,'DistWDE'] <- DistPnigpre[k,"min.dist.cen"] / (DistPnigpre[k,"min.dist.cen"] + DistPnigpre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}
###POSTFIRE
inoutPn <- in.kde(x=DistPnigpost[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutPh)
DistPnigpost <- cbind(DistPnigpost, inoutPn)

DistPnigpost$DistWDE
indice_false <- which(inoutPn == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistPnigpost[k,'DistWDE'] <- DistPnigpost[k,"min.dist.cen"] / (DistPnigpost[k,"min.dist.cen"] - DistPnigpost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutPn == TRUE)
for(k in indice_true){
  
  DistPnigpost[k,'DistWDE'] <- DistPnigpost[k,"min.dist.cen"] / (DistPnigpost[k,"min.dist.cen"] + DistPnigpost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}


##base datos pre y severidad
DPnigpreprova <- ddply(DistPnigpre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
#pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
#pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPnigpreprova$mindpl5 <- apply(DPnigpreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DPnigpreprova$minwpl <- apply(DPnigpreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

sev_pn<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/DataFirePnig_Severity.csv")

DPnigpreprova <-dplyr::full_join(DPnigpreprova,sev_pn, by="IDParcdom")%>% unique()

DPnigpreprova <- cbind(DPnigpreprova,scale(DPnigpreprova[,c(2:21,29:30,40)]))
colnames(DPnigpreprova)[41:ncol(DPnigpreprova)] <- paste0(colnames(DPnigpreprova)[41:ncol(DPnigpreprova)],'_sc')

###base datos post


DPnigpostprova <- ddply(DistPnigpost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                        coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                        edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                        xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                        dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                        drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                        drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                        WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                        ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                        sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
#pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
#pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DPnigpostprova$mindpl5 <- apply(DPnigpostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DPnigpostprova$minwpl <- apply(DPnigpostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

DPnigpostprova <- cbind(DPnigpostprova,scale(DPnigpostprova[,c(2:21,29:30)]))
colnames(DPnigpostprova)[31:ncol(DPnigpostprova)] <- paste0(colnames(DPnigpostprova)[31:ncol(DPnigpostprova)],'_sc')

###union pre y post

DPnigprova <-dplyr::full_join(DPnigpreprova,DPnigpostprova, by="IDParcdom")%>% unique()

DPnigprova$intercambio.x[DPnigprova$intercambio.x != "0"] <- "1"
DPnigprova$intercambio.x <-as.numeric(DPnigprova$intercambio.x)


write.csv2(DPnigprova,"C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Pnigmodel.csv")


mpnp1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) ,DPnigprova)
mpnp2<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) ,DPnigprova)
mpnp3<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DPnigprova)

summ(mpnp3)
pph1 <- dwplot(mpnp1) + ggtitle(paste('progresion PN     ','R2m =',round(unlist(r.squaredGLMM(mpnp1))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp1))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph2 <- dwplot(mpnp2) + ggtitle(paste('sucesion PN   ','R2m =',round(unlist(r.squaredGLMM(mpnp2))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp2))[2],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph3 <- dwplot(mpnp3) + ggtitle(paste('intercambio PN   ','R2m =',round(unlist(r.squaredGLMM(mpnp3))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp3))[2],2))) + scale_color_manual(values=c("orchid")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/dryprovapn.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pph1,pph2,pph3)
dev.off()

###QUERCUS ILEX####
ph_inter <-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/EUF_Qilex_Bio_PCA12Bio.csv")
hpi_inter <- Hscv (x=as.matrix(ph_inter[, c(13,14)]),  binned=T, pilot="samse")# estimating kernel bandwidth with crossvalidation, it can be also selected manually
niche_volume_inter <- kde(x=ph_inter[, c(13,14)], gridsize =  2500, H=hpi_inter)# weighted by specie' abundance #
niche_rast_inter <- raster(niche_volume_inter)# do that in the complex way (niche_volume$estimate, xmx, etc) give a especular plot
niche_rast_inter@data@values[which(niche_rast_inter@data@values<niche_volume_inter$cont[95])] <- 0# here select the desired percentile
niche_inter_spdf <- as(niche_rast_inter, "SpatialPixelsDataFrame")
niche_inter_df <- data.frame(niche_inter_spdf@coords, niche_inter_spdf@data)
niche_inter_df[niche_inter_df$layer< contourLevels(niche_volume_inter, prob=0.05), "layer"] <- 0

cograv <- COGravity(x=niche_inter_df[, 1], y=niche_inter_df[, 2], 
                    wt=niche_inter_df[, 3]^40)# it could also work with lower exponent
cograv# option a centre of mass
centroid <- data.frame(matrix(nrow=1, ncol=2))
centroid[, 1] <- cograv[[1]]
centroid[, 2] <- cograv[[3]]
## coordenades centroide
centroid
xcqi <- -0.7851855
ycqi <-  1.273882

niche_rast1 <- niche_rast_inter
niche_rast1@data@values[which(is.na(niche_rast1@data@values))] <- 0
niche_rast1@data@values[niche_rast1@data@values>0] <- 1 # obtain raster with continuous values to get only one contour
a <- rasterToContour(niche_rast1)
b <-  as(a, "SpatialPointsDataFrame")
b <- data.frame(b@coords)
names(b) <- c("pca1", "pca2")

limitqi <-b

#datbase con los ejes
DQilepreprova <- ddply(DistQilepre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
                       #pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       #pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                       #pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                       #pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DQilepostprova <- ddply(DistQilepost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                        coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                        edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                        ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                        sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
                        #pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                        #pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
                        #pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
                        #pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

## selecció//definició del wet edge (queno s'ha fet servir) i el dry edge
wet_points_limit <- rbind(subset(limitph,pca1<0)%>%subset(.,pca2>(-3)),subset(limitph,pca1>-5.5)%>%subset(.,pca2<(-2.5)))

## !!! dry limit només agafant la projecció de bio04 i la de bio6 ~~ aproximadament a 0.51 de l'eix 1 del PCA
dry_points_limit5 <- rbind(subset(limitph,pca1>1.5)%>%subset(.,pca2<3),subset(limitph,pca1>4)%>%subset(.,pca2<(-2.5)))

wpl <- wet_points_limit[,1:2]
dpl5 <- dry_points_limit5[,1:2]
colnames(wpl) <- c('pc1l','pc2l')
colnames(dpl5) <- c('pc1l','pc2l')

##distance to the wet and dry limits of all observations pre fire
DistQilepre$wetdist <- apply(DistQilepre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistQilepre$drydist5 <- apply(DistQilepre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistQilepre$fedge <- apply(DistQilepre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitqi[,1:2],method='Euclidean')))
DistQilepre$fcore <- apply(DistQilepre[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcqi,ycqi),ncol=2)),method='Euclidean')))

##distance to the wet and dry limits of all observations post fire
DistQilepost$wetdist <- apply(DistQilepost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))
DistQilepost$drydist5 <- apply(DistQilepost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DistQilepost$fedge <- apply(DistQilepost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),limitqi[,1:2],method='Euclidean')))
DistQilepost$fcore <- apply(DistQilepost[,c('Axis1','Axis2')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),as.data.frame(matrix(c(xcqi,ycqi),ncol=2)),method='Euclidean')))

####sACAR LAS OCURRENCIAS FUERA DE LIMITE Y PONDERAR LAS DISTANCIAS BORDE SECO
###PREFIRE
inoutQi <- in.kde(x=DistQilepre[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutQi)
DistQilepre <- cbind(DistQilepre, inoutQi)

DistQilepre$DistWDE
indice_false <- which(inoutQi == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistQilepre[k,'DistWDE'] <- DistQilepre[k,"min.dist.cen"] / (DistQilepre[k,"min.dist.cen"] - DistQilepre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutQi == TRUE)
for(k in indice_true){
  
  DistQilepre[k,'DistWDE'] <- DistQilepre[k,"min.dist.cen"] / (DistQilepre[k,"min.dist.cen"] + DistQilepre[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}
###POSTFIRE
inoutQi <- in.kde(x=DistQilepost[,c('Axis1','Axis2')], fhat=niche_volume_inter,abs.cont=ks::contourLevels(niche_volume_inter, cont=99)) #False fuera, True dentro
table(inoutQi)
DistQilepost <- cbind(DistQilepost, inoutQi)

DistQilepost$DistWDE
indice_false <- which(inoutQi == FALSE)
k <- indice_false[1]
for(k in indice_false){
  # DistPhal[k,'DistPonderada'] <- a / (a-b)  ## a = dist to core ; b == dist edge #crear primero la columna de po derada
  
  DistQilepost[k,'DistWDE'] <- DistQilepost[k,"min.dist.cen"] / (DistQilepost[k,"min.dist.cen"] - DistQilepost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}

indice_true <- which(inoutQi == TRUE)
for(k in indice_true){
  
  DistQilepost[k,'DistWDE'] <- DistQilepost[k,"min.dist.cen"] / (DistQilepost[k,"min.dist.cen"] + DistQilepost[k,"drydist5"])  ## a = dist to core ; b == dist edge
  
}



##base datos pre y severidad
DQilepreprova <- ddply(DistQilepre,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                       coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                       edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                       xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                       dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                       drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                       drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                       WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                       ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                       sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
#pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
#pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DQilepreprova$mindpl5 <- apply(DQilepreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DQilepreprova$minwpl <- apply(DQilepreprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

sev_qi<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/DataFireQile_Severity.csv")

DQilepreprova <-dplyr::full_join(DQilepreprova,sev_qi, by="IDParcdom")%>% unique()

DQilepreprova <- cbind(DQilepreprova,scale(DQilepreprova[,c(2:21,29:30,40)]))
colnames(DQilepreprova)[41:ncol(DQilepreprova)] <- paste0(colnames(DQilepreprova)[41:ncol(DQilepreprova)],'_sc')

###base datos post


DQilepostprova <- ddply(DistQilepost,.(IDParcdom),summarize,pondmin = min(DistPonderada,na.rm=T), pondm = mean(DistPonderada,na.rm=T),  ponds = sd(DistPonderada,na.rm=T),
                        coremin = min(min.dist.cen,na.rm=T), corem = mean(min.dist.cen,na.rm=T),  cores = sd(min.dist.cen,na.rm=T), 
                        edgemin = min(min.dist.lim,na.rm=T), edgem = mean(min.dist.lim,na.rm=T),  edges = sd(min.dist.lim,na.rm=T),
                        xedge = Axis1[which(fedge == min(fedge,na.rm=T))],yedge = Axis2[which(fedge == min(fedge,na.rm=T))],
                        dry = min(drydist5,na.rm=T),wet=min(wetdist,na.rm=T),
                        drym = mean(drydist5,na.rm=T),wetm = mean(wetdist,na.rm=T),
                        drys = sd(drydist5,na.rm=T),wets = sd(wetdist,na.rm=T),
                        WeightedDEmin = min(DistWDE, na.rm=T), WeightedDEm = mean(DistWDE, na.rm=T), WeightedDEs = sd(DistWDE, na.rm=T),
                        ID_fire = mean(fire_id), sup_ha =mean(sup_ha), progresion=mean(progresion),
                        sucesion=mean(sucesion), intercambio=mean(intercambio),diffire3=mean(diffire3), difbiofire=mean(difbiofire))
#pc1c=Axis1[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc2c=Axis2[which(min.dist.cen==min(min.dist.cen,na.rm=T))],
#pc1l=Axis1[which(min.dist.lim==min(min.dist.lim,na.rm=T))],
#pc2l=Axis2[which(min.dist.lim==min(min.dist.lim,na.rm=T))])

DQilepostprova$mindpl5 <- apply(DQilepostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),dpl5,method='Euclidean')))
DQilepostprova$minwpl <- apply(DQilepostprova[,c('xedge','yedge')],1,function(x) min(proxy::dist(as.data.frame(matrix(x,ncol=2)),wpl,method='Euclidean')))

DQilepostprova <- cbind(DQilepostprova,scale(DQilepostprova[,c(2:21,29:30)]))
colnames(DQilepostprova)[31:ncol(DQilepostprova)] <- paste0(colnames(DQilepostprova)[31:ncol(DQilepostprova)],'_sc')

###union pre y post

DQileprova <-dplyr::full_join(DQilepreprova,DQilepostprova, by="IDParcdom")%>% unique()

DQileprova$intercambio.x[DQileprova$intercambio.x != "0"] <- "1"
DQileprova$intercambio.x <-as.numeric(DQileprova$intercambio.x)


write.csv2(DQileprova,"C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Qilemodel.csv")


mqip1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) ,DQileprova)
mqip2<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) ,DQileprova)
mqip3<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + edgemin_sc.x + edges_sc.x+ dry_sc.x+ drym_sc.x+ wet_sc.x + wetm_sc.x+ pondm_sc.y + ponds_sc.y + edgemin_sc.y + edges_sc.y+ dry_sc.y+ drym_sc.y+ wet_sc.y + wetm_sc.y+diffire3.x +severity_sc+ (1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DQileprova)

summ(mqip3)
pqi1 <- dwplot(mqip1) + ggtitle(paste('progresion QI     ','R2m =',round(unlist(r.squaredGLMM(mqip1))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip1))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pqi2 <- dwplot(mqip2) + ggtitle(paste('sucesion QI   ','R2m =',round(unlist(r.squaredGLMM(mqip2))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip2))[2],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pqi3 <- dwplot(mqip3) + ggtitle(paste('intercambio QI   ','R2m =',round(unlist(r.squaredGLMM(mqip3))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip3))[2],2))) + scale_color_manual(values=c("orchid")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/dryprovaqi.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pqi1,pqi2,pqi3)
dev.off()


#####pruebas MODELOS#####

DPhalprova1<-filter(DPhalprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion

mphp4<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DPhalprova1)
mphp5<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DPhalprova1)
mphp6<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+  (1|ID_fire.x), family=binomial,control = glmerControl(optimizer ="bobyqa"),  DPhalprova)
summ(mphp6)


pph4 <- dwplot(mphp4) + ggtitle(paste('Coverture PH     ','R2m =',round(unlist(r.squaredGLMM(mphp4))[1],2),'R2c =',round(unlist(r.squaredGLMM(mphp4))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph6 <- dwplot(mphp6) + ggtitle(paste('Exchange PH   ','R2m =',round(unlist(r.squaredGLMM(mphp6))[1],2),'R2c =',round(unlist(r.squaredGLMM(mphp6))[3],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/PHDEF.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pph6, pph4)
dev.off()


DPnigprova1<-filter(DPnigprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion

mpnp4<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DPnigprova1)
mpnp5<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x+ wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DPnigprova1)
mpnp6<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x+ wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+(1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DPnigprova)


ppn4 <- dwplot(mpnp4) + ggtitle(paste('Coverture PN     ','R2m =',round(unlist(r.squaredGLMM(mpnp4))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp4))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
ppn6 <- dwplot(mpnp6) + ggtitle(paste('Exchange PN   ','R2m =',round(unlist(r.squaredGLMM(mpnp6))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp6))[2],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')


tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/PNDEF.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(ppn6, ppn4)
dev.off()

DQileprova1<-filter(DQileprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion

mqip4<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+  (1|ID_fire.x) ,DQileprova1)
mqip5<-lmer(sucesion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+  (1|ID_fire.x) ,DQileprova1)
mqip6<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DQileprova)
tab_model(mqip6)

summ(mqip6)
pqi4 <- dwplot(mqip4) + ggtitle(paste('Coverture QI     ','R2m =',round(unlist(r.squaredGLMM(mqip4))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip4))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pqi6 <- dwplot(mqip6) + ggtitle(paste('Exchange QI   ','R2m =',round(unlist(r.squaredGLMM(mqip6))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip6))[2],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/QIDEF.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pqi6, pqi4)
dev.off()

###GRAFICA INTERCAMBIO
ggplot(DPhalprova, aes(x=as.factor(intercambio.x), y=dry_sc.x))+ geom_boxplot(fill="slateblue", alpha=0.2)+ xlab("intercambio")
ggplot(DPhalprova, aes(x=as.factor(intercambio.x), y=drym_sc.x))+ geom_boxplot(fill="slateblue", alpha=0.2)+ xlab("intercambio")
plot_model(mphp3, type="pred", terms="dry_sc.x [all]")

####Violin plots
sp<-"PH"
DPhalprova1$especies<-sp
sp<-"PN"
DPnigprova1$especies<-sp
sp<-"QI"
DQileprova1$especies<-sp

provaplot <- rbind(DPhalprova1,DPnigprova1,DQileprova1 )
library(viridis)
Q <- ggplot(provaplot, aes(x=especies, y=progresion.x, fill=especies)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(width=1) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE)

sp<-"PH"
DPhalprova$especies<-sp
sp<-"PN"
DPnigprova$especies<-sp
sp<-"QI"
DQileprova$especies<-sp

provaploti <- rbind(DPhalprova,DPnigprova,DQileprova )
provaploti <- dplyr::filter(provaploti, !is.na(intercambio.x)) 
provaploti %>%
  count(especies, intercambio.x) %>%
  group_by(especies) %>%
  mutate(n = n/sum(n) * 100) %>%
  ggplot() + aes(especies, n, fill = intercambio.x, label = paste0(round(n, 2), "%")) + 
  geom_col() +
  geom_text(position=position_stack(0.5))

###sacar especies para el suplementary material
uviph <- ddply(DistPhalpre,.(IDParcdom),summarize,especie = min(Especie,na.rm=T))
table(uviph$especie)

uvipn <- ddply(DistPnigpre,.(IDParcdom),summarize,especie = min(Especie,na.rm=T))
table(uvipn$especie)

uviqi <- ddply(DistQilepre,.(IDParcdom),summarize,especie = min(Especie,na.rm=T))
table(uviqi$especie)

prop.table(table(uviqi[c('especie')]))


DQileprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Qilemodel.csv")
DPnigprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Pnigmodel.csv")
DPhalprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Phalmodel.csv")
sp<-"PH"
DPhalprova$especies<-sp
sp<-"PN"
DPnigprova$especies<-sp
sp<-"QI"
DQileprova$especies<-sp

provaplotsev <- rbind(DPhalprova,DPnigprova,DQileprova )
library(viridis)
install.packages("viridisLite")
library(viridisLite)
Q <- ggplot(provaplotsev, aes(x=especies, y=severity, fill=especies)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(width=1) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE)
tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/severi.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(Q)
dev.off()


DQileprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Qilemodel.csv")
DPnigprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Pnigmodel.csv")
DPhalprova<-read.csv2("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/dryedge/Phalmodel.csv")
sp<-"PH"
DPhalprova$especies<-sp
sp<-"PN"
DPnigprova$especies<-sp
sp<-"QI"
DQileprova$especies<-sp

####correlaciones
install.packages("ggcorrplot")
library(ggcorrplot)
QInums<-DQileprova[,c("pondm_sc.x", "ponds_sc.x", "dry_sc.x", "drym_sc.x", "drys_sc.x", "wet_sc.x", "wetm_sc.x", "wets_sc.x", "pondm_sc.y", "ponds_sc.y","dry_sc.y","drym_sc.y","drys_sc.y","wet_sc.y","wetm_sc.y","wets_sc.y","diffire3.x","severity_sc","ID_fire.x")]
corr <- round(cor(QInums,use = "pairwise.complete.obs"), 2)
head(corr[, 1:6])
ggcorrplot(corr)



PNnums<-DPnigprova[,c("pondm_sc.x", "ponds_sc.x", "dry_sc.x", "drym_sc.x", "drys_sc.x", "wet_sc.x", "wetm_sc.x", "wets_sc.x", "pondm_sc.y", "ponds_sc.y","dry_sc.y","drym_sc.y","drys_sc.y","wet_sc.y","wetm_sc.y","wets_sc.y","diffire3.x","severity_sc","ID_fire.x")]
corr <- round(cor(PNnums,use = "pairwise.complete.obs"), 2)
head(corr[, 1:6])
ggcorrplot(corr)


PHnums<-DPhalprova[,c("pondm_sc.x", "ponds_sc.x", "dry_sc.x", "drym_sc.x", "drys_sc.x", "wet_sc.x", "wetm_sc.x", "wets_sc.x", "pondm_sc.y", "ponds_sc.y","dry_sc.y","drym_sc.y","drys_sc.y","wet_sc.y","wetm_sc.y","wets_sc.y","diffire3.x","severity_sc","ID_fire.x")]
corr <- round(cor(PHnums,use = "pairwise.complete.obs"), 2)
head(corr[, 1:6])
ggcorrplot(corr)

#####pruebas MODELOS CAMBIANDO EL GLMM POR GLM SI NO HAY MÁS DE UN IFN POR FUEGO#####
head(DPhalprova)
library(dplyr)
table(DQileprova$ID_fire.x)

tabla_observaciones <- table(DPhalprova1$ID_fire.x)
# Contar cuántos incendios tienen exactamente una observación
incendios_unica_observacion <- sum(tabla_observaciones == 1)
# Calcular el porcentaje
total_incedios <- length(tabla_observaciones)
porcentaje_unica_observacion <- (incendios_unica_observacion / total_incedios) * 100
porcentaje_unica_observacion



DPhalprova1<-filter(DPhalprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion
mmph1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+  (1|ID_fire.x) ,DPhalprova1)
mmph2<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+  (1|ID_fire.x), family=binomial,control = glmerControl(optimizer ="bobyqa"),  DPhalprova)

mphp7<-lm(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc ,DPhalprova1)
mphp8<-glm(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc, family=binomial,DPhalprova)
summ(mphp8)
tab_model(mmph2)
AIC(mmph2, mphp8)
###para p halepensis los modelos no mejoran
pph4 <- dwplot(mmph1) + ggtitle(paste('Coverture PH     ','R2m =',round(unlist(r.squaredGLMM(mmph1))[1],2),'R2c =',round(unlist(r.squaredGLMM(mmph1))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pph6 <- dwplot(mmph2) + ggtitle(paste('Exchange PH   ','R2m =',round(unlist(r.squaredGLMM(mmph2))[1],2),'R2c =',round(unlist(r.squaredGLMM(mmph2))[3],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/journal/PHAL.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pph6, pph4)
dev.off()


DPnigprova1<-filter(DPnigprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion

#mpnp7<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + drym_sc.x+ drys_sc.x+ wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + drym_sc.y+ drys_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc +(1|ID_fire.x),DPnigprova1)
#mpnp8<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + drym_sc.x+ drys_sc.x+ wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + drym_sc.y+ drys_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc +(1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"), DPnigprova)
mmpn1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc +(1|ID_fire.x),DPnigprova1, REML = FALSE)
mmpn2<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc + (1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"), DPnigprova)

mpnp7<-lm(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc,DPnigprova1)
mpnp8<-glm(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc , family=binomial,  DPnigprova)
summary(mpnp7)
BIC(mpnp8)

##sigue la singularidad, pero se va si cambio a glm en vez de glmer (no random effect)
isSingular(mpnp8, tol = 1e-4)
summary(mpnp7)
tab_model(mpnp7)


ppn4 <- dwplot(mpnp7) + ggtitle(paste('Coverture PN     ','R2m =',round(unlist(r.squaredGLMM(mpnp7))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp7))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
ppn6 <- dwplot(mpnp8) + ggtitle(paste('Exchange PN   ','R2m =',round(unlist(r.squaredGLMM(mpnp8))[1],2),'R2c =',round(unlist(r.squaredGLMM(mpnp8))[3],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')


tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/journal/PNIG.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(ppn6, ppn4)
dev.off()

DQileprova1<-filter(DQileprova, intercambio.x == 0)####sacar la base de datos en la que no ha habido intercambios para aplicar la sucesion y progresion

#mqip8<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + drym_sc.x+ drys_sc.x+ wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + drym_sc.y+ drys_sc.y + wetm_sc.y+ wets_sc.y + diffire3.x +severity_sc+(1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"),  DQileprova)

mmqi1<-lmer(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) ,DQileprova1)
mmqi2<-glmer(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc+ (1|ID_fire.x) , family=binomial,control = glmerControl(optimizer ="bobyqa"), DQileprova)

mqip7<-lm(progresion.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc,DQileprova1)
mqip8<-glm(intercambio.x ~ pondm_sc.x + ponds_sc.x + dry_sc.x+ drym_sc.x+ drys_sc.x+ wet_sc.x + wetm_sc.x+ wets_sc.x+ pondm_sc.y + ponds_sc.y + dry_sc.y+ drym_sc.y+ drys_sc.y +wet_sc.y + wetm_sc.y+ wets_sc.y+ diffire3.x +severity_sc, family=binomial, DQileprova)

BIC(mqip8)

##sigue la singularidad, pero se va si cambio a glm en vez de glmer (no random effect)
isSingular(mqip8, tol = 1e-4)
tab_model(mqip8)
summ(mqip6)

pqi4 <- dwplot(mmqi1) + ggtitle(paste('Coverture QI     ','R2m =',round(unlist(r.squaredGLMM(mmqi1))[1],2),'R2c =',round(unlist(r.squaredGLMM(mmqi1))[2],2))) + scale_color_manual(values=c("orange")) + geom_vline(xintercept=0,colour='black',linetype='dashed')
pqi6 <- dwplot(mqip8) + ggtitle(paste('Exchange QI   ','R2m =',round(unlist(r.squaredGLMM(mqip8))[1],2),'R2c =',round(unlist(r.squaredGLMM(mqip8))[3],2))) + scale_color_manual(values=c("green")) + geom_vline(xintercept=0,colour='black',linetype='dashed')

tiff("C:/Users/n.jimenez/OneDrive - CREAF/PROYECTO tesis/cap3/journal/QILE.tiff",width=2000,height=2500,res=350,compression='lzw')
grid.arrange(pqi6, pqi4)
dev.off()


