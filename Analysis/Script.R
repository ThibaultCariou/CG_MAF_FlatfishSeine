

load("Data/TabFin_estuseine95-19.RData")
load("Data/cotesseine.RData")
load("Data/polygone.RData")
coast <- rgdal::readOGR("Data/Coast_estuary_detailled.gpkg")
coast <- raster::crop(coast, raster::extent(-0.3,0.3,49.25,49.7))

#Libraries
library(RGeostats)
library(mapdata)
library(rgeos)
library(sp)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(kableExtra)
library(corrplot)
library(ggcorrplot)
library(ade4)
library(vegan)
library(ggalluvial)
source("Analysis/maf_petigas.R")
library(ggforce) #elispsis calcul
library(matlib) #ellipsis calcul
library(igraph) #networks
library(ggraph)
library(tidyr)
library(tidygraph)


#Data formating##### 
Pichon <- toutvcp %>% 
  filter(Nv_Rubbin %in% c("LIMALIM","PLATFLE","PLEUPLA","SOLEVUL")) %>% 
  dplyr::select(Code_Station, Nv_Rubbin, lonmoy, latmoy, annee, Nombre, Groupe, Taille, Division_p, Coeff_post_elev, Surface_chalutée)

#Filling and Homogeneisation of the group variable 
Taille_sp <- toutvcp %>% select(Nv_Rubbin, Groupe, Taille, annee) %>% 
  filter(Nv_Rubbin %in% c("LIMALIM","PLATFLE","PLEUPLA","SOLEVUL"))

Taille_sp.year <- Taille_sp %>% group_by(Nv_Rubbin, Groupe,annee) %>% summarise(Taille.moy=mean(Taille),
                                                                                Taille.min=min(Taille),
                                                                                Taille.max=max(Taille))

ggplot(Taille_sp)+
  geom_boxplot(aes(x= Nv_Rubbin, y=Taille, fill=Groupe  ))
#According to this, we redifine the group variable according to size 
Pichon$Groupe <- as.character(Pichon$Groupe)
Pichon$Groupe[Pichon$Nv_Rubbin=="LIMALIM" & Pichon$Taille<=10] <- "G0"
Pichon$Groupe[Pichon$Nv_Rubbin=="LIMALIM" & Pichon$Taille>10 & Pichon$Taille<=22] <- "G1"
Pichon$Groupe[Pichon$Nv_Rubbin=="LIMALIM" & Pichon$Taille>22] <- "G2p"

Pichon$Groupe[Pichon$Nv_Rubbin=="PLATFLE" & Pichon$Taille<=13] <- "G0"
Pichon$Groupe[Pichon$Nv_Rubbin=="PLATFLE" & Pichon$Taille>13 & Pichon$Taille<=20] <- "G1"
Pichon$Groupe[Pichon$Nv_Rubbin=="PLATFLE" & Pichon$Taille>20] <- "G2p"

Pichon$Groupe[Pichon$Nv_Rubbin=="PLEUPLA" & Pichon$Taille<=17] <- "G0"
Pichon$Groupe[Pichon$Nv_Rubbin=="PLEUPLA" & Pichon$Taille>17 & Pichon$Taille<=27] <- "G1"
Pichon$Groupe[Pichon$Nv_Rubbin=="PLEUPLA" & Pichon$Taille>27] <- "G2p"

Pichon$Groupe[Pichon$Nv_Rubbin=="SOLEVUL" & Pichon$Taille<=14] <- "G0"
Pichon$Groupe[Pichon$Nv_Rubbin=="SOLEVUL" & Pichon$Taille>14 & Pichon$Taille<=22] <- "G1"
Pichon$Groupe[Pichon$Nv_Rubbin=="SOLEVUL" & Pichon$Taille>22] <- "G2p"
Pichon$Groupe <- factor(Pichon$Groupe)

#Calcul of density for species age group
Pichon <- Pichon %>% na.omit()
Pichon <- Pichon %>% group_by(Code_Station, Nv_Rubbin, Groupe) %>% mutate(Nbre_gp= sum(Nombre)) %>% mutate(Densite= Nbre_gp*Division_p*Coeff_post_elev/Surface_chalutée)%>%ungroup()
Pichon <- Pichon %>% dplyr::select(-Coeff_post_elev, -Division_p)

#Density from ind per m2 to ind per ha
Pichon$Densite <- Pichon$Densite * 10000


#Creating table for points followed in time 
Pichon_dens <- Pichon
Pichon_dens <- left_join(Pichon_dens, toutvcp[,c("Numéro_du_Trait","Code_Station")], "Code_Station"="Code_Station")
Pichon_dens <- Pichon_dens[c(ncol(Pichon_dens),1:ncol(Pichon_dens)-1)]
Pichon_dens$Numéro_du_Trait[Pichon_dens$Code_Station==1] <- 1

#Creating corresponding table between hauls of different years
#corres <- read.csv("N:/03_RECHERCHES_ETUDES/17_20_NOURSEINE/These/RProjet/Test_Geostast/correspondance_numero_trait.csv", header = T, sep=";")
corres <- read.csv("Data/correspondance_numero_trait.csv", header = T, sep=";")

colnames(corres) <- c(1995:2002,2008:2010,2017:2019)
corres <- corres[,-1] #Delete 1995 because of less sampling point
corres <- na.omit(corres)



#Keeping data for selected hauls through the years without 1995
Pichon_dens_corres <- Pichon_dens[ which(Pichon_dens$annee==1996 & Pichon_dens$Numéro_du_Trait%in% corres[,1] | 
                                           Pichon_dens$annee==1997 & Pichon_dens$Numéro_du_Trait%in% corres[,2] | 
                                           Pichon_dens$annee==1998 & Pichon_dens$Numéro_du_Trait%in% corres[,3] | 
                                           Pichon_dens$annee==1999 & Pichon_dens$Numéro_du_Trait%in% corres[,4] | 
                                           Pichon_dens$annee==2000 & Pichon_dens$Numéro_du_Trait%in% corres[,5] | 
                                           Pichon_dens$annee==2001 & Pichon_dens$Numéro_du_Trait%in% corres[,6] | 
                                           Pichon_dens$annee==2002 & Pichon_dens$Numéro_du_Trait%in% corres[,7] | 
                                           Pichon_dens$annee==2008 & Pichon_dens$Numéro_du_Trait%in% corres[,8] | 
                                           Pichon_dens$annee==2009 & Pichon_dens$Numéro_du_Trait%in% corres[,9] | 
                                           Pichon_dens$annee==2010 & Pichon_dens$Numéro_du_Trait%in% corres[,10] | 
                                           Pichon_dens$annee==2017 & Pichon_dens$Numéro_du_Trait%in% corres[,11] | 
                                           Pichon_dens$annee==2018 & Pichon_dens$Numéro_du_Trait%in% corres[,12] |
                                           Pichon_dens$annee==2019 & Pichon_dens$Numéro_du_Trait%in% corres[,13]),]


corres$Correspondance <- seq(1:dim(corres)[1])
corres_long <- gather(corres, "annee","Numéro_du_Trait",1:13)
corres_long$annee <- as.numeric(corres_long$annee)


#Join table to get correspondance variable
Pichon_dens_corres <- left_join(Pichon_dens_corres,corres_long, by=c("annee","Numéro_du_Trait"))
Pichon_dens_corres <- Pichon_dens_corres %>% dplyr::select( -Taille, -Nombre,-Surface_chalutée, -Nbre_gp) 
Pichon_dens_corres <- Pichon_dens_corres %>% distinct()
#Creating mean coordinates for corresponding point
Pichon_dens_corres <- Pichon_dens_corres %>% group_by(Correspondance) %>% mutate(loncor=mean(lonmoy), latcor=mean(latmoy)) %>% ungroup()

#Adding the null density to table
toto <- Pichon_dens_corres
ttsel <- toto %>% distinct()
codif <- ttsel %>% dplyr::select(Numéro_du_Trait, Code_Station, annee, lonmoy, latmoy, Groupe, Correspondance, loncor,latcor) %>% distinct()
allstrat <- ttsel %>% expand(Code_Station, Nv_Rubbin)
pipo <- left_join(allstrat, codif) %>% left_join(ttsel) %>% mutate(Densite = ifelse(is.na(Densite), 0, Densite))
Pichon_dens_corres <- left_join(codif, pipo)


#Table basic stats#####
tab <- Pichon_dens_corres %>% dplyr::select(Code_Station, annee, Nv_Rubbin,Groupe,Densite)
tab <- tab %>% filter(annee != 1995)
tab <- unique(tab)
toto <- tab[tab$Densite!=0,]
tab1 <- toto %>% group_by(Nv_Rubbin,Groupe) %>%summarise(Mean=mean(Densite), CV=sd(Densite, na.rm=TRUE)/mean(Densite, na.rm=TRUE)*100)
tab1 <- toto %>% group_by(Nv_Rubbin,Groupe) %>% mutate(Found=length(unique(Code_Station))/length(unique(tab$Code_Station))*100) %>%
  dplyr::select(Nv_Rubbin,Groupe, Found) %>% left_join(tab1,"Nv_Rubbin"="Nv_Rubbin","Groupe"="Groupe") %>% distinct()
tab1
#####


#Figure of mean distribution#####
toto <- Pichon_dens_corres 
toto <- toto %>% group_by(loncor,latcor,Nv_Rubbin,Groupe) %>% summarise(Mean=mean(Densite))
toto <- toto %>% filter(Groupe=="G0") %>% ungroup()
toto$Sci.red <- NA
toto$Sci.red[toto$Nv_Rubbin=="LIMALIM"] <- "L.limanda"
toto$Sci.red[toto$Nv_Rubbin=="PLATFLE"] <- "P.flesus"
toto$Sci.red[toto$Nv_Rubbin=="PLEUPLA"] <- "P.platessa"
toto$Sci.red[toto$Nv_Rubbin=="SOLEVUL"] <- "S.solea"

toto$Common[toto$Nv_Rubbin=="LIMALIM"] <- "Dab"
toto$Common[toto$Nv_Rubbin=="PLATFLE"] <- "Flounder"
toto$Common[toto$Nv_Rubbin=="PLEUPLA"] <- "Plaice"
toto$Common[toto$Nv_Rubbin=="SOLEVUL"] <- "Sole"

toto <- toto %>% filter(Common!="Flounder") %>% ungroup()

mean_distri <- ggplot()+
  geom_point(data=toto[toto$Mean!=0,],aes(x = loncor, y = latcor, col=Mean, size=Mean), alpha=0.7) +  
  geom_point(data=toto[toto$Mean==0,],aes(x = loncor, y = latcor), size=3, pch=4) +  
  theme_bw() + scale_color_viridis()  + 
  ggtitle("") + xlab("Longitude (dec.)")+ylab("Latitude (dec.)") +
  facet_wrap(~Common)+
  geom_polygon(data = coast, aes(x = long, y = lat, group = group),fill="grey",col="black") + coord_fixed() +
  labs(size = "Mean density\n(ind. per hectare)",col = "Mean density\n(ind. per hectare)")
mean_distri

#ggsave("Figure_1.pdf",path="Results/Figures")
#ggsave("Figure_1.tiff",path="Results/Figures")
#####


#Gravity centers of flatfishes#####

#Data formating, spread species x group variable 
data.CG <- Pichon_dens_corres
data.CG$cat <- paste(data.CG$Nv_Rubbin,data.CG$Groupe, sep="_")
data.CG <- data.CG %>% distinct() %>% tidyr::spread(cat, Densite) %>%  dplyr::select(-Numéro_du_Trait, -Code_Station, -lonmoy, -latmoy, -Groupe, -Nv_Rubbin, -Correspondance) 

#Group by year
data.CG <- distinct(data.CG) #unique rows
data.CG <- as.data.frame(data.CG)
data.CG[which(is.na(data.CG),arr.ind = T)] <- 0
data.CG <- data.CG %>% group_by(loncor,latcor,annee) %>% summarise_all(sum)


noms <- data.frame(Rubbin=sort(unique(Pichon_dens_corres$Nv_Rubbin)), Sci=c("Limanda limanda","Platichthys flesus","Pleuronectes platessa","Solea solea"),
                   Species=c("L.limanda","P.flesus","P.platessa","S.solea"), Common=c("Dab","Flounder","Plaice","Sole"))


#Creation of db object for RGeostats
db.dens <- db.create(data.CG)
db.dens <- db.locate(db.dens,c("loncor","latcor"),"x")




#Initiating loop
tax_plats <- unique(Pichon_dens_corres$Nv_Rubbin)
Classe <- character() #Group of selected species
Annee <- numeric() #Year of selected species
Lon <- numeric() #Longitude of gravity center
Lat <- numeric() #Latitude of gravity center
iso <- numeric() #Isotropy index
Tax <- character() #Taxa selected
iner <- numeric()

x1 <- numeric() #Coordinates of the 2 inertia axes
x2 <- numeric()
x3 <- numeric()
x4 <- numeric()

y1 <- numeric()
y2 <- numeric()
y3 <- numeric()
y4 <- numeric()

#Loop
for(sp in tax_plats){
  tmp <- data.CG %>% select("loncor","latcor","annee", contains(sp)) #Group by hauls
  tmp <- tmp %>% group_by(loncor,latcor,annee) %>% summarise_all(sum)
  df <- tmp
  names(df) <- names(tmp)
  t<-0
  
  
  for (i in sort(unique(Pichon_dens_corres[Pichon_dens_corres$Nv_Rubbin==sp,]$Groupe))){ #Select group for selected sp
    t<-t+1
    z <- which(names(db.dens@items)== names(df)[-c(1:3)][t])
    db.dens <- db.locate(db.dens,z,"z")
    
    for (j in sort(unique(data.CG$annee))){ #Select year for selected sp and group
      
      if (colSums(df[df$annee==j,])[names(df)[-c(1:3)][t]]!=0){ #Loop if there are densities recorded
        toto <- db.sel(db.dens, annee==j)
        
        #Area of influnce for weighting factor
        toto <- infl(toto,nodes=c(400,400),origin=c(-0.3,49.2),extend=c(0.7,0.5),dmax=1,polygon=poly.data,plot=t,asp=1)  #recreate influence surface
        
        
        #Keep information of the loop index and results
        Classe <- c(Classe, i)
        Annee <- c(Annee, j)
        tmp <-c(SI.cgi(toto,flag.plot=F)$center[1],SI.cgi(toto,flag.plot=F)$center[2])
        Lon <- c(Lon, tmp[1])
        Lat <- c(Lat, tmp[2])
        iso <- c(iso, SI.cgi(toto,flag.plot=F)$iso)
        iner <- c(iner, SI.cgi(toto,flag.plot=F)$inertia)
        Tax <- c(Tax,sp)
        
        x1 <- c(x1,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[1,1])
        x2 <- c(x2,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[2,1])
        x3 <- c(x3,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[3,1])
        x4 <- c(x4,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[4,1])
        
        y1 <- c(y1,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[1,2])
        y2 <- c(y2,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[2,2])
        y3 <- c(y3,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[3,2])
        y4 <- c(y4,SI.cgi(toto,flag.plot=T, flag.inertia = T)$axes[4,2])
      }
      else{Classe <- c(Classe, i)
      Annee <- c(Annee, j)
      Lon <- c(Lon, NA)
      Lat <- c(Lat, NA)
      iso <- c(iso, NA)
      iner <- c(iner,NA)
      Tax <- c(Tax,sp)
      x1 <- c(x1,NA)
      x2 <- c(x2,NA)
      x3 <- c(x3,NA)
      x4 <- c(x4,NA)
      
      y1 <- c(y1,NA)
      y2 <- c(y2,NA)
      y3 <- c(y3,NA)
      y4 <- c(y4,NA)}
    }
  }
  CG <- data.frame(Classe=Classe, Annee=Annee,Lon=Lon,Lat=Lat, Iso=iso,Iner=iner, x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3,y4=y4)
  db.dens <- db.locate(db.dens,7,"z")
  
}

#Adding scientific names
CG_df <- cbind(CG,factor(Tax))
names(CG_df)[which(names(CG_df)=="factor(Tax)")] <- "Tax"
CG_df <- CG_df %>%
  left_join(noms, by = c("Tax" = "Rubbin"))


Mean_iso <- CG_df %>% dplyr::select(Common, Classe, x1,x2,x3,x4,y1,y2,y3,y4)%>% group_by(Common,Classe) %>% summarise(x1_mean=mean(x1, na.rm = T),
                                                                                                                        x2_mean=mean(x2, na.rm = T),
                                                                                                                        x3_mean=mean(x3, na.rm = T),
                                                                                                                        x4_mean=mean(x4, na.rm = T),
                                                                                                                        y1_mean=mean(y1, na.rm = T),
                                                                                                                        y2_mean=mean(y2, na.rm = T),
                                                                                                                        y3_mean=mean(y3, na.rm = T),
                                                                                                                        y4_mean=mean(y4, na.rm = T))
Mean_iso <- data.frame(Common=rep(Mean_iso$Common,2), Classe=rep(Mean_iso$Classe,2), x=c(Mean_iso$x1_mean,Mean_iso$x3_mean),y=c(Mean_iso$y1_mean,Mean_iso$y3_mean),xend=c(Mean_iso$x2_mean,Mean_iso$x4_mean), yend=c(Mean_iso$y2_mean,Mean_iso$y4_mean))

#####

#Figure inertia/isotropy#####
CG_dfG0 <- CG_df[CG_df$Classe=="G0",]
CG_dfG0 <- CG_dfG0[CG_dfG0$Common!="Flounder",]
Mean_iso <- Mean_iso[Mean_iso$Common!="Flounder",]
names(Mean_iso)[1] <- "Species"
names(CG_dfG0)[17] <- "Sci.red"
names(CG_dfG0)[18] <- "Species"
coastCG <- raster::crop(coast, raster::extent(-0.15,0.3,49.2,49.60))

#↑1st version, with mean isotropy
Grav <- ggplot(data=CG_dfG0,aes(x=Lon, y=Lat)) + 
  geom_point(aes(col=Species),alpha=0.7, pch=19, size=2) +
  geom_segment(data=Mean_iso[Mean_iso$Classe=="G0",],aes(x=x,xend=xend,y=y, yend=yend, col=Species), lwd=1.3)+ 
  geom_path(data=cotes,aes(x=long, y=lat, group=group)) +
  coord_fixed()  + theme_bw() + xlab("Longitude (dec.)") + ylab("Latitude (dec.)") + ggtitle("") +
  xlim(-0.15,0.3) + ylim(49.32,49.60) +
  stat_ellipse(aes(color=Species)) +
  scale_colour_viridis_d() + guides(col=F)

#1st version
Iso <- ggplot(data=CG_dfG0) +
  #geom_violin(aes(x=Species, y=Iso, fill=Species),alpha=0.7) + 
  geom_boxplot(aes(x=Species, y=Iso),width=0.1)+ 
  geom_jitter(shape=21, position=position_jitter(0.2),aes(size=Iner,x=Species,y=Iso,fill=Species),alpha=0.5,col="black")+
  scale_fill_viridis_d() + scale_color_viridis_d() +theme_minimal() + ylab("Isotropy")

plo2 <- ggplot(data=CG_dfG0)+
  geom_point(aes(x=Annee,y=Iso,fill=Species,col=Species,size=Iner),alpha=0.8)+
  facet_wrap(.~Species,nrow=1)+ geom_boxplot(aes(y=Iso,x=2008),width=2)+
  scale_fill_viridis_d() + scale_color_viridis_d() +theme_bw() + ylab("Isotropy") +xlab("Year")+
  theme(axis.text.x = element_text(size=10, angle=45,hjust = 1))+
  scale_x_continuous(breaks=seq(1996,2020,4))+ labs(size="Inertia")

#Plot of all inertia axis
Inertia <- ggplot(data=CG_dfG0) + 
  geom_segment(aes(x=x1,xend=x2,y=y1, yend=y2, col=Species), lwd=1.3)+ 
  geom_segment(aes(x=x3,xend=x4,y=y3, yend=y4, col=Species), lwd=1.3)+ 
  geom_path(data=cotes,aes(x=long, y=lat, group=group)) + theme_bw() + xlab("Longitude (dec.)") + ylab("Latitude (dec.)") + ggtitle("") +
  xlim(-0.2,0.3) + ylim(49.2,49.65) +
  scale_colour_viridis_d() + guides(col=F)

#cowplot::plot_grid(Grav,Iso, nrow=1)


#Violin boxplot of inertia
Iner <- ggplot(data=CG_dfG0) +
  geom_violin(aes(x=Species, y=Iner, fill=Species),alpha=0.7) + 
  geom_boxplot(aes(x=Species, y=Iner),width=0.1)+ 
  scale_fill_viridis_d() + theme_minimal() + ylab("Inertia")



radius <- data.frame(Lon=CG_dfG0$Lon,Lat=CG_dfG0$Lat,
                     a=sqrt((CG_dfG0$x1-CG_dfG0$x2)^2+(CG_dfG0$y1-CG_dfG0$y2)^2),
                     b=sqrt((CG_dfG0$x3-CG_dfG0$x4)^2+(CG_dfG0$y3-CG_dfG0$y4)^2),
                     angle=asin(CG_dfG0$y1-CG_dfG0$y2)/sqrt((CG_dfG0$x1-CG_dfG0$x2)^2+(CG_dfG0$y1-CG_dfG0$y2)^2),
                     Species=CG_dfG0$Species,Annee=CG_dfG0$Annee)


plo1 <-ggplot(data=CG_dfG0) + 
  geom_point(aes(x=Lon,y=Lat,col=Species),alpha=0.8, pch=19, size=2)+
  #geom_segment(aes(x=x1,xend=x2,y=y1, yend=y2, col=Species), lwd=1.3)+ 
  #geom_segment(aes(x=x3,xend=x4,y=y3, yend=y4, col=Species), lwd=1.3)+ 
  geom_ellipse(data=radius,aes(x0=Lon,y0=Lat,a=a/6,b=b/6,angle=angle,col=Species,fill=Species))+
  geom_polygon(data = coastCG, aes(x = long, y = lat, group = group),fill="grey",col="black")+
  theme_bw() + xlab("Longitude (dec.)") + ylab("Latitude (dec.)") + ggtitle("") +
  xlim(-0.15,0.3) + ylim(49.25,49.6) +
  scale_colour_viridis_d() + guides(col=F,fill=F) +scale_fill_viridis_d(alpha=0.3)+coord_fixed()


cowplot::plot_grid(plo1,plo2, nrow=1)

#ggsave("Figure_2.pdf",path="Results/Figures")
#ggsave("Figure_2.tiff",path="Results/Figures")
#####

#Global index of collocation and local index of collocation#####
#Data
dat.IC <- Pichon_dens_corres
dat.IC <- dat.IC[dat.IC$Groupe=="G0",] #G0 selection
dat.IC$var <- paste(dat.IC$Nv_Rubbin,dat.IC$Groupe,dat.IC$annee, sep="_")

#Deleting useless variable and sreading variables
dat.IC <- dat.IC %>% dplyr::select(-Numéro_du_Trait, -Code_Station, -lonmoy, -latmoy, -Nv_Rubbin, -annee, -Groupe)
dat.IC <- spread(dat.IC,var, Densite)
dat.IC <- as.data.frame(dat.IC)
dat.IC[which(is.na(dat.IC),arr.ind = T)] <- 0


#Creation of db object for RGeostats
loc.data <- db.create(dat.IC)
loc.data <- db.locate(loc.data,c("loncor","latcor"),"x")
loc.data <- db.locate(loc.data,"LIMALIM_G0_1996","z")

#Area of influnce for weighting factor
loc.data <- infl(loc.data,nodes=c(400,400),origin=c(-0.3,49.2),extend=c(0.7,0.5),dmax=1,polygon=poly.data,plot=T,asp=1)  #recreate influence surface


#GIC calculation between years

#GIC dab
name <- paste0(substr(names(loc.data@items)[which(grepl("LIMALIM_G0",names(loc.data@items)))],12,16))
GIC.dab.G0 <- matrix(0, 
                     ncol = length(which(grepl("LIMALIM_G0",names(loc.data@items)))), 
                     nrow =  length(which(grepl("LIMALIM_G0",names(loc.data@items)))), 
                     dimnames = list(name,name))

for (i in 1:dim(GIC.dab.G0)[1]){
  for (j in 1:dim(GIC.dab.G0)[1]){
    t <- names(loc.data@items)[which(grepl("LIMALIM_G0",names(loc.data@items)))]
    GIC.dab.G0[i,j] <- SI.gic(loc.data,loc.data,name1=t[i],name2=t[j], flag.plot = F)
  }
  
}

#GIC flounder
name <- paste0(substr(names(loc.data@items)[which(grepl("PLATFLE_G0",names(loc.data@items)))],12,16))
GIC.flounder.G0 <- matrix(0, 
                          ncol = length(which(grepl("PLATFLE_G0",names(loc.data@items)))), 
                          nrow =  length(which(grepl("PLATFLE_G0",names(loc.data@items)))), 
                          dimnames = list(name,name))

for (i in 1:(dim(GIC.flounder.G0)[1]-2)){
  for (j in 1:(dim(GIC.flounder.G0)[1]-2)){
    t <- names(loc.data@items)[which(grepl("PLATFLE_G0",names(loc.data@items)))]
    t <- t[which(colSums(loc.data@items[which(grepl("PLATFLE_G0",names(loc.data@items)))])>0)]
    GIC.flounder.G0[substr(t[i],12,16),substr(t[j],12,16)] <- 
      SI.gic(loc.data,loc.data,name1=t[i],name2=t[j], flag.plot = F)
  }
  
}
GIC.flounder.G0[which(is.nan(GIC.flounder.G0),arr.ind=T)] <- 1
GIC.flounder.G0[c(1,9),] <- NA
GIC.flounder.G0[,c(1,9)] <- NA


#GIC plaice
name <- paste0(substr(names(loc.data@items)[which(grepl("PLEUPLA_G0",names(loc.data@items)))],12,16))
GIC.plaice.G0 <- matrix(0, 
                        ncol = length(which(grepl("PLEUPLA_G0",names(loc.data@items)))), 
                        nrow =  length(which(grepl("PLEUPLA_G0",names(loc.data@items)))), 
                        dimnames = list(name,name))

for (i in 1:dim(GIC.plaice.G0)[1]){
  for (j in 1:dim(GIC.plaice.G0)[1]){
    t <- names(loc.data@items)[which(grepl("PLEUPLA_G0",names(loc.data@items)))]
    GIC.plaice.G0[i,j] <- 
      SI.gic(loc.data,loc.data,name1=t[i],name2=t[j], flag.plot = F)
  }
  
}


#GIC sole
name <- paste0(substr(names(loc.data@items)[which(grepl("SOLEVUL_G0",names(loc.data@items)))],12,16))
GIC.sole.G0 <- matrix(0, 
                      ncol = length(which(grepl("SOLEVUL_G0",names(loc.data@items)))), 
                      nrow =  length(which(grepl("SOLEVUL_G0",names(loc.data@items)))), 
                      dimnames = list(name,name))

for (i in 1:dim(GIC.sole.G0)[1]){
  for (j in 1:dim(GIC.sole.G0)[1]){
    t <- names(loc.data@items)[which(grepl("SOLEVUL_G0",names(loc.data@items)))]
    GIC.sole.G0[i,j] <- 
      SI.gic(loc.data,loc.data,name1=t[i],name2=t[j], flag.plot = F)
  }
  
}



#LIC calculation

#LIC dab
name <- paste0(substr(names(loc.data@items)[which(grepl("LIMALIM_G0",names(loc.data@items)))],12,16))
LIC.dab.G0 <- matrix(0, 
                     ncol = length(which(grepl("LIMALIM_G0",names(loc.data@items)))), 
                     nrow =  length(which(grepl("LIMALIM_G0",names(loc.data@items)))), 
                     dimnames = list(name,name))

for (i in 1:dim(LIC.dab.G0)[1]){
  for (j in 1:dim(LIC.dab.G0)[1]){
    t <- names(loc.data@items)[which(grepl("LIMALIM_G0",names(loc.data@items)))]
    LIC.dab.G0[i,j] <- 
      SI.lic(loc.data,name1=t[i],name2=t[j])
  }
  
}
LIC.dab.G0[which(is.nan(LIC.dab.G0),arr.ind=T)] <- 1


#LIC flounder
name <- paste0(substr(names(loc.data@items)[which(grepl("PLATFLE_G0",names(loc.data@items)))],12,16))
LIC.flounder.G0 <- matrix(0, 
                          ncol = length(which(grepl("PLATFLE_G0",names(loc.data@items)))), 
                          nrow =  length(which(grepl("PLATFLE_G0",names(loc.data@items)))), 
                          dimnames = list(name,name))

for (i in 1:dim(LIC.flounder.G0)[1]){
  for (j in 1:dim(LIC.flounder.G0)[1]){
    t <- names(loc.data@items)[which(grepl("PLATFLE_G0",names(loc.data@items)))]
    LIC.flounder.G0[i,j] <- 
      SI.lic(loc.data,name1=t[i],name2=t[j])
  }
  
}
LIC.flounder.G0[which(is.nan(LIC.flounder.G0),arr.ind=T)] <- 1
LIC.flounder.G0[c(1,9),] <- NA
LIC.flounder.G0[,c(1,9)] <- NA

#LIC plaice
##G0
name <- paste0(substr(names(loc.data@items)[which(grepl("PLEUPLA_G0",names(loc.data@items)))],12,16))
LIC.plaice.G0 <- matrix(0, 
                        ncol = length(which(grepl("PLEUPLA_G0",names(loc.data@items)))), 
                        nrow =  length(which(grepl("PLEUPLA_G0",names(loc.data@items)))), 
                        dimnames = list(name,name))

for (i in 1:dim(LIC.plaice.G0)[1]){
  for (j in 1:dim(LIC.plaice.G0)[1]){
    t <- names(loc.data@items)[which(grepl("PLEUPLA_G0",names(loc.data@items)))]
    LIC.plaice.G0[i,j] <- 
      SI.lic(loc.data,name1=t[i],name2=t[j])
  }
  
}
LIC.plaice.G0[which(is.nan(LIC.plaice.G0),arr.ind=T)] <- 1


#LIC sole
name <- paste0(substr(names(loc.data@items)[which(grepl("SOLEVUL_G0",names(loc.data@items)))],12,16))
LIC.sole.G0 <- matrix(0, 
                      ncol = length(which(grepl("SOLEVUL_G0",names(loc.data@items)))), 
                      nrow =  length(which(grepl("SOLEVUL_G0",names(loc.data@items)))), 
                      dimnames = list(name,name))

for (i in 1:dim(LIC.sole.G0)[1]){
  for (j in 1:dim(LIC.sole.G0)[1]){
    t <- names(loc.data@items)[which(grepl("SOLEVUL_G0",names(loc.data@items)))]
    LIC.sole.G0[i,j] <- 
      SI.lic(loc.data,name1=t[i],name2=t[j])
  }
  
}
LIC.sole.G0[which(is.nan(LIC.sole.G0),arr.ind=T)] <- 1

#####

#Figure GIC/LIC#####
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  



#Correlation matrix
#Dab
index.dab <- matrix(NA,nrow=dim(GIC.dab.G0)[1],ncol=dim(GIC.dab.G0)[2],
                    dimnames = list(letters[1:13],letters[1:13]))
index.dab[lower.tri(index.dab)] <- GIC.dab.G0[lower.tri(GIC.dab.G0)]
index.dab[upper.tri(index.dab)] <- LIC.dab.G0[upper.tri(LIC.dab.G0)]

#Combining GIC/LIC index in matrix plot
cordab <- ggcorrplot(index.dab,
           method="square",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"),
           show.diag = F,
           lab=F)+
  #scale_fill_gradientn(colours=rev(col2(200)),limit=c(0, 1))+
  scale_fill_viridis()+
  scale_x_discrete(breaks=letters[1:13],
                     labels=colnames(GIC.dab.G0))+
  scale_y_discrete(breaks=letters[1:13],
                   labels=colnames(GIC.dab.G0))+theme_bw()+
  xlab("GIC")+ylab("LIC")+ggtitle("Dab") + labs(fill="Index value")+
  theme(axis.text.x = element_text(angle = 90)) + guides(fill=F)

  

#Plaice
index.plaice <- matrix(NA,nrow=dim(GIC.plaice.G0)[1],ncol=dim(GIC.plaice.G0)[2],
                    dimnames = list(letters[1:13],letters[1:13]))
index.plaice[lower.tri(index.plaice)] <- GIC.plaice.G0[lower.tri(GIC.plaice.G0)]
index.plaice[upper.tri(index.plaice)] <- LIC.plaice.G0[upper.tri(LIC.plaice.G0)]


corplaice <- ggcorrplot(index.plaice,
                     method = "square",
                     outline.col = "white",
                     ggtheme = ggplot2::theme_gray,
                     colors = c("#6D9EC1", "white", "#E46726"),
                     show.diag = F,
                     lab=F)+
  #scale_fill_gradientn(colours=rev(col2(200)),limit=c(0, 1))+
  scale_fill_viridis()+
  scale_x_discrete(breaks=letters[1:13],
                   labels=colnames(GIC.sole.G0))+
  scale_y_discrete(breaks=letters[1:13],
                   labels=colnames(GIC.sole.G0))+theme_bw()+
  xlab("GIC")+ylab("LIC")+ggtitle("Plaice")+ labs(fill="Index value")+
  theme(axis.text.x = element_text(angle = 90))+ guides(fill=F)

#Sole
index.sole <- matrix(NA,nrow=dim(GIC.sole.G0)[1],ncol=dim(GIC.sole.G0)[2],
                     dimnames = list(letters[1:13],letters[1:13]))
index.sole[lower.tri(index.sole)] <- GIC.sole.G0[lower.tri(GIC.sole.G0)]
index.sole[upper.tri(index.sole)] <- LIC.sole.G0[upper.tri(LIC.sole.G0)]


corsole <- ggcorrplot(index.sole,
                      method="square",
                      outline.col = "white",
                      ggtheme = ggplot2::theme_gray,
                      colors = c("#6D9EC1", "white", "#E46726"),
                      show.diag = F,
                      lab=F)+
  #scale_fill_gradientn(colours=rev(col2(200)),limit=c(0, 1))+
  scale_fill_viridis()+
  scale_x_discrete(breaks=letters[1:13],
                   labels=colnames(GIC.sole.G0))+
  scale_y_discrete(breaks=letters[1:13],
                   labels=colnames(GIC.sole.G0))+theme_bw()+
  xlab("GIC")+ylab("LIC")+ggtitle("Sole")+ labs(fill="Index value")+
  theme(axis.text.x = element_text(angle = 90)) 


#Plot annexe 1
cor_row <- cowplot::plot_grid(cordab,corplaice,corsole+ theme(legend.position="none"),nrow=1,ncol=3)
#Retrieve legend
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  corsole + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot::plot_grid(cor_row,legend,rel_widths = c(3, .4))

#ggsave("Annexe_2.pdf",path="Results/Figures")
#ggsave("Annexe_2.tiff",path="Results/Figures")





gic.values <- c(GIC.plaice.G0[lower.tri(GIC.plaice.G0, diag = FALSE)],
                GIC.sole.G0[lower.tri(GIC.sole.G0, diag = FALSE)],
                GIC.dab.G0[lower.tri(GIC.dab.G0, diag = FALSE)],
                GIC.flounder.G0[lower.tri(GIC.flounder.G0, diag = FALSE)])

lic.values <- c(LIC.plaice.G0[lower.tri(LIC.plaice.G0, diag = FALSE)],
                LIC.sole.G0[lower.tri(LIC.sole.G0, diag = FALSE)],
                LIC.dab.G0[lower.tri(LIC.dab.G0, diag = FALSE)],
                LIC.flounder.G0[lower.tri(LIC.flounder.G0, diag = FALSE)])

sp <- c(rep("Plaice",length(GIC.plaice.G0[lower.tri(GIC.plaice.G0, diag = FALSE)])),
        rep("Sole",length(GIC.sole.G0[lower.tri(GIC.sole.G0, diag = FALSE)])),
        rep("Dab",length(GIC.dab.G0[lower.tri(GIC.dab.G0, diag = FALSE)])),
        rep("Flounder",length(GIC.flounder.G0[lower.tri(GIC.flounder.G0, diag = FALSE)])))


booxingplot <- data.frame(GIC=gic.values,LIC=lic.values,Species=sp)
booxingplot <- booxingplot %>% pivot_longer(c(1,2), names_to="Index",values_to = "Value")

dodge <- position_dodge(width = 1)
ggplot(booxingplot, aes(x=Index,y=Value))+
  geom_violin(aes(fill=Species),position=dodge, alpha=0.7)+
  geom_boxplot(aes(Group=Species),width=.1, position=dodge)+
  scale_fill_viridis_d() + theme_minimal()

#Violin boxplot isotropy.
Iso <- ggplot(data=CG_dfG0) +
  geom_violin(aes(x=Species, y=Iso, fill=Species),alpha=0.7) + 
  geom_boxplot(aes(x=Species, y=Iso),width=0.1)+ 
  scale_fill_viridis_d() + theme_minimal() + ylab("Isotropy")

#####

#Networks##### 

#DAB
LIC.cor_df <- as.data.frame(LIC.dab.G0)
LIC.cor_df[lower.tri(LIC.cor_df)] <- NA
#Supprime les variables seules au réseau final
rownames(LIC.cor_df) <- colnames(LIC.cor_df)

LIC.cor_df <- tibble::add_column(LIC.cor_df, Var_x = rownames(LIC.cor_df), .before="1996")
LIC.cor_df <- LIC.cor_df %>% gather("Var_y","Value",2:14) 
LIC.cor_df[which(LIC.cor_df$Var_x==LIC.cor_df$Var_y),]$Value<-NA
head(LIC.cor_df)

LIC.cor_df <- na.omit(LIC.cor_df)
#LIC.cor_df <- LIC.cor_df %>%
  #dplyr::filter(LIC.cor_df$Value > quantile(LIC.cor_df$Value)[3])
head(LIC.cor_df)


GIC.cor_df <- as.data.frame(GIC.dab.G0)
GIC.cor_df[lower.tri(GIC.cor_df)] <- NA
#Supprime les variables seules au GIC.cor_df final
rownames(GIC.cor_df) <- colnames(GIC.cor_df)

GIC.cor_df <- tibble::add_column(GIC.cor_df, Var_x = rownames(GIC.cor_df), .before="1996")
GIC.cor_df <- GIC.cor_df %>% gather("Var_y","Value",2:14) 
GIC.cor_df[which(GIC.cor_df$Var_x==GIC.cor_df$Var_y),]$Value<- NA

GIC.cor_df <- na.omit(GIC.cor_df)
#GIC.cor_df <- GIC.cor_df %>%
  #dplyr::filter(GIC.cor_df$Value > quantile(GIC.cor_df$Value)[3])
head(GIC.cor_df)

GIC.cor_df$Index <- "GIC"
LIC.cor_df$Index <- "LIC"
vertex.dab <- rbind(GIC.cor_df,LIC.cor_df)

dab.graph <- vertex.dab %>%
  as_tbl_graph(dab.graph, directed = FALSE)

net_dab_G0 <- dab.graph %>%
  activate(nodes) %>%
  mutate(community = as.factor(group_louvain(weights = Value))) 
quantGIC <- quantile(GIC.cor_df$Value)[3]
quantLIC <- quantile(LIC.cor_df$Value)[3]

net_dab_G0 <-net_dab_G0 %>% activate(edges) %>% filter(Index == "GIC" & Value > quantGIC | Index== "LIC" & Value> quantLIC )
  graph_dab_G0 <- ggraph(net_dab_G0,layout = "lgl") + 
  geom_edge_link(aes(edge_alpha = Value, edge_width = Value, color = Value)) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(0, 1), colors = c("firebrick2", "dodgerblue2")) +
  #scale_edge_colour_gradient2(low = "blue", mid = "#FFFFCC", high = "red", midpoint = 0.5, breaks = c(0,0.25,0.5,0.75,1)) +
  scale_edge_colour_gradientn(limits = c(0, 1), colors =viridis(3))+
  geom_node_point(size = 4,pch=21,col="black", aes(fill=community)) +
  geom_node_text(aes(label = name), repel = TRUE)+
  theme_graph(foreground = 'white', fg_text_colour = 'black', base_family = 'Helvetica')+ 
  facet_edges(.~Index) + scale_fill_viridis_d() + ggtitle("Dab")
graph_dab_G0



#PLAICE

LIC.cor_df <- as.data.frame(LIC.plaice.G0)
LIC.cor_df[lower.tri(LIC.cor_df)] <- NA
#Supprime les variables seules au réseau final
rownames(LIC.cor_df) <- colnames(LIC.cor_df)

LIC.cor_df <- tibble::add_column(LIC.cor_df, Var_x = rownames(LIC.cor_df), .before="1996")
LIC.cor_df <- LIC.cor_df %>% gather("Var_y","Value",2:14) 
LIC.cor_df[which(LIC.cor_df$Var_x==LIC.cor_df$Var_y),]$Value<-NA
head(LIC.cor_df)

LIC.cor_df <- na.omit(LIC.cor_df)
#LIC.cor_df <- LIC.cor_df %>%
 # dplyr::filter(LIC.cor_df$Value > quantile(LIC.cor_df$Value)[3])
head(LIC.cor_df)


GIC.cor_df <- as.data.frame(GIC.plaice.G0)
GIC.cor_df[lower.tri(GIC.cor_df)] <- NA
#Supprime les variables seules au GIC.cor_df final
rownames(GIC.cor_df) <- colnames(GIC.cor_df)

GIC.cor_df <- tibble::add_column(GIC.cor_df, Var_x = rownames(GIC.cor_df), .before="1996")
GIC.cor_df <- GIC.cor_df %>% gather("Var_y","Value",2:14) 
GIC.cor_df[which(GIC.cor_df$Var_x==GIC.cor_df$Var_y),]$Value<- NA

GIC.cor_df <- na.omit(GIC.cor_df)
#GIC.cor_df <- GIC.cor_df %>%
 # dplyr::filter(GIC.cor_df$Value > quantile(GIC.cor_df$Value)[3])
head(GIC.cor_df)

GIC.cor_df$Index <- "GIC"
LIC.cor_df$Index <- "LIC"
vertex.plaice <- rbind(GIC.cor_df,LIC.cor_df)

plaice.graph <- vertex.plaice %>%
  as_tbl_graph(plaice.graph, directed = FALSE)

net_plaice_G0 <- plaice.graph %>%
  activate(nodes) %>%
  mutate(community = as.factor(group_louvain(weights = Value))) 
quantGIC <- quantile(GIC.cor_df$Value)[3]
quantLIC <- quantile(LIC.cor_df$Value)[3]

net_plaice_G0 <-net_plaice_G0 %>% activate(edges) %>% filter(Index == "GIC" & Value > quantGIC | Index== "LIC" & Value> quantLIC )
graph_plaice_G0 <- ggraph(net_plaice_G0,layout = "lgl") + 
  geom_edge_link(aes(edge_alpha = Value, edge_width = Value, color = Value)) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(0, 1), colors = c("firebrick2", "dodgerblue2")) +
  #scale_edge_colour_gradient2(low = "blue", mid = "#FFFFCC", high = "red", midpoint = 0.5, breaks = c(0,0.25,0.5,0.75,1)) +
  scale_edge_colour_gradientn(limits = c(0, 1), colors =viridis(3))+
  geom_node_point(size = 4,pch=21,col="black", aes(fill=community)) +
  geom_node_text(aes(label = name), repel = TRUE)+
  theme_graph(foreground = 'white', fg_text_colour = 'black', base_family = 'Helvetica')+ 
  facet_edges(.~Index) + scale_fill_viridis_d() + ggtitle("Plaice")
graph_plaice_G0



#SOLE

LIC.cor_df <- as.data.frame(LIC.sole.G0)
LIC.cor_df[lower.tri(LIC.cor_df)] <- NA
#Supprime les variables seules au réseau final
rownames(LIC.cor_df) <- colnames(LIC.cor_df)

LIC.cor_df <- tibble::add_column(LIC.cor_df, Var_x = rownames(LIC.cor_df), .before="1996")
LIC.cor_df <- LIC.cor_df %>% gather("Var_y","Value",2:14) 
LIC.cor_df[which(LIC.cor_df$Var_x==LIC.cor_df$Var_y),]$Value<-NA
head(LIC.cor_df)

LIC.cor_df <- na.omit(LIC.cor_df)
#LIC.cor_df <- LIC.cor_df %>%
 # dplyr::filter(LIC.cor_df$Value > quantile(LIC.cor_df$Value)[3])
head(LIC.cor_df)


GIC.cor_df <- as.data.frame(GIC.sole.G0)
GIC.cor_df[lower.tri(GIC.cor_df)] <- NA
#Supprime les variables seules au GIC.cor_df final
rownames(GIC.cor_df) <- colnames(GIC.cor_df)

GIC.cor_df <- tibble::add_column(GIC.cor_df, Var_x = rownames(GIC.cor_df), .before="1996")
GIC.cor_df <- GIC.cor_df %>% gather("Var_y","Value",2:14) 
GIC.cor_df[which(GIC.cor_df$Var_x==GIC.cor_df$Var_y),]$Value<- NA

GIC.cor_df <- na.omit(GIC.cor_df)
#GIC.cor_df <- GIC.cor_df %>%
 # dplyr::filter(GIC.cor_df$Value > quantile(GIC.cor_df$Value)[3])
head(GIC.cor_df)

GIC.cor_df$Index <- "GIC"
LIC.cor_df$Index <- "LIC"
vertex.sole <- rbind(GIC.cor_df,LIC.cor_df)

sole.graph <- vertex.sole %>%
  as_tbl_graph(sole.graph, directed = FALSE)

net_sole_G0 <- sole.graph %>%
  activate(nodes) %>%
  mutate(community = as.factor(group_louvain(weights = Value))) 
quantGIC <- quantile(GIC.cor_df$Value)[3]
quantLIC <- quantile(LIC.cor_df$Value)[3]

net_sole_G0 <-net_sole_G0 %>% activate(edges) %>% filter(Index == "GIC" & Value > quantGIC | Index== "LIC" & Value> quantLIC )
graph_sole_G0 <- ggraph(net_sole_G0,layout = "lgl") + 
  geom_edge_link(aes(edge_alpha = Value, edge_width = Value, color = Value)) +
  guides(edge_alpha = "none", edge_width = "none") +
  #scale_edge_colour_gradientn(limits = c(0, 1), colors = c("firebrick2", "dodgerblue2")) +
  #scale_edge_colour_gradient2(low = "blue", mid = "#FFFFCC", high = "red", midpoint = 0.5, breaks = c(0,0.25,0.5,0.75,1)) +
  scale_edge_colour_gradientn(limits = c(0, 1), colors =viridis(3))+
  geom_node_point(size = 4,pch=21,col="black", aes(fill=community)) +
  geom_node_text(aes(label = name), repel = TRUE)+
  theme_graph(foreground = 'white', fg_text_colour = 'black', base_family = 'Helvetica')+ 
  facet_edges(.~Index) + scale_fill_viridis_d() + ggtitle("Sole")
graph_sole_G0

networks <-cowplot::plot_grid(graph_dab_G0,graph_plaice_G0,graph_sole_G0,nrow=3)
networks

#ggsave("Figure_3.pdf",path="Results/Figures", height = 28, width=20, units = "cm")
#ggsave("Figure_3.tiff",path="Results/Figures", height = 28, width=20, units = "cm")

#Alluvial
dab.com <- net_dab_G0$data[,c(3,4)]
plaice.com <- net_plaice_G0$data[,c(3,4)]
sole.com <- net_sole_G0$data[,c(3,4)]

tmp <- dplyr::left_join(dab.com,plaice.com, by="name")

com <- left_join(tmp,sole.com, by="name")

dab.com$sp <- "Dab"
plaice.com$sp <- "Plaice"
sole.com$sp <- "Sole"
com <- rbind(dab.com,plaice.com,sole.com)

gg <- ggplot(com,
             aes(x = sp, stratum = community, alluvium = name,
                 fill = community)) +
  geom_stratum()
gg + geom_flow(stat = "alluvium", color = "black") +
  geom_text(aes(label = name), stat = "alluvium") + scale_fill_viridis_d()+ theme_bw()

#####

#MAF construction#####
#Initialize loop
vp <- numeric() #eigenvalue
maf1 <- numeric() #maf value
maf2 <- numeric()
maf3 <- numeric()
maf4 <- numeric()
Lon <- numeric() #longitude
Lat <- numeric() #latitude   
taxa.maf <- character() #taxa selected
group.maf <- character() #group selected
period <- character()


sp <- unique(Pichon_dens_corres$Nv_Rubbin)
gp <- unique(Pichon_dens_corres$Groupe)

#Data formatting, spread species x group variable 
data.maf <- Pichon_dens_corres
data.maf$cat <- paste(data.maf$Nv_Rubbin,data.maf$Groupe,data.maf$annee, sep="_")
data.maf <- data.maf %>% distinct() %>% tidyr::spread(cat, Densite) %>%  dplyr::select(-Numéro_du_Trait, -Code_Station, -lonmoy, -latmoy, -Groupe, -Nv_Rubbin, -Correspondance, -annee) 
data.maf <- as.data.frame(data.maf)
data.maf[which(is.na(data.maf),arr.ind = T)] <- 0
data.maf <- data.maf %>% group_by(loncor, latcor) %>% summarise_each(funs(sum))

#Loop all

for(x in sp){
  for(i in gp){
    z <- which(grepl(paste0(x,"_",i),names(data.maf)))
    maf.res <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)
    vp <- c(vp,maf.res$eig)
    MAF_factors <- maf.res$x
    maf1 <- c(maf1,MAF_factors[,1])
    maf2 <- c(maf2,MAF_factors[,2])
    maf3 <- c(maf3,MAF_factors[,3])
    maf4 <- c(maf4,MAF_factors[,4])
    taxa.maf <- c(taxa.maf,rep(x,length(MAF_factors[,1])))
    group.maf <- c(group.maf,rep(i,length(MAF_factors[,1])))
    period <- c(period,rep("1996-2019",length(MAF_factors[,1])))
    Lon <- c(Lon,data.maf$loncor)
    Lat <- c(Lat,data.maf$latcor)
  }
}

#Data frames 
maf.df <- data.frame(Taxa=taxa.maf,
                     Group=group.maf,
                     Period=period,
                     Lon=Lon,
                     Lat=Lat,
                     MAF1=maf1,
                     MAF2=maf2,
                     MAF3=maf3,
                     MAF4=maf4)

# ev.df <- data.frame(Taxa=c(rep(sp,each=8),rep(sp,each=8)),
#                     Group=rep(gp,each=4,times=8),
#                     MAF=c(rep(c("MAF1","MAF2","MAF3","MAF4"), times=8),rep(c("MAF1","MAF2","MAF3","MAF4"), times=8)),
#                     Period=rep(c("1996-2002","2008-2019"),each=32),
#                     EV=vp)

maf.df <- maf.df %>% pivot_longer(c(6:9), names_to ="MAF", values_to ="Value")
#maf.df <- left_join(maf.df,ev.df,"Taxa"="Taxa","Group"="Group","MAF"="MAF","Period"="Period")

#Graphs for G0
maf.df.dab <- maf.df %>% filter(Taxa %in% c("LIMALIM") & Group %in% c("G0") )
maf.df.plaice <- maf.df %>% filter(Taxa %in% c("PLEUPLA") & Group %in% c("G0") )
maf.df.sole <- maf.df %>% filter(Taxa %in% c("SOLEVUL") & Group %in% c("G0") )
maf.df.flounder <- maf.df %>% filter(Taxa %in% c("PLATFLE") & Group %in% c("G0") )

coastMAF <- raster::crop(coast, raster::extent(-0.2,0.25,49.25,49.7))
### G0
toto <- maf.df.dab %>% filter(MAF %in% c("MAF1","MAF2","MAF3"))
toto$Color <- ifelse(toto$Value>0,1,2)
dab.maf <-ggplot(toto)+
  geom_point(aes(x = Lon, y = Lat, col=Color, size=abs(Value)), pch=19, alpha=0.5) +  
  theme_bw()+ 
  xlab("Longitude (dec.)")+ylab("Latitude (dec.)") +
  geom_polygon(data = coastMAF, aes(x = long, y = lat, group = group),fill="grey",col="black")+
  facet_wrap(.~MAF) +ggtitle("Dab") + guides(col=F) + coord_fixed() + labs(size="absolute\nMAF value") +
  scale_color_viridis(option="D")


toto <- maf.df.plaice %>% filter(MAF %in% c("MAF1","MAF2","MAF3"))
toto$Color <- ifelse(toto$Value>0,1,2)
plaice.maf <- ggplot(toto)+
  geom_point(aes(x = Lon, y = Lat, col=Color, size=abs(Value)), pch=19, alpha=0.5) +  
  theme_bw() + 
  xlab("Longitude (dec.)")+ylab("Latitude (dec.)") +
  geom_polygon(data = coastMAF, aes(x = long, y = lat, group = group),fill="grey",col="black")+
  facet_wrap(.~MAF) +ggtitle("Plaice") + guides(col=F)+ coord_fixed() + labs(size="absolute\nMAF value") +
  scale_color_viridis(option="D")

toto <- maf.df.sole %>% filter(MAF %in% c("MAF1","MAF2","MAF3"))
toto$Color <- ifelse(toto$Value>0,1,2)
sole.maf <-ggplot(toto)+
  geom_point(aes(x = Lon, y = Lat, col=Color, size=abs(Value)), pch=19, alpha=0.5) +  
  theme_bw() + 
  xlab("Longitude (dec.)")+ylab("Latitude (dec.)") +
  geom_polygon(data = coastMAF, aes(x = long, y = lat, group = group),fill="grey",col="black")+
  facet_wrap(.~MAF) +ggtitle("Sole") + guides(col=F)+ coord_fixed() + labs(size="absolute\nMAF value") +
  scale_color_viridis(option="D")

# toto <- maf.df.flounder %>% filter(MAF %in% c("MAF1","MAF2","MAF3"))
# toto$Color <- ifelse(toto$Value>0,1,2)
# flounder.maf <- ggplot(toto)+
#   geom_point(aes(x = Lon, y = Lat, col=Color, size=abs(Value)), pch=19, alpha=0.5) +  
#   theme_bw() + 
#   xlab("Longitude (dec.)")+ylab("Latitude (dec.)") +
#   geom_path(data = cotes, aes(x = long, y = lat, group = group)) +
#   facet_wrap(.~MAF) +ggtitle("Flounder") + guides(col=F)+ coord_fixed()+
#   scale_color_viridis(option="D")

cowplot::plot_grid(dab.maf,plaice.maf, sole.maf,
                   nrow=3)

#ggsave("Annexe_3.pdf",path="Results/Figures")
#ggsave("Annexe_3.tiff",path="Results/Figures")

#####

#MAF variogram#####

z <- which(grepl(paste0("LIMALIM_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)


vario.db <- res.all$x
vario.db <- cbind(as.matrix(data.maf[,c(1,2)]),vario.db)

vario.db <-db.create(vario.db)
vario.db <-db.locate(vario.db,c("loncor","latcor"),"x")
vario.db <-db.locate(vario.db,c(4),"z")
varMAF1 <- RGeostats::vario.calc(vario.db)
plot(varMAF1)

vario.df <- data.frame(
  Distance=c(varMAF1@vardirs[[1]]$hh),
  Variance=c(varMAF1@vardirs[[1]]$gg),
  Size=c(varMAF1@vardirs[[1]]$sw),
  Species=rep("Dab",length(varMAF1@vardirs[[1]]$gg))
)

z <- which(grepl(paste0("PLEUPLA_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)


vario.db <- res.all$x
vario.db <- cbind(as.matrix(data.maf[,c(1,2)]),vario.db)

vario.db <-db.create(vario.db)
vario.db <-db.locate(vario.db,c("loncor","latcor"),"x")
vario.db <-db.locate(vario.db,c(4),"z")
varMAF1 <- RGeostats::vario.calc(vario.db)
plot(varMAF1)

tmp <- data.frame(
  Distance=c(varMAF1@vardirs[[1]]$hh),
  Variance=c(varMAF1@vardirs[[1]]$gg),
  Size=c(varMAF1@vardirs[[1]]$sw),
  Species=rep("Plaice",length(varMAF1@vardirs[[1]]$gg))
)
vario.df <- rbind(vario.df,tmp)

z <- which(grepl(paste0("SOLEVUL_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)


vario.db <- res.all$x
vario.db <- cbind(as.matrix(data.maf[,c(1,2)]),vario.db)

vario.db <-db.create(vario.db)
vario.db <-db.locate(vario.db,c("loncor","latcor"),"x")
vario.db <-db.locate(vario.db,c(4),"z")
varMAF1 <- RGeostats::vario.calc(vario.db)
plot(varMAF1)

tmp <- data.frame(
  Distance=c(varMAF1@vardirs[[1]]$hh),
  Variance=c(varMAF1@vardirs[[1]]$gg),
  Size=c(varMAF1@vardirs[[1]]$sw),
  Species=rep("Sole",length(varMAF1@vardirs[[1]]$gg))
)
vario.df <- rbind(vario.df,tmp)


maf_vario <- ggplot(vario.df)+
  geom_point(aes(x=Distance,y=Variance,size=Size,col=Species))+
  geom_line(aes(x=Distance,y=Variance,col=Species))+
  scale_color_viridis_d()+theme_bw()

maf_vario
#ggsave("Figure_4.pdf", path="Results/Figures")
#ggsave("Figure_4.tiff", path="Results/Figures")
#####

#MAF clustering#####
#https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}
set_labels_params <- function(nbLabels,direction = c("tb", "bt", "lr", "rl"),fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}
plot_ggdendro <- function(hcdata,direction   = c("lr", "rl", "tb", "bt"),
                          fan = FALSE,scale.color = NULL,branch.size = 1,label.size  = 3,nudge.label = 0.01,expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}



#Cluster LIMANDE
z <- which(grepl(paste0("LIMALIM_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)

library(dendextend)
toto <- data.frame(res.all$c1)
rownames(toto) <- c(substr(rownames(toto),12,15))

nb_maf <- 3

dist <- matrix(0,ncol=dim(toto)[1], nrow=dim(toto)[1])
for(i in 1:dim(dist)[1]){
  for(j in 1:dim(dist)[1]){
    dist[i,j] <- sum((toto[i,1:nb_maf]-toto[j,1:nb_maf])^2)
  }
}
rownames(dist) <- rownames(toto)
colnames(dist) <- rownames(toto)

arbre <- hclust(as.dist(dist), method="ward.D2")
ddata <- dendro_data(as.dendrogram(arbre), type = "rectangle")
dab.dend <- ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + ggtitle("Dab dendrogram using Ward")
dab.dend

clust <-cascadeKM(toto[,c(1:2)],2,6,iter=999, criterion = "ssi")
plot(clust)
part.dab <- clust$partition


dendro_cut <- dendro_data_k(arbre, 5)

dab_dendro <- plot_ggdendro(dendro_cut,
                            direction   = "lr",
                            expand.y    = 0.2, 
                            scale.color=viridis(6))+theme_minimal()+ggtitle("Dab")




toto <- data.maf[,c(1,2,z)]
colnames(toto)[3:15] <- c(substr(colnames(toto)[3:15],12,15))
toto <- toto %>% pivot_longer(3:15,names_to="Annee",values_to="Densite")

toto<- toto %>% left_join(data.frame(Annee=rownames(clust$partition),Group=clust$partition[,4]),"Annee"="Annee")
for(i in 1:max(toto$Group)){
  toto$Group[toto$Group==i] <- paste0(as.character(unique(toto$Annee[toto$Group==i])),collapse="-")
  
}
toto <- toto %>% group_by(Group,loncor,latcor) %>% summarise(Mean=mean(Densite))

dab.group <-ggplot()+
  geom_point(data=toto[toto$Mean !=0,], aes(x=loncor,y=latcor, size=abs(Mean), fill=Mean),alpha=0.7, pch=21)+
  geom_point(data=toto[toto$Mean ==0,], aes(x=loncor,y=latcor, size=abs(Mean)), pch=4)+
  scale_fill_viridis()+ 
  geom_path(data = cotes, aes(x = long, y = lat, group = group))+
  facet_wrap(~Group) + theme_minimal() + xlab("Longitude(°)") + ylab("Latitude(°)")+ coord_fixed()



#Cluster FLET
# z <- which(grepl(paste0("PLATFLE_G0"),names(data.maf)))
# res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)
# 
# library(dendextend)
# toto <- data.frame(res.all$c1)
# rownames(toto) <- c(substr(rownames(toto),12,15))
# 
# nb_maf <- 2
# dist <- matrix(0,ncol=dim(toto)[1], nrow=dim(toto)[1])
# for(i in 1:dim(dist)[1]){
#   for(j in 1:dim(dist)[1]){
#     dist[i,j] <- sum((toto[i,1:nb_maf]-toto[j,1:nb_maf])^2)
#   }
# }
# rownames(dist) <- rownames(toto)
# colnames(dist) <- rownames(toto)
# 
# arbre <- hclust(as.dist(dist), method="ward.D2")
# ddata <- dendro_data(as.dendrogram(arbre), type = "rectangle")
# flounder.dend <- ggplot() + 
#   geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
#   geom_text(data=label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         panel.background=element_rect(fill="white"),
#         panel.grid=element_blank()) + ggtitle("Flounder dendrogram using Ward")
# flounder.dend
# clust <-cascadeKM(toto[,c(1:2)],2,6,iter=999, criterion = "ssi")
# plot(clust)
# clust$partition
# 
# 
# toto <- data.maf[,c(1,2,z)]
# colnames(toto)[3:15] <- c(substr(colnames(toto)[3:15],12,15))
# toto <- toto %>% pivot_longer(3:15,names_to="Annee",values_to="Densite")
# 
# toto<- toto %>% left_join(data.frame(Annee=rownames(clust$partition),Group=clust$partition[,2]),"Annee"="Annee")
# for(i in 1:max(toto$Group)){
#   toto$Group[toto$Group==i] <- paste0(as.character(unique(toto$Annee[toto$Group==i])),collapse="-")
#   
# }
# toto <- toto %>% group_by(Group,loncor,latcor) %>% summarise(Mean=mean(Densite))
# 
# flounder.group <-ggplot()+
#   geom_point(data=toto[toto$Mean !=0,], aes(x=loncor,y=latcor, size=abs(Mean), fill=Mean),alpha=0.7, pch=21)+
#   geom_point(data=toto[toto$Mean ==0,], aes(x=loncor,y=latcor, size=abs(Mean)), pch=4)+
#   scale_fill_viridis()+ 
#   geom_path(data = cotes, aes(x = long, y = lat, group = group))+
#   facet_wrap(~Group, ncol=2) + theme_minimal() + xlab("Longitude(°)") + ylab("Latitude(°)")+ coord_fixed()


#Cluster PLAICE
z <- which(grepl(paste0("PLEUPLA_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)

library(dendextend)
toto <- data.frame(res.all$c1)
rownames(toto) <- c(substr(rownames(toto),12,15))

nb_maf <- 3

dist <- matrix(0,ncol=dim(toto)[1], nrow=dim(toto)[1])
for(i in 1:dim(dist)[1]){
  for(j in 1:dim(dist)[1]){
    dist[i,j] <- sum((toto[i,1:nb_maf]-toto[j,1:nb_maf])^2)
  }
}
rownames(dist) <- rownames(toto)
colnames(dist) <- rownames(toto)

arbre <- hclust(as.dist(dist), method="ward.D2")
ddata <- dendro_data(as.dendrogram(arbre), type = "rectangle")
plaice.dend <- ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + ggtitle("Plaice dendrogram using Ward")
plaice.dend


clust <-cascadeKM(toto[,c(1:2)],2,6,iter=100, criterion = "ssi")
plot(clust)
part.plaice <- clust$partition


dendro_cut <- dendro_data_k(arbre, 5)

plaice_dendro <- plot_ggdendro(dendro_cut,
                               direction   = "lr",
                               expand.y    = 0.2, 
                               scale.color=viridis(6))+theme_minimal()+ggtitle("Plaice")



toto <- data.maf[,c(1,2,z)]
colnames(toto)[3:15] <- c(substr(colnames(toto)[3:15],12,15))
toto <- toto %>% pivot_longer(3:15,names_to="Annee",values_to="Densite")

toto<- toto %>% left_join(data.frame(Annee=rownames(clust$partition),Group=clust$partition[,5]),"Annee"="Annee")
for(i in 1:max(toto$Group)){
  toto$Group[toto$Group==i] <- paste0(as.character(unique(toto$Annee[toto$Group==i])),collapse="-")
  
}
toto <- toto %>% group_by(Group,loncor,latcor) %>% summarise(Mean=mean(Densite))

plaice.group <-ggplot()+
  geom_point(data=toto[toto$Mean !=0,], aes(x=loncor,y=latcor, size=abs(Mean), fill=Mean),alpha=0.7, pch=21)+
  geom_point(data=toto[toto$Mean ==0,], aes(x=loncor,y=latcor, size=abs(Mean)), pch=4)+
  scale_fill_viridis()+ 
  geom_path(data = cotes, aes(x = long, y = lat, group = group))+
  facet_wrap(~Group) + theme_minimal() + xlab("Longitude(°)") + ylab("Latitude(°)")+ coord_fixed()







#Cluster SOLE
z <- which(grepl(paste0("SOLEVUL_G0"),names(data.maf)))
res.all <- maf.f(data.maf[,c(z)],as.matrix(data.maf[,c(1,2)]),hmin=0,hmax=0.05)

library(dendextend)
toto <- data.frame(res.all$c1)
rownames(toto) <- c(substr(rownames(toto),12,15))

nb_maf <- 3

dist <- matrix(0,ncol=dim(toto)[1], nrow=dim(toto)[1])
for(i in 1:dim(dist)[1]){
  for(j in 1:dim(dist)[1]){
    dist[i,j] <- sum((toto[i,1:nb_maf]-toto[j,1:nb_maf])^2)
  }
}
rownames(dist) <- rownames(toto)
colnames(dist) <- rownames(toto)

arbre <- hclust(as.dist(dist), method="ward.D2")
ddata <- dendro_data(as.dendrogram(arbre), type = "rectangle")
sole.dend <- ggplot() + 
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ddata), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + ggtitle("Sole dendrogram using Ward")
sole.dend

clust <-cascadeKM(toto[,c(1:2)],2,6,iter=999, criterion = "ssi")
plot(clust)
part.sole <-clust$partition


dendro_cut <- dendro_data_k(arbre, 5)

sole_dendro <- plot_ggdendro(dendro_cut,
                             direction   = "lr",
                             expand.y    = 0.2, 
                             scale.color=viridis(6))+theme_minimal()+ggtitle("Sole")




toto <- data.maf[,c(1,2,z)]
colnames(toto)[3:15] <- c(substr(colnames(toto)[3:15],12,15))
toto <- toto %>% pivot_longer(3:15,names_to="Annee",values_to="Densite")

toto<- toto %>% left_join(data.frame(Annee=rownames(clust$partition),Group=clust$partition[,4]),"Annee"="Annee")
for(i in 1:max(toto$Group)){
  toto$Group[toto$Group==i] <- paste0(as.character(unique(toto$Annee[toto$Group==i])),collapse="-")
  
}
toto <- toto %>% group_by(Group,loncor,latcor) %>% summarise(Mean=mean(Densite))

sole.group <-ggplot()+
  geom_point(data=toto[toto$Mean !=0,], aes(x=loncor,y=latcor, size=abs(Mean), fill=Mean),alpha=0.7, pch=21)+
  geom_point(data=toto[toto$Mean ==0,], aes(x=loncor,y=latcor, size=abs(Mean)), pch=4)+
  scale_fill_viridis()+ 
  geom_path(data = cotes, aes(x = long, y = lat, group = group))+
  facet_wrap(~Group, ncol=3) + theme_minimal() + xlab("Longitude(°)") + ylab("Latitude(°)") + coord_fixed()




#####

#Combined plot#####
cowplot::plot_grid(dab.dend,plaice.dend,sole.dend, nrow=2,ncol=2)
cowplot::plot_grid(dab.group,plaice.group,sole.group, nrow=3,ncol=1)
cowplot::plot_grid(dab_dendro,plaice_dendro,sole_dendro, nrow=3,ncol=1)
#ggsave("Figure_5.pdf", path="Results/Figures")
#ggsave("Figure_5.tiff", path="Results/Figures")

dab.com.maf <- data.frame(name=rownames(part.dab),community=part.dab[,4])
plaice.com.maf <- data.frame(name=rownames(part.plaice),community=part.plaice[,4])
sole.com.maf <- data.frame(name=rownames(part.sole),community=part.sole[,4])

dab.com.maf$sp <- "Dab"
plaice.com.maf$sp <- "Plaice"
sole.com.maf$sp <- "Sole"
com.maf <- rbind(dab.com.maf,plaice.com.maf,sole.com.maf)
com.maf$community <- factor(com.maf$community)

gg <- ggplot(com.maf,
             aes(x = sp, stratum = community, alluvium = name,
                 fill = community)) +
  geom_stratum()
gg + geom_flow(stat = "alluvium", color = "black") +
  geom_text(aes(label = name), stat = "alluvium") + scale_fill_viridis_d()+ theme_bw()
#####

#Clustering comparing#####
com.maf$cluster <- "maf"
com$cluster <- "index"
com.comp <- rbind(com.maf,com)
com.comp$community <- factor(com.comp$community)

gg <- ggplot(com.comp,
             aes(x = cluster, stratum = community, alluvium = name,
                 fill = community)) +
  geom_stratum()+ geom_flow(stat = "alluvium", color = "black") +
  geom_text(aes(label = name), stat = "alluvium") + scale_fill_viridis_d()+ theme_bw()+
  facet_wrap(.~sp)
gg
#####

# MnM figs #####

load("Data/cotesseine.RData")
load("Data/stratepoly.RData")


library(sf)
library(ggplot2)
library(cowplot)
library(spData)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

library(sp)
library(dplyr)
library(tidyr)
library(ggmap)
library(viridis)

##V1####
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

pt <- Pichon_dens_corres %>% dplyr::select(lonmoy,latmoy) %>% unique()
pt <- st_as_sf(pt, coords = c("lonmoy", "latmoy"), crs= 4326, agr="constant")

seine_bbox <- st_as_sfc(st_bbox(pt))



ggm1 <- ggplot() + 
  geom_sf(data = world, fill = "grey") + 
  geom_sf(data = seine_bbox, fill = NA, color = "red", size = 1.2)+
  coord_sf(xlim = c(-5, 5), ylim = c(42, 52), expand = FALSE)+
  theme_void() + xlab("Longitude") +ylab("Latitude")


ggm1




poly<- rgdal::readOGR("Data/secteur_poly.shp")
poly <- spTransform(poly, CRS("+proj=longlat +datum=WGS84"))

pt <- Pichon_dens_corres %>% dplyr::select(loncor,latcor) %>% unique()
bc_bbox <- make_bbox(lat = latcor, lon = loncor, data = pt)

#carte <-get_map(location =bc_bbox, maptype="toner-lite", color="bw")
#ggmap(carte)

poly <- poly[poly$Ifremer_id != "B",]
centroids <- coordinates(poly)
colnames(centroids) <- c("long","lat")
centroids <- as.data.frame(centroids)
centroids$id <- as.character(unique(stratepoly.df$Ifremer_id))

rivers <- data.frame(long=c(-0.25,0.3),lat=c(49.2,49.4),legend=c("Orne","Seine"))
city <- data.frame(long=c(-0.32,0.12,0.2),lat=c(49.28,49.52,49.68),legend=c("Ouistreham","Le Havre","Antifer"))

coast <- rgdal::readOGR("Data/Coast_estuary_detailled.gpkg")
coast <- raster::crop(coast, raster::extent(-0.4,0.4,49.15,49.8))

ggm2 <- ggplot()+geom_point(data=pt, aes(x=loncor,y=latcor),alpha=0.7,col="black",size=2.4)+
  geom_polygon(data = coast, aes(x = long, y = lat, group = group),fill="grey", col="black")+
  geom_polygon(data=stratepoly.df, aes(x=long,y=lat,group=group),col="black",fill=NA)+
  #geom_text(data=centroids,aes(x=long,y=lat,label=id, fontface=2))+
  geom_text(data=rivers,aes(x=long,y=lat,label=legend, fontface=2))+
  geom_text(data=city,aes(x=long,y=lat,label=legend, fontface=3))+
  theme_minimal() + xlab("Longitude (dec.)") + ylab("Latitude (dec.)") +
  ggtitle("") + coord_fixed()+
  theme(legend.position="bottom") +
  ggsn::scalebar(dist=10, transform=T, dist_unit ="km",location="bottomright",
                 x.min=(-0.4), x.max=0.35, y.min=49.2, y.max=49.8, st.size=3)
ggm2

gg_inset_map <- ggdraw() +
  draw_plot(ggm2) +
  draw_plot(ggm1, x = 0.06 ,y = 0.585, width = 0.30, height = 0.35)

gg_inset_map
#ggsave("Figure_0.pdf", path="Results/Figures")
#ggsave("Figure_0.tiff", path="Results/Figures")
####

##V2####

world <- ne_countries(scale = "large", returnclass = "sf")
#Element to add
rivers <- data.frame(long=c(-0.25,0.3),lat=c(49.2,49.4),legend=c("Orne","Seine"))
city <- data.frame(long=c(-0.32,0.12,0.2),lat=c(49.28,49.52,49.68),legend=c("Ouistreham","Le Havre","Antifer"))
pt <- Pichon_dens_corres %>% dplyr::select(loncor,latcor) %>% unique()

#Bathymetry
bathy <- rgdal::readOGR("C:/Users/thibcari/Desktop/Brouillon/Angie/Angie costatis/bathy.gpkg")
bathy25 <- bathy[bathy@data$DEPTH==c(-25),]
bathy25 <- fortify(bathy25)
bathy50 <- bathy[bathy@data$DEPTH==c(-50),]
bathy50 <- fortify(bathy50)
bathy100 <- bathy[bathy@data$DEPTH==c(-100),]
bathy100 <- fortify(bathy100)

bathy25 <- cbind(bathy25,Depth=25)
bathy50 <- cbind(bathy50,Depth=50)
bathy100 <- cbind(bathy100,Depth=100)
bathy2 <- rbind(bathy25,bathy50,bathy100)
bathy2$Depth <- factor(bathy2$Depth)

#Plot Channel
EastChan <- ggplot() + 
  geom_path(data = cotes, aes(x = long, y = lat, group = group))+
  #geom_polygon(data=stratepoly.df, aes(x=long,y=lat,group=group),col="black",fill=NA)+
  coord_sf(xlim = c(-0.4, 0.5), ylim = c(49.15, 49.7), expand = FALSE)+
  geom_text(data=rivers,aes(x=long,y=lat,label=legend, fontface=2))+
  geom_text(data=city,aes(x=long,y=lat,label=legend, fontface=3))+
  geom_point(data=pt,aes(x=loncor,y=latcor))+
  geom_path(data=bathy25,aes(x=long,y=lat, group=group),alpha=0.3)+
  theme_classic() +xlab("Longitude") +ylab("Latitude")

#Plot France or Europe  
library(sf)
box <- bc_bbox <- make_bbox(lat = latcor, lon = loncor, data = pt)
pt <- st_as_sf(box, coords = c("Coord_x", "Coord_y"), crs= 4326, agr="constant")
bbox <- st_as_sfc(st_bbox(pt))

FR <- ggplot() + 
  geom_sf(data = world, fill = "white") + 
  geom_sf(data = bbox, fill = NA, color = "red", size = 1.2)+
  coord_sf(xlim = c(-5, 5), ylim = c(42, 52), expand = FALSE)+
  theme_void()

WestEU <- ggplot() + 
  geom_sf(data = world, fill = "white") + 
  geom_sf(data = bbox, fill = NA, color = "red", size = 1.2)+
  coord_sf(xlim = c(-10, 10), ylim = c(37, 57), expand = FALSE)+
  theme_void()



#Combine plot
library(cowplot)
ggdraw() +
  draw_plot(EastChan) +
  draw_plot(FR, x = 0.70,y = 0.17, width = 0.30, height = 0.30)






















