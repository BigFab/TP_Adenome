
# Initialisations

library(rJava)
library(ggplot2)
library(corrplot)
library(cluster)
library(factoextra)
library(xlsx)
library(reshape2)
library(rpart)
library(rpart.plot)
theme_set(theme_classic())

# Plot Functions

adenome_scatterplot <- function(df, feature_1, feature_2, f1, f2) {
  x_axis = f1
  y_axis = f2
  title = "RTUPB/VBPPS"
  subtitle = paste(x_axis, " Vs. ", y_axis)
  g <- ggplot(df, aes(x=feature_1, y=feature_2))
  g <- g + geom_point()
  g <- g + labs(title=title,y=y_axis,x=x_axis)
  print(g)
}

# Obtenir le triangle inférieur d'une matrice
get_lower_tri <- function(mat){
    mat[upper.tri(mat)] <- NA
    return(mat)
}

# Obtenir le triangle supérieur d'une matrice
get_upper_tri <- function(mat){
    mat[lower.tri(mat)] <- NA
    return(mat)
}

# Re-ordonne la matrice de corrélation
reorder_cormat <- function(cormat){
    # Utiliser la corrélation entre les variables
    # comme mesure de distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
}

adenome_heatmap <- function(cormat) { 
    # Reordonner la matrice de corrélation
    cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)
    # Fondre la matrice de corrélation
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
    # Créer le heatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
        midpoint = 0, limit = c(-1,1), space = "Lab",
        name="Pearson\nCorrelation") +
        theme_minimal()+ # minimal theme
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
        size = 12, hjust = 1))+
        coord_fixed()
    
    ggheatmap <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))
    
    # Afficher heatmap
    print(ggheatmap)
}

# Import du dataset complet XLSX
rtupb_vbpps_complete <- read.xlsx(file = "datasets/RTUPB-VBPPS.xlsx", sheetIndex = 1, header = 1, startRow = 2)[,-1]
attributes(rtupb_vbpps_complete)$names = c('Age','Comorbidite','Duree_Traitement_Medical','Porteur_Sonde','IPSS','QoL','Qmax','PSA','Volume_Prostatique','Residu_Post_Mictionnel','Indication','Anesthesie','Evenement','Technique','Transfusion','Temps_Operation','Volume_Reseque','Delai_Ablation','Caillotage','Reprise_Bloc',
    'M1_IPSS','M1_QoL','M1_Qmax','M3_IPSS','M3_QoL','M3_Qmax','M6_IPSS','M6_QoL','M6_Qmax','M9_IPSS','M9_QoL','M9_Qmax','M12_IPSS','M12_QoL','M12_Qmax','M15_IPSS','M15_QoL','M15_Qmax','M18_IPSS','M18_QoL','M18_Qmax')

# Cast des variables catégoriques
rtupb_vbpps_complete$Comorbidite<-as.factor(rtupb_vbpps_complete$Comorbidite)
rtupb_vbpps_complete$Porteur_Sonde<-as.factor(rtupb_vbpps_complete$Porteur_Sonde)
rtupb_vbpps_complete$Caillotage<-as.factor(rtupb_vbpps_complete$Caillotage)
rtupb_vbpps_complete$Transfusion<-as.factor(rtupb_vbpps_complete$Transfusion)
rtupb_vbpps_complete$Reprise_Bloc<-as.factor(rtupb_vbpps_complete$Reprise_Bloc)
rtupb_vbpps_complete$Anesthesie<-as.factor(rtupb_vbpps_complete$Anesthesie)
rtupb_vbpps_complete$Indication<-as.factor(rtupb_vbpps_complete$Indication)
rtupb_vbpps_complete$Evenement<-as.factor(rtupb_vbpps_complete$Evenement)
rtupb_vbpps_complete$Technique<-as.factor(rtupb_vbpps_complete$Technique)

rtupb_vbpps_complete[c(17,30,36),]

rtupb_vbpps_complete[c(2,32),]

# Préparation du jeu de données

# on garde tous les individus
# on supprime les features que nous avons choisies d'ignorer
rtupb_vbpps <- subset(rtupb_vbpps_complete,select=-c(Volume_Reseque,Residu_Post_Mictionnel,Qmax,Evenement,Transfusion,Reprise_Bloc))

# 14 premières colonnes du dataset -> pré-opératoire, les dernières -> post-opératoire
rtupb_vbpps_pre <- rtupb_vbpps[,1:14]
rtupb_vbpps_post <- rtupb_vbpps[,c(11,15:35)]

# On crée des dataframes ne comportant que les variables numériques (pour corrélation linéaire de Pearson)
rtupb_vbpps_num <- rtupb_vbpps[,sapply(rtupb_vbpps, function(x) is.numeric(x))]
rtupb_vbpps_pre_num = rtupb_vbpps_pre[,sapply(rtupb_vbpps_pre, function(x) is.numeric(x))]
rtupb_vbpps_post_num = rtupb_vbpps_post[,sapply(rtupb_vbpps_post, function(x) is.numeric(x))]

# On crée des dataframes ne comportant que les données booléennes, ordinales ou catégoriques
rtupb_vbpps_cat <- rtupb_vbpps[,sapply(rtupb_vbpps, function(x) is.ordered(x) | is.factor(x) | is.logical(x))]
rtupb_vbpps_pre_cat = rtupb_vbpps_pre[,sapply(rtupb_vbpps_pre, function(x) is.ordered(x) | is.factor(x) | is.logical(x))]

t <- table(rtupb_vbpps_complete$Technique)
pie(x = t, labels = c("RTUPB","VBPPS"),main = "RTUPB/VBPPS - Pré-opératoire\nDistribution des patients selon la technique opératoire", radius = 0.5)

t <- table(rtupb_vbpps_complete$Porteur_Sonde)
pie(x = t, labels = c("Non porteur","Porteur"),main = "RTUPB/VBPPS - Pré-opératoire\nDistribution des patients par porteur de sonde", radius = 0.5)

# on divise le dataset en 2 pour voir les techniques séparément
rtupb_pre <- rtupb_vbpps_pre[which(rtupb_vbpps_pre$Technique == 1),]
vbpps_pre <- rtupb_vbpps_pre[which(rtupb_vbpps_pre$Technique == 2),]

rtupb_pre_num <- rtupb_pre[,sapply(rtupb_pre, function(x) is.numeric(x))]
vbpps_pre_num <- vbpps_pre[,sapply(vbpps_pre, function(x) is.numeric(x))]

# Barplots des variables catégoriques
fill <- c("#d3d3d3", "#a8a8a8", "#7e7e7e", "#545454", "#2a2a2a")

for (i in 1:ncol(rtupb_vbpps_pre)) {
    # Boxplots des variables numériques
    if (is.numeric(rtupb_vbpps_pre[,i]))
        boxplot(rtupb_vbpps_pre[,i], rtupb_pre[,i], vbpps_pre[,i],
            names=c("Both","RTUPB","VBPPS"),horizontal=TRUE,main=paste("RTUPB/VBPPS - Pré-opératoire\n(", colnames(rtupb_vbpps_pre)[i],")"), col=c("blue","red","#009900"))
    else
        if (colnames(rtupb_vbpps_pre)[i] != "Technique"){
            # Stacked barplots des variables catégoriques
            p4 <- ggplot() + geom_bar(aes(y = factor(1), x = Technique, fill = rtupb_vbpps_pre[,i]), data = rtupb_vbpps_pre,
                stat="identity") +
                theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
                scale_fill_manual(values=fill) + labs(x="", y="") + ggtitle(paste("RTUPB/VBPPS - Pré-opératoire\n(",colnames(rtupb_vbpps_pre)[i],")"))
            print(p4)
    }
}

boxplot(rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 1),'Age'], rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 0),'Age'], names=c("Porteur","Sans sonde"),horizontal=TRUE,main="RTUPB/VBPPS - Pré-opératoire\n(Age selon sonde)", col=c("blue","red"))
boxplot(rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 1),'IPSS'], rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 0),'IPSS'], names=c("Porteur","Sans sonde"),horizontal=TRUE,main="RTUPB/VBPPS - Pré-opératoire\nIPSS selon sonde", col=c("blue","red"))
boxplot(rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 1),'QoL'], rtupb_vbpps_pre[which(rtupb_vbpps_pre$Porteur_Sonde == 0),'QoL'], names=c("Porteur","Sans sonde"),horizontal=TRUE,main="RTUPB/VBPPS - Pré-opératoire\nQoL selon sonde", col=c("blue","red"))

table(rtupb_vbpps$Porteur_Sonde, rtupb_vbpps$Technique,dnn = c('Sonde','Technique'))
table(rtupb_vbpps$Porteur_Sonde, rtupb_vbpps$Indication,dnn = c('Sonde','Indication'))

mat_cor = cor(rtupb_vbpps_pre_num, method = c("pearson"))
adenome_heatmap(mat_cor)

adenome_scatterplot(rtupb_vbpps_pre_num,rtupb_vbpps_pre_num$Age,rtupb_vbpps_pre_num$IPSS,"Age","IPSS")

adenome_scatterplot(rtupb_vbpps_pre_num,rtupb_vbpps_pre_num$Duree_Traitement_Medical,rtupb_vbpps_pre$Delai_Ablation,"Durée Traitement Médical","Délai Ablation")
adenome_scatterplot(rtupb_vbpps_pre,rtupb_vbpps_pre$Duree_Traitement_Medical,rtupb_vbpps_pre$Volume_Prostatique,"Durée Traitement Médical","Volume Prostatique")

mat_cor = cor(rtupb_pre_num, method = c("pearson"))
adenome_heatmap(mat_cor)

adenome_scatterplot(rtupb_pre_num,rtupb_pre_num$Temps_Operation,rtupb_pre_num$Delai_Ablation,"Temps Opération","Délai Ablation")

mat_cor = cor(vbpps_pre_num, method = c("pearson"))
adenome_heatmap(mat_cor)

pca <- stats::prcomp(x=rtupb_vbpps_pre_num)
summary(pca)

sum(100 * (pca$sdev^2)[1:2] / sum(pca$sdev^2))

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 90))

library("factoextra")
theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()
loadings <- data.frame(pca$rotation, 
.names = row.names(pca$rotation))
p + geom_text(data=loadings, 
mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2")

fviz_pca_var(pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
             )

# on divise le dataset en 2 pour voir les techniques séparément
rtupb_post <- rtupb_vbpps_post[which(rtupb_vbpps_post$Technique == 1),]
vbpps_post <- rtupb_vbpps_post[which(rtupb_vbpps_post$Technique == 2),]

require(graphics)

# On choisit trois individus au hasard
individu = sample(1:nrow(rtupb_vbpps_post), 3)
i1=individu[1]
i2=individu[2]
i3=individu[3]

# On crée 2 individus moyens
imean_vbpps <- colMeans(subset(vbpps_post,select=-c(Technique)))
imean_rtupb <- colMeans(subset(rtupb_post,select=-c(Technique)))
# que l'on ajoute au dataset
mat <- rbind(rtupb_vbpps_post, c(1,t(imean_rtupb)), c(2,t(imean_vbpps)))

# On crée des séries temporels pour nos variable Q_max, IPSS et QoL
rtupb_vbpps_post_qmax=t(subset(mat,select=c(M1_Qmax,M3_Qmax,M6_Qmax,M9_Qmax,M12_Qmax,M15_Qmax,M18_Qmax)))
rtupb_vbpps_post_ipss=t(subset(mat,select=c(M1_IPSS,M3_IPSS,M6_IPSS,M9_IPSS,M12_IPSS,M15_IPSS,M18_IPSS)))
rtupb_vbpps_post_QoL=t(subset(mat,select=c(M1_QoL,M3_QoL,M6_QoL,M9_QoL,M12_QoL,M15_QoL,M18_QoL)))

rtupb_vbpps_post_series_qmax=ts(rtupb_vbpps_post_qmax, start=1/12,deltat=3/12)
rtupb_vbpps_post_series_ipss=ts(rtupb_vbpps_post_ipss, start=1/12,deltat=3/12)
rtupb_vbpps_post_series_QoL=ts(rtupb_vbpps_post_QoL, start=1/12,deltat=3/12)

# On plot l'évolution de Q_max pour les 3 individus tirés au sort et les individus moyens pour chaque technique
Q_max_Individu1=rtupb_vbpps_post_series_qmax[,i1]
Q_max_Individu2=rtupb_vbpps_post_series_qmax[,i2]
Q_max_Individu3=rtupb_vbpps_post_series_qmax[,i3]
Q_max_Individu_rtupb=rtupb_vbpps_post_series_qmax[,70]
Q_max_Individu_vbpps=rtupb_vbpps_post_series_qmax[,69]
Q_max=cbind(Q_max_Individu1,Q_max_Individu2,Q_max_Individu3,Q_max_Individu_rtupb,Q_max_Individu_vbpps)

ts.plot(Q_max,gpars= list(col=rainbow(6)),ylab= "Q_max",main="RTUPB/VBPPS - Post-opératoire\nEvolution de Q_max pour différents individus")
legend("topright", colnames(Q_max), col=1:ncol(Q_max), lty=1, cex=.65)

#On trace l'évolution de QMax en fonction du temps
Boxplot=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
for (i in 1:nrow(Boxplot)){
    Boxplot[i,2]=rtupb_vbpps_post$M1_Qmax[i]
    Boxplot[i,1]=1
}
Boxplot2=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
Qmax=subset(rtupb_vbpps_post,select=c(M3_Qmax,M6_Qmax,M9_Qmax,M12_Qmax,M15_Qmax,M18_Qmax))
for (i in 1:ncol(Qmax)){
for (k in 1:nrow(Boxplot2)){
    Boxplot2[k,2]=Qmax[k,i]
    Boxplot2[k,1]=i*3
}
Boxplot=rbind(Boxplot,Boxplot2)
}
boxplot(Boxplot[,2]~Boxplot[,1], main="RTUPB/VBPPS - Post-opératoire\nEvolution de Q_max en fonction du temps", col=6 ,xlab="Nombre de mois", ylab="Valeur de Q_max")

# On plot l'évolution de IPSS pour les 3 individus tirés au sort et les individus moyens pour chaque technique
IPSS_Individu1=rtupb_vbpps_post_series_ipss[,i1]
IPSS_Individu2=rtupb_vbpps_post_series_ipss[,i2]
IPSS_Individu3=rtupb_vbpps_post_series_ipss[,i3]
IPSS_Individu_rtupb=rtupb_vbpps_post_series_ipss[,70]
IPSS_Individu_vbpps=rtupb_vbpps_post_series_ipss[,69]
IPSS=cbind(IPSS_Individu1,IPSS_Individu2,IPSS_Individu3,IPSS_Individu_rtupb,IPSS_Individu_vbpps)

ts.plot(IPSS,gpars= list(col=rainbow(6)),ylab= "IPSS",main="RTUPB/VBPPS - Post-opératoire\nEvolution de IPSS pour différents individus")
legend("topright", colnames(IPSS), col=1:ncol(IPSS), lty=1, cex=.65)

cat ("Moyenne pré-opératoire - toute population / RTUPB / VBPPS : ...........", mean(rtupb_vbpps$IPSS), mean(rtupb_pre$IPSS), mean(vbpps_pre$IPSS),"\n")
cat ("Moyenne post-opératoire 1er mois - toute population /RTUPB / VBPPS : ..", mean(rtupb_vbpps$M1_IPSS), mean(rtupb_post$M1_IPSS), mean(vbpps_post$M1_IPSS),"\n")

#On trace l'évolution d'IPSS en fonction du temps
Boxplot=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
for (i in 1:nrow(Boxplot)){
    Boxplot[i,2]=rtupb_vbpps_post$M1_IPSS[i]
    Boxplot[i,1]=1
}
Boxplot2=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
IPSS=subset(rtupb_vbpps_post,select=c(M3_IPSS,M6_IPSS,M9_IPSS,M12_IPSS,M15_IPSS,M18_IPSS))
for (i in 1:ncol(IPSS)){
for (k in 1:nrow(Boxplot2)){
    Boxplot2[k,2]=IPSS[k,i]
    Boxplot2[k,1]=i*3
}
Boxplot=rbind(Boxplot,Boxplot2)
}
boxplot(Boxplot[,2]~Boxplot[,1], main="RTUPB/VBPPS - Post-opératoire\nEvolution de IPSS en fonction du temps", col=6 ,xlab="Nombre de mois", ylab="Valeur de IPSS")

# On plot l'évolution de QoL pour les 3 individus tirés au sort et les individus moyens pour chaque technique
QoL_Individu1=rtupb_vbpps_post_series_QoL[,i1]
QoL_Individu2=rtupb_vbpps_post_series_QoL[,i2]
QoL_Individu3=rtupb_vbpps_post_series_QoL[,i3]
QoL_Individu_rtupb=rtupb_vbpps_post_series_QoL[,70]
QoL_Individu_vbpps=rtupb_vbpps_post_series_QoL[,69]
QoL=cbind(QoL_Individu1,QoL_Individu2,QoL_Individu3,QoL_Individu_rtupb,QoL_Individu_vbpps)

ts.plot(QoL,gpars= list(col=rainbow(6)),ylab= "QoL",main="RTUPB/VBPPS - Post-opératoire\nEvolution de QoL pour différents individus")
legend("topright", colnames(QoL), col=1:ncol(QoL), lty=1, cex=.65)

#On trace l'évolution de QoL en fonction du temps
Boxplot=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
for (i in 1:nrow(Boxplot)){
    Boxplot[i,2]=rtupb_vbpps_post$M1_QoL[i]
    Boxplot[i,1]=1
}
Boxplot2=matrix(NA,nrow=nrow(rtupb_vbpps_post),ncol=2)
QoL=subset(rtupb_vbpps_post,select=c(M3_QoL,M6_QoL,M9_QoL,M12_QoL,M15_QoL,M18_QoL))
for (i in 1:ncol(QoL)){
for (k in 1:nrow(Boxplot2)){
    Boxplot2[k,2]=QoL[k,i]
    Boxplot2[k,1]=i*3
}
Boxplot=rbind(Boxplot,Boxplot2)
}
boxplot(Boxplot[,2]~Boxplot[,1], main="RTUPB/VBPPS - Post-opératoire\nEvolution de QoL en fonction du temps", col=6 ,xlab="Nombre de mois", ylab="Valeur de QoL")

mat <- subset(rtupb_vbpps_post, select = c(M1_Qmax,M1_IPSS,M3_Qmax,M3_IPSS,M6_Qmax,M6_IPSS,M9_Qmax,M9_IPSS,M12_Qmax,M12_IPSS,M15_Qmax,M15_IPSS,M18_Qmax,M18_IPSS))
mat_cor = cor(mat, method = c("pearson"))
adenome_heatmap(mat_cor)

# On supprime les 2 individus moyens de nos time series qui seront utilisées plus tard
rtupb_vbpps_post_series_qmax <- subset(rtupb_vbpps_post_series_qmax, select = -c(69,70))
rtupb_vbpps_post_series_ipss <- subset(rtupb_vbpps_post_series_ipss, select = -c(69,70))
rtupb_vbpps_post_series_QoL <- subset(rtupb_vbpps_post_series_QoL, select = -c(69,70))

require(cluster)
library(fpc)

# On ne prend pas en compte la variable Technique dans le clustering
df_pour_clustering <- subset(rtupb_vbpps_pre,select=-c(Technique))

# On crée notre matrice de dissimilarité utilisant daisy car matrice hétérogène
gower_dist = daisy(df_pour_clustering, metric = "gower")

summary(gower_dist)
gower_mat <- as.matrix(gower_dist)

## Tracé du dendrogramme
cah_pre_op <- hclust(gower_dist)

plot(as.dendrogram(cah_pre_op))

# Calculate silhouette width for many k using PAM
sil_width <- c(NA)

for(i in 2:10){
  pam_test <- pam(gower_dist, diss = TRUE, k = i)
  sil_width[i] <- pam_test$silinfo$avg.width
}

# Plot sihouette width (higher is better)

plot(1:10, sil_width, xlab = "Nombre de clusters", ylab = "Valeur Silhouette")
lines(1:10, sil_width)

k.patient.optimal = 3

pam.patient <- pam(gower_dist, diss = TRUE, k = k.patient.optimal, metric = "gower")
pam.patient$clusinfo

clusplot(pam.patient, main = paste("Cluster plot, k = ", k.patient.optimal), color = TRUE,lines=0,labels=5,span=TRUE,stand=FALSE)

df_pour_clustering[pam.patient$medoids,]

# Ajout du cluster aux datasets post-opératoire
rtupb_vbpps_post_clus <- data.frame(rtupb_vbpps_post, pam.patient$clustering)
attributes(rtupb_vbpps_post_clus)$names[23] <- 'patient.cluster'

rtupb_post_clus <- rtupb_vbpps_post_clus[which(rtupb_vbpps_post_clus$Technique == 1),]
vbpps_post_clus <- rtupb_vbpps_post_clus[which(rtupb_vbpps_post_clus$Technique == 2),]

for (i in 1:k.patient.optimal) {
    boxplot(rtupb_post_clus[which(rtupb_post_clus$patient.cluster == i),'M12_Qmax'],
            vbpps_post_clus[which(vbpps_post_clus$patient.cluster == i),'M12_Qmax'],
            names=c("RTUPB","VBPPS"),horizontal=TRUE,main=paste("M12 Qmax - Cluster ", i), col=c("blue","red"),
           ylim=c(0,55))
}

require(cluster)
library(fpc)

# On crée une matrice de séries temporelles
df_series = as.data.frame(cbind(t(rtupb_vbpps_post_series_ipss),t(rtupb_vbpps_post_series_QoL),t(rtupb_vbpps_post_series_qmax)))

# On crée notre matrice de dissimilarité avec distance euclidienne
euclidean_dist = dist(df_series, method = "euclidean")
attributes(df_series)$names = c("M1_IPSS","M3_IPSS","M6_IPSS","M9_IPSS","M12_IPSS","M15_IPSS","M18_IPSS",
                                "M1_QoL","M3_QoL","M6_QoL","M9_QoL","M12_QoL","M15_QoL","M18_QoL",
                                "M1_Qmax","M3_Qmax","M6_Qmax","M9_Qmax","M12_Qmax","M15_Qmax","M18_Qmax"
                               )

summary(euclidean_dist)
euclidean_mat <- as.matrix(euclidean_dist)

## Tracé du dendrogramme
cah_post_op <- hclust(euclidean_dist)
plot(as.dendrogram(cah_post_op))

# Calculate silhouette width for many k using PAM
sil_width <- c(NA)

for(i in 2:10){
  pam_test <- pam(euclidean_dist, diss = TRUE, k = i)
  sil_width[i] <- pam_test$silinfo$avg.width
}

# Plot sihouette width (higher is better)
plot(1:10, sil_width, xlab = "Nombre de clusters", ylab = "Valeur Silhouette")
lines(1:10, sil_width)

k.guerison.optimal = 2

pam.guerison <- pam(euclidean_dist, diss = TRUE, k = k.guerison.optimal, metric = "euclidean")
pam.guerison$clusinfo

clusplot(pam.guerison, main = paste("Cluster plot, k = ", k.guerison.optimal), color = TRUE,lines=0,labels=5,span=TRUE,stand=FALSE)

k.guerison.optimal = 3

pam.guerison <- pam(euclidean_dist, diss = TRUE, k = k.guerison.optimal, metric = "euclidean")
pam.guerison$clusinfo

clusplot(pam.guerison, main = paste("Cluster plot, k = ", k.guerison.optimal), color = TRUE,lines=0,labels=5,span=TRUE,stand=FALSE)

df_series[pam.guerison$medoids,1:7]
df_series[pam.guerison$medoids,8:14]
df_series[pam.guerison$medoids,15:21]

# Calcul de l'arbre

ipss = subset(rtupb_vbpps_post, select = c(M12_IPSS))
matrix_for_tree_ipss = cbind(rtupb_vbpps_pre, ipss)

tree_ipss = rpart(method = "anova", formula = matrix_for_tree_ipss$M12_IPSS ~.,data=matrix_for_tree_ipss,control=rpart.control(minsplit=5,cp=0))

pSimple <- prune(tree_ipss,cp=0)

# Plot de l'arbre
rpart.plot(pSimple)

# Evaluation de la performance de l'arbre avec 10 individus tirés au hasard
rnd_individus = sample(1:nrow(matrix_for_tree_ipss), 10)
df <- data.frame(predict(pSimple,matrix_for_tree_ipss[rnd_individus,]), matrix_for_tree_ipss[rnd_individus,15])
colnames(df) <- c("Prédiction","Observé")
df

# Calcul de l'arbre
qmax = subset(rtupb_vbpps_post, select = c(M12_Qmax))
matrix_for_tree_qmax = cbind(rtupb_vbpps_pre, qmax)

tree_qmax = rpart(method = "anova", formula = matrix_for_tree_qmax$M12_Qmax ~.,data=matrix_for_tree_qmax,control=rpart.control(minsplit=10,cp=0))

pSimple <- prune(tree_qmax,cp=0)

# Plot de l'arbre
rpart.plot(pSimple)

# Evaluation de la performance de l'arbre avec 10 individus tirés au hasard
rnd_individus = sample(1:nrow(matrix_for_tree_qmax), 10)
df <- data.frame(predict(pSimple,matrix_for_tree_qmax[rnd_individus,]), matrix_for_tree_qmax[rnd_individus,15])
colnames(df) <- c("Prédiction","Observé")
df

# Calcul de l'arbre
qol = subset(rtupb_vbpps_post, select = c(M12_QoL))
matrix_for_tree_qol = cbind(rtupb_vbpps_pre, qol)

tree_qol = rpart(method = "anova", formula = matrix_for_tree_qol$M12_QoL ~.,data=matrix_for_tree_qol,control=rpart.control(minsplit=5, cp=0))

pSimple <- prune(tree_qol,cp=0)

# Plot de l'arbre
rpart.plot(pSimple)

# Evaluation de la performance de l'arbre avec 10 individus tirés au hasard
rnd_individus = sample(1:nrow(matrix_for_tree_qol), 10)
df <- data.frame(predict(pSimple,matrix_for_tree_qol[rnd_individus,]), matrix_for_tree_qol[rnd_individus,15])
colnames(df) <- c("Prédiction","Observé")
df

# Calcul de l'arbre
Profil_Guerison <- pam.guerison$clustering
matrix_for_tree_guerison = cbind(rtupb_vbpps_pre, Profil_Guerison)

tree_guerison = rpart(method = "class", formula = matrix_for_tree_guerison$Profil_Guerison ~.,data=matrix_for_tree_guerison,control=rpart.control(minsplit=5,cp=0))

pSimple <- prune(tree_guerison,cp=0)

# Plot de l'arbre
rpart.plot(pSimple)

# Evaluation de la performance de l'arbre avec 10 individus tirés au hasard

rnd_individus = sample(1:nrow(matrix_for_tree_guerison), 10)

# Calcul de la prédiction
p <- predict(pSimple,matrix_for_tree_guerison[rnd_individus,])
# et facilitation de la lecture de prédiction
p[p == 1] <- 'X'

res <- cbind(matrix_for_tree_guerison[rnd_individus,15], p)
colnames(res) <- c("Observé", "Prédiction classe = 1", "Prédiction classe = 2", "Prédiction classe = 3")
res
