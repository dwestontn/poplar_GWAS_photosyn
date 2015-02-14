

#Starting with the duplicate analysis
data1<-read.csv("unique_dup.csv", stringsAsFactors= FALSE)
data2<-read.csv("clone_dup.csv", stringsAsFactors=FALSE)
head(data1)
tail(data1)
head(data2)
tail(data2)

str(data1)

data1


dim(data1)
dim(data2)

A<-data1[2:40, 2:27]
B<-data2[2:40, 2:27]

x<-diag(cor(A[sapply(A, is.numeric)], B[sapply(B, is.numeric)]))

#cnbf.28.4     dena.17.1     fnyi.28.3     glcb.26.1     homd.21.4     hopf.27.1     hopf.27.5     hopg.27.5     hrsp.27.4     klnd.20.3     klng.20.3     klng.20.7 
#0.8395983            NA            NA     0.9657183     0.9480020     0.9792961     0.7020884     0.9454716     0.9827976     0.9886425     0.9788342     0.9498837 
#mcmn.27.3     sqma.25.1     sqmb.25.4 chwh.27.1.csv  gw.11053.csv skwe.24.1.csv     whte.28.4     klng.20.4     klnc.20.2      gw.11027       besc.36      besc.881 
#0.9471502     0.8694467     0.9204827     0.9269767     0.9421240            NA     0.8406483     0.9725034            NA     0.9795961     0.9859390     0.9907436 
#besc.886      besc.183 
#0.9839980            NA 


############################################
######start of the phenotypic variations####
############################################

list.files()

location<-read.csv("BESC_LocationInfoAll.csv")
data1<-read.csv("R_phys_phenotypedata.csv")

head(data1)
names(data1)
head(location)

summary(location)

########### modify data so file string ids are same format
##skip by using R_phys_phenotypedata.csv file
#############

data1$id<-toupper(data1$id) #upper case
data1$id<-gsub("*.CSV", "", data1$id) #replace .csv to blank
head(data1)
str(data1)

##### merge files
data2<-merge(data1, location)

## not many ids mapped, due to luke not using same genotypes as wellington. 

#########################################
#####clustering our own kniship matix####
###########################################

kin<-read.csv("kinship_matrix_dec_2014.csv")

head(kin)

dissimilarity <- 1 - kin
distance <- as.dist(dissimilarity)

#we did not need to do this, but good practise
#kin <- na.omit(kin) # listwise deletion of missing
#mydata <- scale(mydata) # standardize variables 

kin <- na.omit(kin) # listwise deletion of missing
d <- dist(kin, method = "euclidean")
fit <- hclust(d, method="ward.D2") 
plot(fit)
groups <- cutree(fit, k=12) # cut tree into 12 clusters
rect.hclust(fit, k=12, border="green") # draw dendogram with red borders around the 5 clusters
str(fit)

## clustering on PCA
p = princomp(kin)







library(ggplot2)
qplot(Jmax25, data=data1)


p2 <- ggplot(xy, aes(x = xvar)) + geom_histogram(aes(y = ..density..), color = "black", fill = NA) + geom_density(color = "blue")

g3 <- ggplot(data1, aes(Jmax25)) + geom_histogram(aes(y = ), color = "blue",fill = NA, binwidth = 10) + theme_bw() + geom_rug(col = "darkred", alpha = 0.5) + geom_density(color = "green")
g3


#merged histogram/density plot
g3 <- ggplot(data1, aes(Jmax25)) + geom_histogram(aes(y = ..density..), color = "black",fill = NA, binwidth = 8) + theme_bw() + xlab("Jmax at 25 ˚C")+ geom_rug(col = "darkred", alpha = 0.5) + geom_density(color = "dark gray") + theme(axis.text = element_text(size = 18)) + theme(axis.title = element_text(size = 18))

g3 <- ggplot(data1, aes(Asat_ALight)) + geom_histogram(aes(y = ..density..), color = "black",fill = NA, binwidth = 2) + theme_bw() + xlab("Light Saturated CO2 Assimilation (Asat)")+ geom_rug(col = "darkred", alpha = 0.5) + geom_density(color = "dark gray") + theme(axis.text = element_text(size = 18)) + theme(axis.title = element_text(size = 18)) + xlim(0,30)
g3
#lets see if there is a relationship to lattitude

list.files()
lat<-read.csv("BESC_latitude_longitude.csv")
head(lat)
lat2<-lat[,1:2]
data2<-merge(data1,lat2)





#modify data so file string ids are same format
data1$id<-toupper(data1$id) #upper case
data1$id<-gsub("*.CSV", "", data1$id) #replace .csv to blank
head(data1)
head(data1)

qplot(Latitude,Jmax25, geom = c("point", "lm"), data=data2)
head(data2)


#use 500m or lower. 





head(data1)
#port 1 = S. mag, 2 = S. fallax, 3 & 4 = S. mag, 5 & 6 = S. fallax

#subset if wated
data_mid<-subset(data1, Binned.Time.of.Day=="Midday")
data_night<-subset(data1, Binned.Time.of.Day=="Night")





data_sub_1_2<-subset(data_outlier, Binned.Time.of.Day=="Midday" & Port.<=2)
data_sub_3__6<-subset(data1, Binned.Time.of.Day=="Midday" & Port.==3,6)

#using ggplot2
p1 <- ggplot(data1, aes(x = DFOY, y = Exp_Flux_Corr))
p1 + geom_point() + + geom_point(aes(y=), color="blue") (geom_smooth()
#p1+ geom_point() + facet_grid(.~Port.) + geom_smooth(method = lm)
#p1+ geom_point(color='steelblue', alpha=1/2) + facet_grid(.~Port.) + geom_smooth(method = lm)
#p1 + geom_point(aes(color=Port., alpha=1/2))
p1 + geom_point(aes(color=Port., alpha=1/2)) + theme_bw(base_family="Times")

p1 <- ggplot(data1, aes(x = DFOY, y = Exp_Flux_Corr))
p1 + geom_point(color="steelblue", alpha=0.3) + facet_wrap(Binned.Time.of.Day ~ species) + labs(title = "Binned by PAR & Species") + labs(x = "Day of Year") + labs(y = "CO2 Flux")

#scatterplot with mulitple axis
p1 <- ggplot(data_mid, aes(x = DFOY, y = Exp_Flux_Corr))
p1 + geom_point(aes(color = factor(species))) + geom_point(aes(y=Mean.Sph.Temp.0cm), color="steelblue")


#density and rug plot
ggplot(data1, aes(x=Exp_Flux_Corr, fill=species)) + geom_density(alpha=.3) + geom_rug(col = "darkred", alpha = 0.1)

#box plot on the shoulder season as denoted by temp quantiles
ggplot(data1, aes(temp_level, Exp_Flux_Corr)) + geom_boxplot(aes(fill = species)) + theme(legend.position = "none")


#if you want to keep outliers, don't just change ylim do this:
 #p1+geom_line() + coord_cartesian(ylim=c(low, max))

#try to facet by PAR level
p1 <- ggplot(data1, aes(x = DFOY, y = Exp_Flux_Corr))
p1 + geom_point(color="steelblue", alpha=0.3) + facet_wrap(Binned.Time.of.Day ~ temp_level) + labs(title = "Binned by PAR & Temp") + 
  labs(x = "Day of Year") + labs(y = "CO2 Flux")


> #using ggplot2
  > p1 <- ggplot(data_night, aes(x = DFOY, y = Exp_Flux_Corr))
> p1 + geom_point(aes(color = factor(species))) + geom_point(aes(y=Mean.Sph.Temp.0cm), color="steelblue", alpha=I(1/4)) + geom_point(aes(y=EM1_AirTC_2M), color="green",alpha=I(1/4))


#some intersting plots
qplot(DFOY, Exp_Flux_Corr, data=data1, color=Binned.Time.of.Day)

#facets
qplot(DFOY, Exp_Flux_Corr, data=data1, facets=.~Binned.Time.of.Day)

#density plots
qplot(Exp_Flux_Corr, data=data1, geom = "density", color = Binned.Time.of.Day)

#If you want to categorize a continuous vairable:
temp_cut<-quantile(data1$EM1_Hollow.0cm, seq(0,1, length = 6), na.rm=TRUE) 
head(temp_cut)
data1$temp_level<-cut(data1$EM1_Hollow.0cm, temp_cut)
levels(data1$temp_level)



#remove outliers
data_outlier<-data1
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}