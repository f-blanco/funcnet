######## Blanco et al. 2020 script


######Load the packages########
setwd("")
library(bipartite)
library(picante)
library(googlesheets)
library(hrbrthemes)
library(RColorBrewer)
library(viridis)
library(googlesheets4)
library(igraph)
library(gridExtra)
library(ggplot2)
library(openxlsx)
library(nlme)
library(paleoTS)
library(plotrix)
library(wesanderson)
library(munsell)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(doParallel)

options(stringsAsFactors = FALSE)
options(scipen=999)



####### Load and handling data #######

#####Select the sampling

Functional_entities <- T
Taxonomic <- T

MPO <- T
Quantile <- F
      quantile <- 0.75
Localities <- T
      loc <- 4
Bins <- F
      bin_l <- 0.5


occ <- read.xlsx("")

sites <- read.xlsx("")

fdata <- read.xlsx("")
rownames(fdata) <- paste(fdata$genus, fdata$species, sep="_")
fdata <- fdata[c("body_size_cat_corrected", "diet_cat","locomotion_cat")]
names(fdata) <- c("body_size_cat", "diet", "locomotion")
               

if(Taxonomic==T){
  ###We check that we have the same sites in the matrix and the sites database  
  sites_def <- sites$site
  occ <- occ[,sites_def]
 
  
  ###We use the occ as the matrix for the network analysis
  occ_analysis <- occ
}


if(MPO==T){
  ####Sampling one carnivore and 2 of 3 artiodactyl-perissodactyl-proboscidean
  
  ###Load the functional traits data
  order_sampling <- read.xlsx("")
  
  ##### Vector with the same order as occ with the names of the orders
  order_vector<- as.vector(order_sampling$order)
  
  ##### Change species names in occ by order names in the vector with aggregate
  names(order_vector) <- rownames(occ)
  
  order_occ <- aggregate(occ, by=list(order_vector),sum)
  rownames(order_occ)=order_occ$Group.1
  order_occ<- order_occ[,-1]
  
  ######## Take only sites with at least 1 carnivore and 2 of 3 inside artiodactyl-perissodactyl-proboscidean
  order_final <- order_occ[row.names(order_occ)=="Carnivora"|row.names(order_occ)=="Artiodactyla"|row.names(order_occ)=="Perissodactyla"|row.names(order_occ)=="Proboscidea",]
  order_final[order_final!=0] <- 1

  ####First with the 3 herbivores groups
  rest_order <- order_final[row.names(order_final)=="Artiodactyla"|row.names(order_final)=="Perissodactyla"|row.names(order_final)=="Proboscidea",]
  rest_order <- rest_order[colSums(rest_order)>=2]
  

  ####Second with the carnivores
  carnivora <- order_final[row.names(order_final)=="Carnivora",]
  carnivora <- carnivora[colSums(carnivora)!=0]

  
  #####Combine the two matrix to get the sites with  our condition for the sampling
  order_sites <- colnames(carnivora)[match(colnames(rest_order),colnames(carnivora))]
  order_sites <- order_sites[!is.na(order_sites)]


  ####we only take sites from the order_final
  occ_sampling_order<- occ[,order_sites]
  occ_sampling_order<- occ_sampling_order[colSums(occ_sampling_order)>4]
  ####remove rows that sum=0
  occ_sampling_order<- occ_sampling_order[rowSums(occ_sampling_order)!=0,]
  
  
  ###We add the current localities, because we know that they are well sample
  current_data <- occ[,colnames(occ)=="Barcelona" | colnames(occ)=="Sanlucar" ]

  
  occ_MPO_total <- merge(occ_sampling_order, current_data,by="row.names",all=T,sort=T)
  rownames(occ_MPO_total) <- occ_MPO_total[,1]
  occ_MPO_total <- occ_MPO_total[,-1]
  occ_MPO_total[is.na(occ_MPO_total)] <-0
  

  
  ######We take only sites in the MPO sampling
  sites <- sites[match(colnames(occ_order_actuales),sites$site),]
 
  
 ####We transform to get our matrix for the analysis
  occ_analysis <- occ_order_actuales

  
}


if(Quantile==T){
  
  
  ######Sampling quartile
  mid_ages <- (sites$first+sites$last)/2
  names(mid_ages) <- sites$site
  occ <- occ[,sites$site]
  
  bins <- as.numeric(cut(mid_ages, breaks=c(-1,0.001,0.0117,0.126,0.781,1.8,2.58,3.6,5.3,7.246,11.63,13.82,15.97,20.44,23.03))) # seq(0,1.2,0.2)
  names(bins) <-sites$site
  cbind(mid_ages,bins)
  
  
  
  quantiles <- aggregate(colSums(occ), by=list(bins), function(x){quantile(x,probs = quantile)})
  
  
  rich_locs_l <- list()
  b <- 14
  for(b in 1:14){
    occ_bin <- occ[,bins%in%b] 
    rich_locs_l[[b]] <- colnames(occ_bin)[which(colSums(occ_bin) >= quantiles[b,2])]
  }
  rich_locs_l
  
  quantil_sampling <- unlist(rich_locs_l)
  
  
  ####we only take sites from the quantile
  occ_quantile<- occ[,quantil_sampling]
  ####remove rows that sum=0
  occ_quantile<- occ_quantile[rowSums(occ_quantile)!=0,]
  
  ###We transformed into our data for the analysis
  occ_analysis <- occ_quantile
  
}


if(Localities==T){
  ######To check that we do not have localities with less than 5 species
  
  occ_analysis<- occ_analysis[colSums(occ_analysis)>loc]
  occ_analysis<- occ_analysis[rowSums(occ_analysis)!=0,]
  
}


if(Functional_entities==T){
  
  ##Check that we do not have NAs
  fdata <- fdata[(!is.na(fdata[,1]) & !is.na(fdata[,2]) & !is.na(fdata[,3])),] 

  ###Order the categories
  bs_ordered <- c("< 1","1-10","10-45","45-90","90-180","180-360","360-1000",">1000")
  fdata$body_size_cat <- ordered(fdata$body_size_cat, levels= bs_ordered)
  fdata$locomotion <- as.factor(fdata$locomotion)
  fdata$diet <- as.factor(fdata$diet)
  

  ####### To take only species in our selected sampling method
  fdata <- fdata[match(rownames(occ_analysis),rownames(fdata)),]
 
  ###To create the functional entities
  func_ent <- paste(fdata$body_size_cat,fdata$diet, fdata$locomotion, sep="_")
  names(func_ent) <- rownames(fdata)
  fe <- unique(func_ent)

  
  ######## Combine the matrix to have fe in the rows 
  occ_analysis <- aggregate(occ_analysis, by=list(func_ent),sum)
  rownames(occ_analysis)=occ_analysis$Group.1
  occ_analysis<- occ_analysis[,-1]
}

if(Bins==T){
  
  ####Create the bins from the sites information
  sites$mean_age <-(sites$first+sites$last)/2
  bin_lab <- seq(bin_l,ceiling(max(sites$mean_age)),bin_l)
  mid_ages_bins <- bin_lab-bin_l/2
  bin_lab <- paste(bin_lab-bin_l,"-",bin_lab)
  sites$bin <- as.numeric(cut(sites$mean_age,breaks=seq(0,ceiling(max(sites$mean_age)),bin_l),include.lowest = T))
  
  #####nº of localities per bin
  n_localities<- table(bin_lab[sites$bin])
  bin_lab <- unique(bin_lab[sites$bin])
  mid_ages_bins <- unique(mid_ages_bins[sites$bin])
  
  occ_fe_bins <- t(occ_analysis)
  bins_sites <- sites$bin[match(rownames(occ_fe_bins),sites$site)]
  occ_fe_bins <- aggregate(occ_fe_bins, by=list(bins_sites), sum)
  rownames(occ_fe_bins) <- occ_fe_bins$Group.1
  occ_fe_bins <- occ_fe_bins[,-1]
  rownames(occ_fe_bins) <- bin_lab
  occ_fe_bins <- t(occ_fe_bins)
  
  ###We transformed into our data for the analysis
  occ_analysis <- occ_fe_bins
  
  
  
}



#### Network analysis  ######



#####Load the function to run Infomap for the network analysis
mat.link.input<-function(mat){
  edg<-list()
  
  for(i in 1:ncol(mat)){
    edg[[i]]<-cbind(i,(which(mat[,i]>0)+ncol(mat)),mat[which(mat[,i]>0),i])
  }
  do.call(rbind,edg)
}

write.input.info<-function(mat,path,file){
  link.node<-mat.link.input(mat)
  vert<-cbind(c(1:(ncol(mat)+nrow(mat))),c(colnames(mat),rownames(mat)))
  input<-list(cbind("*Vertices",nrow(vert)),vert,cbind("*Edges",nrow(link.node)),link.node)
  write.table(input[[1]], paste(path,file,sep=""),row.names=F,col.names=F,quote=F )
  write.table(input[[2]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=c(2))
  write.table(input[[3]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F )
  write.table(input[[4]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F)
  return(input)
}

run.infomap <- function(path.info, path.in, name, extension, path.out, parameters){
  system(paste("bash -c", shQuote(paste("cd ", path.info ," && ./Infomap ", path.in, name , extension  ," ", path.out, " ", parameters, sep=""),type="sh")),wait=T)
  out.win<-gsub("/mnt/c","C:",path.out)
  par<-read.table(paste(out.win,name,".tree",sep=""))
  quality <- read.csv(paste(out.win,name,".tree",sep=""),nrows=1,sep="",check.names=F)
  quality <- as.numeric(names(quality)[grep("codelength",names(quality))+1])
  quality <- 1- quality [2]/quality [1]
  list(par,quality)
}

sam.val.fun <- function(x,path.unix,grep,treshold,Ntrain,Ntest,Nsam){
  path.win<-gsub("/mnt/c","C:",path.unix)
  setwd(path.win)
  write.table(x,"input_partitions.txt",sep=" ",row.names=F,col.names=F,quote=F)
  k<- paste(Ntrain,Ntest,Nsam,sep=" ")
  system(paste("bash -c", shQuote(paste("cd ",path.unix," && ./partition-validation -s ", runif(1,0,1000)," -t ", treshold, " --validation-sampling ", k, " input_partitions.txt ", grep, sep=""),type="sh")),wait=T)
  read.table(paste(grep,"_validation",sep=""))
}

jac.fun<-function(x,y){sum(x%in%y)/length(unique(c(x,y)))}


# Input #
# ::::: #
file <-  ""
path <-  ""
write.input.info(occ_analysis,path,file)


# To run infomap #
# :::::::::::::: #
path.info<-"" 
path.in<-"" 
name<-""
extension<-""
path.out<-""


par <- list()
q <- c()
for( i in 1:100){
  seed <- sample(1:1000,1)
  parameters<-paste("-N 100", "-s",seed)
  res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
  par[[i]] <- res.i[[1]]
  q[[i]] <- res.i[[2]]
  print(rep(i,100))
}


##### Robustness of the module searching, using "partition-valdiation" from: https://github.com/mapequation/partition-validation
# ::::::::::::::::::::: #

par.c <- lapply(par,function(x)as.character(x[order(x[,3]),1]))
par.c <- do.call(cbind,lapply(par.c,function(x)sapply(strsplit(x,":"),function(x)x[1])))
par.c <- par.c[,order(1-q)]


path.unix <- "" # route to the "partition-validation" software
grep <- paste("samplig_val",i,sep="_")
Ntrain <- 75
Ntest <- 25
Nsam <- 100
p <- sam.val.fun(par.c,path.unix,grep,0.25,Ntrain,Ntest,Nsam)
mean(p[,1]/Ntest)	
##More than 0.9 indicates that the search for modules is complete

####Modules robustness

id.nodes <- as.character(par[[1]][order(par[[1]][,3]),3])
mod.ref <- tapply(id.nodes,par.c[,1],list)
mod.alt <- apply(par.c[,-1],2,function(x)tapply(id.nodes,x,list))
id.mod <- names(mod.ref)

#similarity to the most similar
#::::::::::::::::::::::::::::::

jac.max<-list()
for(i in id.mod){
  sim.avg<-c()
  for(j in 1:length(mod.alt)){
    sim <- unlist(lapply(mod.alt[[j]],function(x){jac.fun(mod.ref[[i]],x)}))
    sim.avg[j] <- max(sim)
    names(sim.avg)[j]	<- names(sim)[sim==max(sim)]
    print(c(i,j))
  }
  jac.max[[i]]<- sim.avg
}



#proportion more similar than a threshold 0.5
#::::::::::::::::::::::::::::::::::::::::
# p.sim is the modules robustness

p.sim <- sapply(jac.max,function(x)sum(x>=0.5)/length(x)) 



#Node probability module assignation
#::::::::::::::::::::::::::
#res.prob is the probability that a node is asigned to is module

res.prob <- list()
for(i in 1:length(jac.max)){
  id.mod.i <- names(jac.max[[i]])
  nodes.mod.i <-c()
  for(j in 1:length(id.mod.i)){
    nodes.mod.i<-c(nodes.mod.i,mod.alt[[j]][[id.mod.i[j]]])
  }
  prob <- data.frame(table(nodes.mod.i)/length(id.mod.i))			
  res.prob[[i]] <- prob[prob[,1]%in% mod.ref[[names(jac.max)[i]]],]
  print(i)
}
names(res.prob) <- names(jac.max)




# Affinity, Specificity and indval index #
# :::::::::::::::::::::::::::::::: #


par <- cbind(id.nodes,par.c[,1])
id.gf <- rownames(occ_analysis)
id.lc <- colnames(occ_analysis)
par.gf <- par[par[,1]%in%id.gf,]
par.lc <- par[par[,1]%in%id.lc,]


specificity.lc <- c()
afinity.lc <- c()
for(i in 1:nrow(par.lc)){
  mod.i <- par.lc[i,2]
  gf.mod.i <- as.character(par.gf[par.gf[,2]==mod.i,1])
  lc.mod.i <- as.character(par.lc[par.lc[,2]==mod.i,1])
  val.mod <- sum(occ_analysis[rownames(occ_analysis)%in% gf.mod.i,colnames(occ_analysis)==par.lc[i,1]])
  val.tot.lc.i <- sum(occ_analysis[,colnames(occ_analysis)==par.lc[i,1]])
  val.mod.gf <- sum(occ_analysis[rownames(occ_analysis)%in% gf.mod.i,colnames(occ_analysis)==par.lc[i,1]]>0)
  specificity.lc[i]<-val.mod/val.tot.lc.i
  afinity.lc [i] <- val.mod.gf/length(gf.mod.i)
}


specificity.gf <- c()
afinity.gf <- c()
for(i in 1:nrow(par.gf)){
  mod.i <- par.gf[i,2]
  gf.mod.i <- as.character(par.gf[par.gf[,2]==mod.i,1])
  lc.mod.i <- as.character(par.lc[par.lc[,2]==mod.i,1])
  val.mod <- sum(occ_analysis[rownames(occ_analysis)==par.gf[i,1] ,colnames(occ_analysis)%in%lc.mod.i])
  val.tot.gf.i <- sum(occ_analysis[rownames(occ_analysis)==par.gf[i,1],])
  val.mod.lc <- sum(occ_analysis[rownames(occ_analysis)==par.gf[i,1],colnames(occ_analysis)%in%lc.mod.i]>0)
  specificity.gf[i]<-val.mod/val.tot.gf.i
  afinity.gf [i] <- val.mod.lc/length(lc.mod.i)
}

res_fin <- data.frame(ID = c(as.character(par.lc[,1]),as.character(par.gf[,1])),specificidad = c(specificity.lc,specificity.gf), afinidad= c(afinity.lc,afinity.gf))
res_fin$indval <- res_fin[,2]* res_fin[,3]
res_fin[,2:4]<-apply(res_fin[,2:4],2,as.numeric)



# Final table with the data for the plot  #
# :::::::::::::::::::::::::::::::: #


res.prob<-do.call(rbind,res.prob)
colnames(res.prob) <- c("ID","P_class")
colnames(par)<-c("ID","module")
RES <- merge(par,res_fin,by="ID")
RES <- merge(RES,res.prob,by="ID")

#############Plots for the different sampling methods##########

Plot_Bins <- F

Plot_sampling <- T


if(Plot_Bins==T){
  
  mod_loc <- RES[RES$ID %in% colnames(occ_analysis),]
  sites_bins <- as.data.frame(colnames(occ_analysis))
  names(sites_bins) <- c("ID")
  total_mid_ages_bins <- mid_ages_bins
  sites_bins$mean_age <- total_mid_ages_bins

  mod_fe <- RES[RES$ID %in% rownames(occ_analysis),]
  
  sites_bins$module <- RES$module[match(sites_bins$ID,RES$ID)]
  sites_bins$indval <- RES$indval[match(sites_bins$ID,RES$ID)]
  sites_bins$sp_richness <- colSums(occ_analysis)[match(sites_bins$ID, colnames(occ_analysis))]


  order_age <- aggregate(sites_bins$mean_age,by=list(sites_bins$module), median)
  RES$module <- match(RES$module,order(order_age$x, decreasing = T))
  sites_bins$module <- match(sites_bins$module,order(order_age$x, decreasing = T))
  
  
  mod_fe$module <- match(mod_fe$module,order(order_age$x, decreasing = T))

  mod_fun_ent_site <- matrix(0,nrow = nrow(sites_bins),ncol = max(unique(sites_bins$module)))
  mod_fun_ent_site <- mod_fun_ent_site/sites_bins$sp_richness
  
  y_position <- sites_bins$module*10-5
  y_position <- y_position+((sites_bins$module_fuc_ent)/(sites_bins$sp_richness))*10
  y_position_2 <- rep(1:2*10-5, each=nrow(sites_bins))+as.vector(mod_fun_ent_site)*8+1
  

  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  y_position_indval <-  sites_bins$module*10-5+range01(sites_bins$indval, na.rm=T)*8+1
  
  
  module_colours <- c("")
  
  
  module_colours_2 <- module_colours[sites_bins$module]
  
  
  
  modFE <- as.numeric(mod_fe$module)
  
  inner_cex <- log(modFE)/3+0.2
  outer_cex <- log(sites_bins$sp_richness)*1.5#/3+0.3*12
  
  lwd=2
  line <- 2
  tcl=0.15
  line.ylabs <- -0.25
  
  
  
  pch=16
  
  
  
  n_module <- max(sites_bins$module)
  y.lab_line <- 1.3
  cex.axis <- 1
  col="#313131"
  col_var <- "3414141"
  point.transp <- 0.1
  
  
  quartz(width = 16, height = 12)
 
  plot(NULL, xlim=rev(range(sites_bins$mean_age)),ylim = c(5,(n_module+1)*10), axes = F, ylab = "", xlab="")

  abline(h=1:n_module*10-5, col=module_colours, lty=1, lwd=0.5)
  
  points(sites_bins$mean_age,y_position_indval, pch=pch, col=alpha(module_colours_2,0.6), cex=outer_cex)
  axis(1,labels=F,cex.axis=cex.axis,col=col,col.axis=col, las=1, line=0, tcl=tcl)
  axis(1,cex.axis=cex.axis,col=col,col.axis=col, las=1, line= line.ylabs-0.25, tick=F,tcl=tcl)
  
  axis(2,at=1:n_module*10,labels = c(paste("M",c(1:n_module), sep="")),cex.axis=cex.axis,col=module_colours,col.axis=col, las=1, line= line.ylabs, tick=F,tcl=tcl)
  mtext("time (Ma)",side=1,col= col,line= 3,cex=cex.axis, font=2)

  
}



if(Plot_sampling==T){
  
  
  mod_loc <- RES[RES$ID %in% sites$site,]
  

  sites <- sites[sites$site %in% mod_loc$ID,]
  sites[order(sites$module, decreasing = T),]

  mod_fe <- RES[RES$ID %in% rownames(occ_analysis),]

  sites$module <- RES$module[match(sites$site,RES$ID)]
  sites$indval <- RES$indval[match(sites$site,RES$ID)]
  sites$sp_richness <- colSums(occ_analysis)[match(sites$site, colnames(occ_analysis))]

  sites$mean_age <-(sites$first+sites$last)/2
  order_age <- aggregate(sites$mean_age,by=list(sites$module), median)

  RES$module <- match(RES$module,order(order_age$x, decreasing = T))
  
  sites$module <- match(sites$module,order(order_age$x, decreasing = T))
  
  
  mod_fe$module <- match(mod_fe$module,order(order_age$x, decreasing = T))
  

  mod_fun_ent_site <- matrix(0,nrow = nrow(sites),ncol = max(unique(sites$module)))
  

  mod_fun_ent_site <- mod_fun_ent_site/sites$sp_richness
  

  y_position <- sites$module*10-5
  y_position <- y_position+((sites$module_fuc_ent)/(sites$sp_richness))*10
  y_position_2 <- rep(1:2*10-5, each=nrow(sites))+as.vector(mod_fun_ent_site)*8+1

  
  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  y_position_indval <-  sites$module*10-5+range01(sites$sp_richness, na.rm=T)*8+1
  

  module_colours <- c()
  
  module_colours_2 <- module_colours[sites$module]
  
  unique(module_colours_2)
  
  modFE <- as.numeric(mod_fe$module)
  
  inner_cex <- log(modFE)/3+0.2
  outer_cex <- log(sites$sp_richness)*1.3#/3+0.3*12
  
  
  lwd=2
  line <- 2
  tcl=0.15
  line.ylabs <- -0.25
  
  pch=16
  
  
  n_module <- max(sites$module)
  y.lab_line <- 1.3
  cex.axis <- 1
  col="#313131"
  col_var <- "3414141"
  point.transp <- 0.1
  
  
  
  quartz(width = 16, height = 12)

  plot(NULL, xlim=rev(range(sites$mean_age)),ylim = c(5,(n_module+1)*10), axes = F, ylab = "", xlab="")
  abline(h=1:n_module*10-5, col=module_colours, lty=1, lwd=0.5)
  points(sites$mean_age,y_position_indval, pch=pch, col=alpha(module_colours_2,0.6), cex=outer_cex)
  axis(1,labels=F,cex.axis=cex.axis,col=col,col.axis=col, las=1, line=0, tcl=tcl)
  axis(1,cex.axis=cex.axis,col=col,col.axis=col, las=1, line= line.ylabs-0.25, tick=F,tcl=tcl)
  axis(2,at=1:n_module*10,labels = c(paste("M",c(1:n_module), sep="")),cex.axis=cex.axis,col=module_colours,col.axis=col, las=1, line= line.ylabs, tick=F,tcl=tcl)
  mtext("time (Ma)",side=1,col= col,line= 3,cex=cex.axis, font=2)
  
  
}




###### Functional richness and functional diversity########


fR <- vector("list",length(unique(sites$module)))
spR <- vector("list",length(unique(sites$module)))
fdiv <- vector("list",length(unique(sites$module)))


occ_density <- occ_analysis
occ_density[occ_density!=0] <- 1



for(m in 1:length(unique(sites$module))){
  
  module_sites<- sites$site[sites$module==m]
  module_fe <- mod_fe$ID[mod_fe$module==m]
  
  
  
  if(length(module_sites)==1){
    fe_module <- occ_density[module_sites]
    fe_total <- fe_module[rownames(fe_module) %in% module_fe,]
    occ_module <- occ_analysis[module_sites]
    occ_module_total <- occ_module[rownames(occ_module) %in% module_fe,]
    freqs_1fe <- length(occ_module_total)
    
    
    corrected_sp_richness <- sum(occ_module_total)
    corrected_fe_richness <- sum(fe_total)
    corrected_FDiv<- sum(table(corrected_fe_richness)/length(corrected_fe_richness) * log(table(corrected_fe_richness)/length(corrected_fe_richness))) * -1
    
    
  }else{
    fe_module <- occ_density[,module_sites]
    fe_total <- fe_module[rownames(fe_module) %in% module_fe,]
    occ_module <- occ_analysis[module_sites]
    occ_module_total <- occ_module[rownames(occ_module) %in% module_fe,]
    freqs_fe <- apply(occ_module_total,2,table)
    freqs_1fe <- unlist(lapply(freqs_fe,function(x) x["1"]))
    freqs_1fe[is.na(freqs_1fe)] <- 0
    
    
    
    corrected_sp_richness <- colSums(occ_module_total)
    corrected_fe_richness <- colSums(fe_total)
    corrected_FDiv<- sum(table(corrected_fe_richness)/length(corrected_fe_richness) * log(table(corrected_fe_richness)/length(corrected_fe_richness))) * -1
    
    
  }
  
  
  sp_rich <- corrected_sp_richness
  feRIch<- corrected_fe_richness
  feRedun<- corrected_redundancy
  fediv <- corrected_FDiv
  
  names(sp_rich) <- module_sites
  names(feRIch) <- module_sites
  names(fediv) <- module_sites

  
  spR[[m]] <- sp_rich
  fR[[m]]<- feRIch
  fDiv[[m]]<- fediv
  
  
  
}


measures_violin <- do.call(rbind, Map(data.frame, sp_richness=spR,fe_richness=fR, fe_diversity=fDiv))


measures_violin$module <- as.factor(sites$module[match(rownames(measures_violin),sites$site)])
measures_violin$mid_age <- as.factor(sites$mean_age[match(rownames(measures_violin),sites$site)])

violin_plots <- measures_violin[measures_violin$module=="3" | measures_violin$module=="4" | measures_violin$module=="6",]



####Plots


mod_cols <- c("#EE8959", "#E9C46A","#2A9D8F")

function_for_plot <- function(x){
  y <- median(x)
  ymin=hdi(x, 0.95)[1]
  ymax=hdi(x, 0.95)[2]
  res <- c(y=y, ymin=ymin ,ymax=ymax)
  names(res) <- c("y", "ymin","ymax")
  return(res)
}

function_ymin <- function(x){
  ymin=hdi(x, 0.95)[1]
  return(ymin)
}

function_ymax <- function(x){
  ymax=hdi(x, 0.95)[2]
  return(ymax)
}


quartz(height = 3.5, width = 5)

violin_fe_richness <- ggplot(violin_plots, aes(x=violin_plots$module, y=violin_plots$fe_rich, fill=module)) +
  geom_violin(trim=FALSE, colour=F, alpha=0.4) +
  labs(title="",x="module", y = expression("fe richness")) +
  theme_minimal() +
  scale_fill_manual(values=mod_cols) +
  stat_summary(fun.y="median",fun.ymin=function_ymin, fun.ymax = function_ymax,geom = "pointrange", colour=mod_cols,show.legend = F) 
violin_fe_richness 

quartz(height = 3.5, width = 5)

violin_fdiv <- ggplot(violin_plots, aes(x=violin_plots$module, y=violin_plots$fe_diversity, fill=module)) +
  geom_violin(trim=FALSE, colour=F, alpha=0.4) +
  labs(title="",x="module", y = expression("fe richness")) +
  theme_minimal() +
  scale_fill_manual(values=mod_cols) +
  stat_summary(fun.y="median",fun.ymin=function_ymin, fun.ymax = function_ymax,geom = "pointrange", colour=mod_cols,show.legend = F) 
violin_fdiv 


require(gridExtra)
quartz(height = 3, width = 10)
par(mfcol=c(1,1))

g <- grid.arrange(violin_fe_richness + theme_minimal() + stat_compare_means(label = "p.format",ref.group = ".all."),
                  violin_fdiv + theme_minimal() + stat_compare_means(method = "anova") + stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all."),
                  ncol=2)

######Null model#########


random_q <- list()

for (n in 1:100) {
  
  occ_null <- aggregate(occ_analysis, by=list(sample(func_ent)),sum)
  rownames(occ_null)=occ_null$Group.1
  occ_null<- occ_null[,-1]
  
  
  
  file <-  ""
  path <-  ""
  write.input.info(occ_null,path,file)
  
  
  # run infomap #
  # :::::::::::::: #
  path.info<-"" 
  path.in<-"" 
  name<-""
  extension<-""
  path.out<-""
  
  
  par <- list()
  q <- c()
  for( i in 1:100){
    seed <- sample(1:1000,1)
    parameters<-paste("-N 100", "-s",seed)
    res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
    par[[i]] <- res.i[[1]]
    q[[i]] <- res.i[[2]]
    print(rep(i,100))
  }
  
  random_q[[n]] <- q
  
}












####### Sensitivity analysis ########

### load Infonmap functions

mat.link.input<-function(mat){
  edg<-list()
  
  for(i in 1:ncol(mat)){
    edg[[i]]<-cbind(i,(which(mat[,i]>0)+ncol(mat)),mat[which(mat[,i]>0),i])
  }
  do.call(rbind,edg)
}

write.input.info<-function(mat,path,file){
  link.node<-mat.link.input(mat)
  vert<-cbind(c(1:(ncol(mat)+nrow(mat))),c(colnames(mat),rownames(mat)))
  input<-list(cbind("*Vertices",nrow(vert)),vert,cbind("*Edges",nrow(link.node)),link.node)
  write.table(input[[1]], paste(path,file,sep=""),row.names=F,col.names=F,quote=F )
  write.table(input[[2]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=c(2))
  write.table(input[[3]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F )
  write.table(input[[4]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F)
  return(input)
}

run.infomap <- function(path.info, path.in, name, extension, path.out, parameters){
  system(paste("bash -c", shQuote(paste("cd ", path.info ," && ./Infomap ", path.in, name , extension  ," ", path.out, " ", parameters, sep=""),type="sh")),wait=T)
  out.win<-gsub("/mnt/c","C:",path.out)
  par<-read.table(paste(out.win,name,".tree",sep=""))
  quality <- read.csv(paste(out.win,name,".tree",sep=""),nrows=1,sep="",check.names=F)
  quality <- as.numeric(names(quality)[grep("codelength",names(quality))+1])
  quality <- 1- quality [2]/quality [1]
  list(par,quality)
}

sam.val.fun <- function(x,path.unix,grep,treshold,Ntrain,Ntest,Nsam){
  path.win<-gsub("/mnt/c","C:",path.unix)
  setwd(path.win)
  write.table(x,"input_partitions.txt",sep=" ",row.names=F,col.names=F,quote=F)
  k<- paste(Ntrain,Ntest,Nsam,sep=" ")
  system(paste("bash -c", shQuote(paste("cd ",path.unix," && ./partition-validation -s ", runif(1,0,1000)," -t ", treshold, " --validation-sampling ", k, " input_partitions.txt ", grep, sep=""),type="sh")),wait=T)
  read.table(paste(grep,"_validation",sep=""))
}

jac.fun<-function(x,y){sum(x%in%y)/length(unique(c(x,y)))}





##################################Loop for the sensitivity analysis 

##Load the data for the sensitivity analysis
occ_analysis <- read.xlsx("",rowNames = T)

fdata <- read.xlsx("",rowNames = T)


#####Create the loop to do the sensitivity analysis 

diet_error <- seq(0,0.1,by=0.01)
locomotion_error <- seq(0,0.1,by=0.01)
randomizations <- 100

###Potencial categories of change for diet and locomotion following a biological meaning 

option_carnivores <- c("carnivore_invert","omnivore","hypercarnivore")
option_hypercar <- c("meat_bone","carnivore")
option_meat_bone <- c("hypercarnivore")
option_car_invert <- c("carnivore","omnivore")
option_omnivore <- c("carnivore","carnivore_invert","frugivore","browser")
option_frugivore <- c("omnivore","browser")
option_browser <- c("omnivore","mixed_feeder")
option_mixed_feeder <- c("browser","grazer")
option_grazer <- c("mixed_feeder")

diet_option <- list(carnivore=option_carnivores,hypercarnivore=option_hypercar,meat_bone=option_meat_bone,carnivore_invert=option_car_invert,omnivore=option_omnivore,frugivore=option_frugivore,browser=option_browser,mixed_feeder=option_mixed_feeder,grazer=option_grazer)

option_cursorial <- c("ambulatory")
option_ambulatory <- c("cursorial","aquatic","scansorial")
option_aquatic <- c("ambulatory")
option_scansorial <- c("ambulatory","arboreal")
option_arboreal <- c("scansorial")

locomotion_option <- list(cursorial=option_cursorial,ambulatory=option_ambulatory,aquatic=option_aquatic,scansorial=option_scansorial,arboreal=option_arboreal)


#####Original functional entities
original_func_ent <- paste(fdata$body_size_cat,fdata$diet, fdata$locomotion, sep="_")
names(original_func_ent) <- rownames(fdata)
fe_original <- unique(original_func_ent)


##### Loop to change randomly, with a predefined error, the functional categories of the species 
diet_list <- list()
locomotion_list <- list()
ran_list <- list()
spp_change_list <- list()

for (ran in 1:randomizations){
  
  for (diet_loop in 1:length(diet_error)) {
    
    
    
    for (locomotion_loop in 1:length(locomotion_error)) {
      
      diet_error_value <- diet_error[diet_loop]
      spp_changed_diet <-  rownames(fdata)[runif(length(rownames(fdata)))<diet_error_value]
      locmotion_error_value <- locomotion_error[locomotion_loop]
      spp_changed_loc <-  rownames(fdata)[runif(rownames(fdata))<locmotion_error_value]
      
      fdata_error <- fdata
      
      
      if(length(spp_changed_diet)!=0){
        for (sp_diet in 1:length(spp_changed_diet)) {
          sp_new_value_diet <- spp_changed_diet[sp_diet]
          current_diet_values<- fdata$diet[rownames(fdata) %in% sp_new_value_diet]
          new_diet_value<- sample(diet_option[[current_diet_values]],1)
          fdata_error$diet[rownames(fdata_error) %in% sp_new_value_diet]  <- new_diet_value
        }
      }
      
      
      if(length(spp_changed_loc)!=0){
        for (sp_loc in 1:length(spp_changed_loc)) {
          sp_new_value_loc <- spp_changed_loc[sp_loc]
          current_locomotion_values<- fdata$locomotion[rownames(fdata) %in% sp_new_value_loc]
          new_locomotion_value<- sample(locomotion_option[[current_locomotion_values]],1)
          fdata_error$locomotion[rownames(fdata_error) %in% sp_new_value_loc]  <- new_locomotion_value
        }
      }
      
      
      
      
      
      func_ent <- paste(fdata$body_size_cat,fdata_error$diet, fdata_error$locomotion, sep="_")
      names(func_ent) <- rownames(fdata_error)
      occ_sensitivity <- aggregate(occ_analysis, by=list(func_ent),sum)
      rownames(occ_sensitivity)=occ_sensitivity$Group.1
      occ_sensitivity<- occ_sensitivity[,-1]
      
      
      file_name <-  paste(diet_error_value,"_",body_size_error_precentage,"_",locmotion_error_value,"_",ran,"_","Infomap.net",sep = "")
      path_folder <-  ""
      write.input.info(occ_sensitivity,file=file_name,path=path_folder)
      
      
      path.info<-"" 
      path.in<-""
      name<-paste(diet_error_value,"_",body_size_error_precentage,"_",locmotion_error_value,"_",ran,"_","Infomap",sep = "")
      extension<-".net"
      path.out<-""
      
      
      
      spp_change_percentage <- (sum(func_ent!= original_func_ent)/length(original_func_ent))*100
      fe_nuevas <- unique(func_ent)
      fe_change <- length(fe_nuevas)-length(fe_original)
      
      diet_list <- c(diet_list,diet_error_value)
      ran_list <- c(ran_list,ran)
      locomotion_list <- c(locomotion_list,locmotion_error_value)
      spp_change_list <- c(spp_change_list,spp_change_percentage)
      fe_change_list <- c(fe_change_list,fe_change)
      
      
      
      
      seed <- sample(1:1000,1)
      parameters<-paste("-N 100", "-s",seed)
      res.i <- run.infomap (path.info, path.in, name, extension, path.out, parameters)
      
      
    }
    
    
  }
  
  
}

###Create a data frame for the plots with the info of all rans 
change_matrix <- do.call(rbind, Map(data.frame, locomotion=locomotion_list,diet=diet_list,randomization=ran_list,percentage_change=spp_change_list))
change_matrix <- aggregate(change_matrix,by=list(change_matrix$diet,change_matrix$locomotion),mean)







####################################Calculate similarity between the results of the sensitivity analysis and the original network





##### Jaccard function to calculate the similaraty of functional faunas composition


jac_mul <- function(p1,p2){
  p1l<-tapply(p1[,1],p1[,2],list)
  p2l<-tapply(p2[,1],p2[,2],list)
  
  cp1<-c()
  for(i in 1:length(p1l)){
    com <- unlist(lapply(p2l,function(x)sum(x%in%p1l[[i]])))
    tot <- unlist(lapply(p2l,function(x)length(unique(c(x,p1l[[i]])))))
    cp1[i]<-max(com/tot)
  }
  
  cp2<-c()
  for(i in 1:length(p2l)){
    com <- unlist(lapply(p1l,function(x)sum(x%in%p2l[[i]])))
    tot <- unlist(lapply(p1l,function(x)length(unique(c(x,p2l[[i]])))))
    cp2[i]<-max(com/tot)
  }
  
  p1_to_p2 <- mean(cp1)
  p2_to_p1 <- mean(cp2)
  jacmul   <- mean(	c(p1_to_p2,p2_to_p1) )
  
  return(cbind(p1_to_p2,p2_to_p1,jacmul))
}


##### Load the original network data 


par.1 <- read.xlsx("")
par.1 <- par.1[,c("Node.ID","module")]
par.1 <- par.1[par.1$Node.ID %in% colnames(occ_analysis),]
colnames(par.1) <- c("V3","p2")
###Chose the modules that form the functional faunas in the original network
par.1 <- par.1[par.1$p2=="3" | par.1$p2=="4" | par.1$p2=="8",]



files<- list.files("")


######Select the values for the error and rans to create a dataframe
diet<-  sapply(strsplit(files,"_"),function(x)x[1])
loc<- sapply(strsplit(files,"_"),function(x)x[3])
ran<- sapply(strsplit(files,"_"),function(x)x[4])



####Create a dataframe to introduce the results of the analysis
sensitivity_matrix <- data.frame(diet=diet,loc=loc,ran=ran,sim_network=NA)


####Loop to calculate the similarity between the original network partition and all the files from the sensitivity analysis
for (i in 1:length(files)) {
  
  par.files  <- read.table(paste("",files[i],sep=""))
  par.files$p2 <- sapply(strsplit(par.files[,1],":"),function(x)x[1])
  pk <- par.files[par.files[,3] %in% par.1$V3,c("V3","p2")]
  sensitivity_matrix[i,4]<- jac_mul(par.1,pk)[3]
  
}





####Create a data frame for the plots

sensitivity_matrix_heat_map <- aggregate(sensitivity_matrix,by=list(sensitivity_matrix$diet,sensitivity_matrix$loc),max)
error_diet_plot<- as.numeric(sensitivity_matrix_heat_map$Group.1)*100
error_loc_plot<- as.numeric(sensitivity_matrix_heat_map$Group.2)*100





#################################Plots sensitivity analysis

####### Heatmaps plots

change_matrix$Group.1 <- change_matrix$Group.1*100 ##diet
change_matrix$Group.2 <- change_matrix$Group.2*100 ##locomotion


##To create a custimezed palette
col_1 <- c("#facf5a","#f85959","#7c203a")
new_palette <- colorRampPalette(col_1,bias=0.7)
col <- new_palette(1000)


col_2 <- c("#fcbf1e","#40bad5","#035aa6")
new_palette <- colorRampPalette(col_2,bias=0.3)
col_3 <- rev(new_palette(1000))
transparency <- 0.8
font_type <- 1
font_family <- "sans"

theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 12),
          axis.line.x = element_line(color="#525252"), 
          axis.line.y = element_line(color="#525252"),
          panel.border = element_blank(),
          panel.grid.major.x = element_line(colour = alpha("#525252",0.3),size=0.2),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(colour = alpha("#525252",0.3),size=0.2),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 14, vjust = 2, hjust = 0.5),
          legend.text = element_text(size = 10),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "#525252", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

data <- data.frame(sp_error=change_matrix$percentage_change,similarity_module=sensitivity_matrix_heat_map$sim_network*100)
similarity_modules_plot <- sensitivity_matrix_heat_map$sim_network*100



sp_change <- levelplot(as.numeric(change_matrix$percentage_change) ~ as.numeric(change_matrix$Group.1)*as.numeric(change_matrix$Group.2),font.axis=list(fontfamily=font_family,font=font_type),main=list("Fe change %",fontfamily=font_family,font=font_type),ylab=list("diet error",fontfamily=font_family,font=font_type),xlab=list("locomotion error",fontfamily=font_family,font=font_type),col.regions = alpha(col,transparency),cuts = 1000)

similarity_ff_change_error <- levelplot(similarity_modules_plot~ error_diet_plot*error_loc_plot,font.axis=list(fontfamily=font_family,font=font_type),main=list("Functional faunas similarity",fontfamily=font_family,font=font_type),ylab=list("diet error",fontfamily=font_family,font=font_type),xlab=list("locomotion error",fontfamily=font_family,font=font_type),col.regions = alpha(col_3,transparency),cuts = 1000)

similarity_vs_sp_change <- ggplot(data, aes(sp_error, similarity_module)) +
  geom_point(color=alpha("#525252",0.3),size=4) +
  geom_smooth(method="loess", color=alpha("#774898",0.8), fill=alpha("#774898",0.3), se=TRUE) +
  theme_niwot() +
  labs(y = "FF similarity \n", x = "\n Fe %") +
  xlim(0,10) +
  ylim(0,100)+
  ggtitle("Species % vs. FF similarity \n")




quartz(width = cm(4.75),height = cm(2),family = font_family)
grid.arrange(fe_error_change,similarity_ff_change_error, ncol=2)











