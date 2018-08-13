### Code for Section 7.2
###################################################
# You neeed a folder called "Air pollution in Piemonte" inside your working directory (my.dir)
# with the data downloaded from
# https://sites.google.com/a/r-inla.org/stbook/datasets
remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

require(INLA)
inla.setOption(scale.model.default=FALSE)

require(splancs)
require(sp)
require(fields)
require(maptools)
require(lattice)
require(abind)

# Load the data for the 24 stations and 182 days
Piemonte_data <- read.table(paste(my.dir,"Air pollution in Piemonte/Piemonte_data_byday.csv",sep=""),header=TRUE,sep=",")
head(Piemonte_data)

View(Piemonte_data)

coordinates <-read.table(paste(my.dir,"Air pollution in Piemonte/coordinates.csv",sep=""),header=TRUE,sep=",")
rownames(coordinates) <- coordinates[,"Station.ID"]
# Borders of Piemonte (in km)
borders <-read.table(paste(my.dir,"Air pollution in Piemonte/Piemonte_borders.csv",sep=""),header=TRUE,sep=",")
View(borders)

n_stations <- length(coordinates$Station.ID) #24 stations
n_data <- length(Piemonte_data$Station.ID) #4368 space-time data
n_days <- n_data/n_stations #182 time points

Piemonte_data$time <- rep(1:n_days, each=n_stations)
coordinates.allyear <- as.matrix(coordinates[Piemonte_data$Station.ID, c("UTMX","UTMY")])
dim(coordinates.allyear)

Piemonte_data$logPM10 <- log(Piemonte_data$PM10)
mean_covariates <- apply(Piemonte_data[,3:10],2,mean)
sd_covariates <- apply(Piemonte_data[,3:10],2,sd)
Piemonte_data[,3:10] <- scale(Piemonte_data[,3:10],center=mean_covariates, scale=sd_covariates)

# Load the covariate arrays (each array except for A is 56x72x182)
load(paste(my.dir,"Air pollution in Piemonte/Covariates/Altitude_GRID.Rdata",sep="")) #A; AltitudeGRID
load(paste(my.dir,"Air pollution in Piemonte/Covariates/WindSpeed_GRID.Rdata",sep="")) #WS; WindSpeedGRID
load(paste(my.dir,"Air pollution in Piemonte/Covariates/HMix_GRID.Rdata",sep="")) #HMIX; HMixMaxGRID
load(paste(my.dir,"Air pollution in Piemonte/Covariates/Emi_GRID.Rdata",sep="")) #EMI; EmiGRID
load(paste(my.dir,"Air pollution in Piemonte/Covariates/Temp_GRID.Rdata",sep="")) #TEMP; Mean_Temp
load(paste(my.dir,"Air pollution in Piemonte/Covariates/Prec_GRID.Rdata",sep="")) #PREC; Prec
# Load the Piemonte grid c(309,529),c(4875,5159),dims=c(56,72)
load(paste(my.dir,"Air pollution in Piemonte/Covariates/Piemonte_grid.Rdata",sep=""))


# Extract the standardized covariate for day i_day (you get a 56X72X8 matrix)
i_day <- 122
which_date <- unique(Piemonte_data$Date)[i_day]
print(paste("**---- You will get a prediction for ", which_date, "---**"))

# Standardise the covariates for the selected day
source(paste(my.dir,"Air pollution in Piemonte/Covariates/covariates_selector.R",sep=""))

covariate_array_std <- covariates_selector_funct(i_day, mean_covariates, sd_covariates)

# Set to NA the (standardized) altitude values >7 (1000 n)
elevation <- covariate_array_std[,,1]
index_mountains <- which(elevation>7)
elevation[elevation>7] <- NA
covariate_array_std[,,1] <- elevation

# Reshape the 3D array (56x72x8) into a dataframe (4032x8) with the 8 covariates on the columns
covariate_matrix_std <- data.frame(apply(covariate_array_std,3,function(X) c(t(X))))
colnames(covariate_matrix_std) <- colnames(Piemonte_data[,3:10])


# *** Code for Figure 7.8 top
plot(Piemonte_grid,col="grey",pch=18, asp=1, xlim=range(Piemonte_grid$UTMX_km))
lines(borders, lwd=3, asp=1)
points(coordinates$UTMX, coordinates$UTMY, pch=20, cex=2)
# ***
plot(borders$UTM_X, borders$UTM_Y)

# *** Code for Figure 7.8 bottom
Piemonte_mesh <- inla.mesh.2d(loc=cbind(coordinates$UTMX,coordinates$UTMY), 
                              loc.domain=borders, offset=c(10, 140), max.edge=c(50, 1000))
plot(Piemonte_mesh,asp=1,main="")
lines(borders, lwd=3)
points(coordinates$UTMX, coordinates$UTMY, pch=20, cex=2)
# ***


# *** Toronto
load("CD_Toronto_X.RDA")
load("TO_Border.RDA")
coordinates <- CD_Toronto_X[c(3,4,8,11,12,13),]
borders <- TO_border
save(x = Toronto_mesh, file = "Toronto_mesh.RDA")
View(Toronto_mesh)
Toronto_mesh <- inla.mesh.2d(loc=cbind(coordinates$x,coordinates$y),
                             loc.domain=borders, max.edge= 1000) #, offset=c(10, 140) 

plot(Toronto_mesh, asp=1, main = "Toronto Mesh 1")
points(borders, lwd=1)
points(coordinates$x, coordinates$y, pch=20, cex=2, col = "blue")
# points(coordinates$x[c(3,4,8,11,12,13)], coordinates$y[c(3,4,8,11,12,13)], pch=20, cex=2, col = "red")

?inla.mesh.2d
# ***

Piemonte_spde <- inla.spde2.matern(mesh = Piemonte_mesh, alpha=2)
A_est <- inla.spde.make.A(mesh = Piemonte_mesh,
                          loc = coordinates.allyear,
                          group = Piemonte_data$time,
                          n.group = n_days)
# *** Toronto
Toronto_spde <- inla.spde2.matern(mesh = Toronto_mesh, alpha=2)
A_est <- inla.spde.make.A(mesh = Toronto_mesh,
                          loc = coordinates.allyear,
                          group = Piemonte_data$time,
                          n.group = n_days)
dim(A_est)

s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = Piemonte_spde$n.spde,
                                n.group = n_days)
names(s_index)

stack_est <- inla.stack(data = list(logPM10 = Piemonte_data$logPM10),
                        A = list(A_est, 1),
                        effects = list(c(s_index,list(Intercept = 1)), list(Piemonte_data[,3:10])), tag = "est")

A_pred <- inla.spde.make.A(mesh=Piemonte_mesh,
                           loc=as.matrix(Piemonte_grid),
                           group=i_day,  #selected day for prediction
                           n.group=n_days)
stack_pred <- inla.stack(data=list(logPM10=NA),
                         A=list(A_pred,1),
                         effects=list(c(s_index,list(Intercept=1)), list(covariate_matrix_std)),
                         tag="pred")

stack <- inla.stack(stack_est, stack_pred)

inla(y ~ 1, data = data.frame(y = 1:10))
formula <- logPM10 ~ -1 + Intercept + A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI + 
  f(spatial.field, model=Piemonte_spde,group=spatial.field.group, control.group=list(model="ar1"))

# ATTENTION: the run is computationally intensive!
# Code crashes here
output <- inla(formula,
               data = inla.stack.data(stack, spde = Piemonte_spde),
               family ="gaussian", verbose=TRUE,
               control.predictor=list(A=inla.stack.A(stack), compute = TRUE))   

# Fixed effects betas
fixed.out <- round(output$summary.fixed,3)
# Hyperparameters sigma2eps and AR(1) a
rownames(output$summary.hyperpar)

sigma2e_marg <- inla.tmarginal(function(x) 1/x,output$marginals.hyperpar[[1]])
sigma2e_m1 <- inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 <- inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev <- sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)

ar <- output$summary.hyperpar["GroupRho for spatial.field",]

# Spatial parameters sigma2 and range
mod.field <- inla.spde2.result(output, name="spatial.field", Piemonte_spde)

var.nom.marg <- mod.field$marginals.variance.nominal[[1]]
var.nom.m1 <- inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 <- inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev <- sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

range.nom.marg <- mod.field$marginals.range.nominal[[1]]
range.nom.m1 <- inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 <- inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev <- sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

index_pred <- inla.stack.index(stack,"pred")$data
lp_marginals <- output$marginals.linear.predictor[index_pred]

lp_mean <- unlist(lapply(lp_marginals, function(x) inla.emarginal(exp, x)))
lp_grid_mean <- matrix(lp_mean, 56, 72, byrow=T)

# Select only points inside Piemonte and set NA to the outer points 
lp_grid_mean[index_mountains] <- NA
library(splancs)
inside_Piemonte <- matrix(inout(Piemonte_grid, borders), 56, 72, byrow=T)
inside_Piemonte[inside_Piemonte==0] <- NA
inside_lp_grid_mean <- inside_Piemonte *  lp_grid_mean

seq.x.grid <- seq(range(Piemonte_grid[,1])[1],range(Piemonte_grid[,1])[2],length=56)
seq.y.grid <- seq(range(Piemonte_grid[,2])[1],range(Piemonte_grid[,2])[2],length=72)

# *** Code for Figure 7.9
print(levelplot(x=inside_lp_grid_mean,
                row.values=seq.x.grid,
                column.values=seq.y.grid,
                ylim=c(4875,5159), xlim=c(309,529),
                col.regions=gray(seq(.9,.2,l=100)),
                aspect="iso",
                contour=TRUE, labels=FALSE, pretty=TRUE, 
                xlab="",ylab=""))
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(borders,col=1,cex=.25)
lpoints(coordinates$UTMX, coordinates$UTMY,col=1,lwd=2,pch=21)
trellis.unfocus()
# ***

# *** Code for Figure 7.10
threshold <- log(50)
prob  <- lapply(X=lp_marginals, FUN=function(x) inla.pmarginal(marginal=x,threshold))
tailprob_grid <- matrix(1-unlist(prob),56,72, byrow=T)

tailprob_grid[index_mountains] <- NA
inside_tailprob_grid <- inside_Piemonte *  tailprob_grid

print(levelplot(x=inside_tailprob_grid,
                row.values=seq.x.grid,
                column.values=seq.y.grid,
                ylim=c(4875,5159), xlim=c(309,529),
                at=seq(0,1,by=.1),
                col.regions=gray(seq(.9,.2,l=100)),
                aspect="iso",
                contour=TRUE, labels=FALSE, pretty=TRUE, 
                xlab="",ylab=""))
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(borders,col=1,cex=.25)
lpoints(coordinates$UTMX, coordinates$UTMY,col=1,lwd=2,pch=21)
trellis.unfocus()