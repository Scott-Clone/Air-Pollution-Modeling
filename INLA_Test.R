# CHeck this code

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
