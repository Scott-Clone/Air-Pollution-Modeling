
#library(rlist)

# Make data frame to store results
NO2_2005_day1_30 <- as.data.frame(matrix(nrow = 6364, ncol = 30))
colnames(NO2_2005_day1_30) <- 1:30
View(NO2_2005_day1_30)
for(i_day in 1:30){
  i_day <- i_day
  #load("CD_Toronto_X.RDA")
  #load("TO_Border.RDA")
  #load("Toronto_mesh.RDA")
  #load("NO2_2005.RDA")
  load("Toronto_Data.RDA")
  #save(x = cov_matrix_std, file = "cov_matrix_std.RDA")
  #load("cov_matrix_std.RDA")
  #load("TO_Out.RDA")
  
  Fix_TO_Data <- Toronto_Data[6:7]
  Fix_TO_Data[] <- lapply(Fix_TO_Data, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  Toronto_Data$WS <- Fix_TO_Data$WS
  Toronto_Data$WD <- Fix_TO_Data$WD
  
  # Data for calc
  Toronto_Data <- Toronto_Data[1:180,]
  
  # Prepare Borders for Toronto
  TO_borders <- TO_border 
  
  
  # Prepare coordinates for stations
  TO_coords <- CD_Toronto_X[c(3,4,8,11,12,13),1:3]
  colnames(TO_coords) <- c("Station.ID", "UTMX", "UTMY")
  
  # Prepare Data for Pollutant
  n_stations <- length(TO_coords[,1]) #number of stations NB data sparced
  n_data <- length(Toronto_Data[,1]) #number of space-time data
  n_days <- n_data/n_stations #number of time points
  
  #rownames(TO_coords) <- TO_coords[,"Station.ID"]
  rownames(TO_coords) <- 1:6
  TO_coords$Station.ID <- 1:6
  Toronto_Data$Station.ID <- rep(1:6, n_days) 
  
  Toronto_Data$time <- rep(1:n_days, each = n_stations)
  coords.allyear <- as.matrix(TO_coords[Toronto_Data$Station.ID, c("UTMX","UTMY")])
  dim(coords.allyear)
  
  # Standardize covariates and calculate log pollutant
  Toronto_Data$logNO2 <- log(Toronto_Data$NO2)
  mean_covariates <- apply(Toronto_Data[, 3:9], 2, mean, na.rm = TRUE)
  sd_covariates <- apply(Toronto_Data[, 3:9], 2, sd, na.rm = TRUE)
  Toronto_Data[, 3:9] <- scale(Toronto_Data[, 3:9],center = mean_covariates, 
                               scale = sd_covariates)
  # Need to standardize Altitude separate as there isno sd (all stations are same elev)
  Toronto_Data$A <- 173.4/100
 
  # Extract the standardized covariate for day i_day (you get a 56X72X8 matrix)
  
  which_date <- unique(Toronto_Data$Date)[i_day]
  print(paste("**---- You will get a prediction for ", which_date, "---**"))
  
  # Toronto_mesh <- inla.mesh.2d(loc=cbind(TO_coords$UTMX, TO_coords$UTMY),
  #                             loc.domain = TO_border, max.edge= c(2000, 3500), 
  #                             offset = c(400, 5000))  
  
  Toronto_mesh <- inla.mesh.2d(loc=cbind(TO_coords$UTMX, TO_coords$UTMY),
                               loc.domain = TO_border, max.edge= c(4000, 10000), 
                               offset = c(300, 6000))  
  # Make Toronto_grid <- grid
  #load("TO_map.RDA")
  #dd <- subset(TO_map, TO_map$FCODE_DESC == "Major Shoreline" | TO_map$FCODE_DESC == "Geostatistical line")
  
  #Toronto_grid1 <- makegrid(dd, cellsize = 500, offset = c(-0.5, -4), pretty = TRUE)
  #length(Toronto_grid1$x1)
  # 74* 86
  Toronto_grid <- makegrid(TO_Out, cellsize = 500, offset = c(-0.5, -4), pretty = TRUE)
  
  Toronto_grid <- as.data.frame(cbind(Toronto_grid$x1, Toronto_grid$x2))
  colnames(Toronto_grid) <- c("UTMX_km", "UTMY_km")
  
  plot(Toronto_grid, col="grey",pch = 18, asp = 1, xlim = range(Toronto_grid$UTMX_km)) #, 
  points(TO_border, lwd = 1, pch=20, cex=1)
  points(TO_coords$UTMX, TO_coords$UTMY, pch=20, cex=2, col = "blue")
  
  Toronto_spde <- inla.spde2.matern(mesh = Toronto_mesh, alpha=2)
  A_est <- inla.spde.make.A(mesh = Toronto_mesh,
                            loc = coords.allyear,
                            group = Toronto_Data$time,
                            n.group = n_days)
  # *** Toronto
  Toronto_spde <- inla.spde2.matern(mesh = Toronto_mesh, alpha=2)
  A_est <- inla.spde.make.A(mesh = Toronto_mesh,
                            loc = coords.allyear,
                            group = Toronto_Data$time,
                            n.group = n_days)
  dim(A_est)
  
  s_index <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = Toronto_spde$n.spde,
                                  n.group = n_days)
  names(s_index)
  
  stack_est <- inla.stack(data = list(logNO2 = Toronto_Data$logNO2),
                          A = list(A_est, 1),
                          effects = list(c(s_index, list(Intercept = 1)), 
                                         list(Toronto_Data[,3:9])), tag = "est")
  
  
  A_pred <- inla.spde.make.A(mesh = Toronto_mesh,
                             loc = as.matrix(Toronto_grid), #?????
                             group = i_day,  #selected day for prediction
                             n.group = n_days)
  #Toronto_Data[i_day,]
  #View(Toronto_Data[3:9])
  
  # cov_matrix_std <- as.data.frame(cbind(rep(Toronto_Data[i_day,3], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,4], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,5], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,6], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,7], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,8], length(Toronto_grid$x1)),
  #                                      rep(Toronto_Data[i_day,9], length(Toronto_grid$x1))))
  
  cov_matrix_std <- as.data.frame(cbind(rep(Toronto_Data[i_day,3], length(Toronto_grid[,1])),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,4], Toronto_Data[i_day,4]+0.1),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,5], Toronto_Data[i_day,5]+0.4),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,6], Toronto_Data[i_day,6]+0.6),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,7], Toronto_Data[i_day,7]+0.8),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,8], Toronto_Data[i_day,8]+0.2),
                                        runif(length(Toronto_grid[,1]), Toronto_Data[i_day,9], Toronto_Data[i_day,9]+0.8)))
  
  colnames(cov_matrix_std) <- colnames(Toronto_Data[,3:9])
  
  stack_pred <- inla.stack(data = list(logNO2 = NA),
                           A = list(A_pred, 1),
                           effects = list(c(s_index, list(Intercept = 1)), 
                                          list(cov_matrix_std)), tag="pred") #?????
  
  stack <- inla.stack(stack_est, stack_pred)
  
  inla(y ~ 1, data = data.frame(y = 1:9))
  formula <- logNO2 ~ -1 + Intercept + UTMX + UTMY + WS + A + WD + PREC + TEMP + #A+ 
    f(spatial.field, model = Toronto_spde, group = spatial.field.group, 
      control.group = list(model = "ar1"))
  
  # ATTENTION: the run is computationally intensive!
  output <- inla(formula, data = inla.stack.data(stack, spde = Toronto_spde),
                 family = "gaussian", verbose = T,
                 control.predictor = list(A = inla.stack.A(stack), compute = TRUE))  
  
  # After running model check length cov matrix for plotting window
  #length(cov_matrix_std$A)
  
  # Fixed effects betas
  fixed.out <- round(output$summary.fixed, 3)
  
  # Hyperparameters sigma2eps and AR(1) a
  rownames(output$summary.hyperpar)
  
  sigma2e_marg <- inla.tmarginal(function(x) 1/x,output$marginals.hyperpar[[1]])
  sigma2e_m1 <- inla.emarginal(function(x) x, sigma2e_marg)
  sigma2e_m2 <- inla.emarginal(function(x) x^2, sigma2e_marg)
  sigma2e_stdev <- sqrt(sigma2e_m2 - sigma2e_m1^2)
  sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
  
  ar <- output$summary.hyperpar["GroupRho for spatial.field",]
  
  # Spatial parameters sigma2 and range
  mod.field <- inla.spde2.result(output, name="spatial.field", Toronto_spde)
  
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
  
  #lp_grid_mean <- matrix(lp_mean, 86, 74, byrow = F)
  NO2_2005_day1_30$i_day <- lp_mean
}

