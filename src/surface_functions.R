# Functions for sampling surface

specify_corr <- function(surface_data){
  distances <- as.matrix(dist(surface_data))
  corr <- exp(-distances/3)
  return(corr)
}

distance_func <- function(x, center_coord){
  dist <- sqrt((x[1] - center_coord[1])^2 + (x[2] - center_coord[2])^2)
  return(dist)
}

get_sampling_surface_simple <- function(k){
  # Define the window of interest
  dim <- c(k, k)
  win <- owin(c(0,dim[1]), c(0,dim[2]))
  
  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  multiplier <- 1/dim[2]
  # Get coordinate combinations to make a dataframe,
  # in the same format as other data gen functions
  surface_data <- expand.grid(x0, y0)
  names(surface_data) <- c("x", "y")
  
  aux_x <- outer(x0,y0, function (x,y) multiplier*y + 0*x)
  aux_z <- outer(x0,y0, function (x,y) 0*y + multiplier*x)
  #plot(im(aux_x))
  #plot(im(aux_z))
  surface_data$aux_x <- c(t(aux_x))
  surface_data$aux_z <- c(t(aux_z))
  
  # Intensity and bias
  p <- ggplot(surface_data, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, name="X") +
    ggtitle("Initial sampling surface, X covariate") +
    #theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  print(p)
  
  p <- ggplot(surface_data, aes(x, y, fill=aux_z)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, "Z") +
    ggtitle("Initial sampling surface, Z covariate") +
    #theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  print(p)
  
  return(surface_data)
}

get_sampling_surface_donut <- function(k){
  # Function generates X and correlated Z auxiliary data
  x <- seq(1:k)
  y <- seq(1:k)
  center_coord <- c(k/2, k/2)
  # Generate all possible coordinate points
  sampling_grid <- expand.grid(x, y)
  sampling_grid
  # Then compute distance from each location to center of the grid.
  r <- apply(sampling_grid, 1, distance_func, center_coord=center_coord)
  r
  # Then create x covariate
  x <- exp(-10*((r-7)/5)^2)
  # This gives us the data for the auxiliary surveys x_i in the paper, which will be the 
  # covariates used in the paper X_i^T = [1, X_i1]
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_x <- x_1
  # This gives us auxiliary data for the intensity parameter lambda.
  # Now generate the plot of auxiliary data for the grid
  colnames(sampling_grid) <- c("x", "y", "aux_x")
  # Heatmap 
  p <- ggplot(sampling_grid, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Initial sampling surface, X covariate")
  print(p)
  return(sampling_grid)
}

get_bias_surface_correlated <- function(sampling_grid, aux_cor){
  # Various ways to generate correlated vector from one already existing:
  # https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  V <- matrix(c(1, aux_cor, aux_cor, 1), nrow=2, ncol=2)
  R <- chol(V)
  aux_z <- rnorm(k*k, 0, 1)
  X <- cbind(sampling_grid$aux_x, aux_z)
  cor_X <- X %*% R 
  print("Correlation between X and Z")
  print(cor(cor_X))
  #cor_X[1:5, 1]
  sampling_grid$aux_z <- cor_X[,2]
  p <- ggplot(sampling_grid, aes(x, y, fill=aux_z)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Initial sampling surface, Z covariate")
  print(p)
  return(sampling_grid)
}

get_bias_surface_exponential <- function(sampling_grid, centroid){
  x <- seq(1:k)
  y <- seq(1:k)
  # Generate all possible coordinate points
  grid_points <- expand.grid(x, y)
  # Then compute distance from centroid to every other location
  r <- apply(grid_points, 1, function(
    x, center_coord){sqrt((center_coord[1] - x[1])^2 + (center_coord[2] - x[2])^2)}, center_coord=centroid)
  # Then create z covariate
  x <- exp(-2*((r)/5))
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_z <- x_1
  p <- ggplot(sampling_grid, aes(x, y, fill=aux_z)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Initial sampling surface, Z covariate")
  print(p)
  return(sampling_grid)
}