
# script to analyze Reddit bias data

# clean environment and initialize
rm(list = ls())
gc()
set.seed(123)
options(scipen = 999)

# download necessary libraries
install.packages(c("devtools","reticulate","igraph","rstudioapi","TDA"))
devtools::install_github("shaelebrown/TDApplied")
reticulate::py_install("ripser")

# set working directory to analysis directory
path <- rstudioapi::getSourceEditorContext()$path
path <- strsplit(path,split = "analysis.R")[[1]][[1]]
setwd(path)

# verify that python can be found
reticulate::py_config()
reticulate::py_available()

# load TDApplied
library(TDApplied)

# read in four datasets
race_df <- read.csv("../data/race_df.csv")
race_df$X <- NULL
gender_df <- read.csv("../data/gender_df.csv")
gender_df$X <- NULL
orientation_df <- read.csv("../data/orientation_df.csv")
orientation_df$X <- NULL
religion_df <- read.csv("../data/religion_df.csv")
religion_df$X <- NULL
all_df <- do.call(rbind,list(race_df[sample(1:nrow(race_df),size = round(0.25*nrow(race_df))),],
                             gender_df[sample(1:nrow(gender_df),size = round(0.25*nrow(gender_df))),],
                             orientation_df[sample(1:nrow(orientation_df),size = round(0.25*nrow(orientation_df))),],
                             religion_df[sample(1:nrow(religion_df),size = round(0.25*nrow(religion_df))),]))

# compute persistence diagrams with boostrapping
ripser <- import_ripser()
race_PH <- bootstrap_persistence_thresholds(race_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
gender_PH <- bootstrap_persistence_thresholds(gender_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.3,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
orientation_PH <- bootstrap_persistence_thresholds(orientation_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
religion_PH <- bootstrap_persistence_thresholds(religion_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.2,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)
all_PH <- bootstrap_persistence_thresholds(all_df,FUN_diag = "PyH",FUN_boot = "PyH",maxdim = 2,thresh = 1.45,ripser = ripser,num_samples = 100,return_subsetted = T,return_diag = T,alpha = 0.1)

# save files
saveRDS(race_PH,"race_PH.rds")
saveRDS(gender_PH,"gender_PH.rds")
saveRDS(orientation_PH,"orientation_PH.rds")
saveRDS(religion_PH,"religion_PH.rds")
saveRDS(all_PH,"all_PH.rds")

# load files
race_PH <- readRDS("race_PH.rds")
gender_PH <- readRDS("gender_PH.rds")
orientation_PH <- readRDS("orientation_PH.rds")
religion_PH <- readRDS("religion_PH.rds")
all_PH <- readRDS("all_PH.rds")

# plot diagrams with and without thresholds
plot_diagram(D = race_PH$diag,title = "Race",max_radius = 1.3)
plot_diagram(D = race_PH$diag,title = "Race with thresholds",thresholds = race_PH,max_radius = 1.3)
plot_diagram(D = gender_PH$diag,title = "Gender",max_radius = 1.3)
plot_diagram(D = gender_PH$diag,title = "Gender with thresholds",thresholds = gender_PH,max_radius = 1.3)
plot_diagram(D = orientation_PH$diag,title = "Orientation",max_radius = 1.3)
plot_diagram(D = orientation_PH$diag,title = "Orientation with thresholds",thresholds = orientation_PH,max_radius = 1.3)
plot_diagram(D = religion_PH$diag,title = "Religion",max_radius = 1.3)
plot_diagram(D = religion_PH$diag,title = "Religion with thresholds",thresholds = religion_PH,max_radius = 1.3)
plot_diagram(D = all_PH$diag,title = "All data")
plot_diagram(D = all_PH$diag,title = "All data with thresholds",thresholds = all_PH)

# how many significant features existed in each dimension?
table(race_PH$subsetted_diag$dimension)
table(gender_PH$subsetted_diag$dimension)
table(orientation_PH$subsetted_diag$dimension)
table(religion_PH$subsetted_diag$dimension)
table(all_PH$subsetted_diag$dimension)

# one 2-sphere in race df and one loop in religion df
# most clusters in orientation df then gender df

# which diagrams are more similar to each other?
# compute embeddings in each dimension and plot
parallel_approx_distance_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,sigma = 1,rho = 1e-3,num_workers = parallelly::availableCores(omit = 1)){
  
  # create cluster
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  
  # calculate distances in parallel
  # clusters are closed if there is an error
  tryCatch(expr = {
    
    if(is.null(other_diagrams))
    {
      # not cross distance matrix, only need to compute the upper diagonal
      # since the matrix is symmetric
      d <- matrix(data = 0,nrow = length(diagrams),ncol = length(diagrams))
      u <- which(upper.tri(d),arr.ind = T)
      R <- lapply(X = 1:nrow(u),FUN = function(X){
        
        return(list(diagrams[[u[[X,1]]]],diagrams[[u[[X,2]]]]))
        
      })
      
      # remove diagrams to preserve memory
      rm(diagrams)
      
      # calculate distances in parallel, export TDApplied to nodes
      d_off_diag <- foreach::`%dopar%`(obj = foreach::foreach(r = R,.combine = c,.packages = c("TDApplied")),ex = {TDApplied::diagram_distance(D1 = r[[1]],D2 = r[[2]],dim = dim,distance = "fisher",sigma = sigma,rho = rho)})
      
      # store results in matrix
      d[upper.tri(d)] <- d_off_diag
      d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- d_off_diag
      diag(d) <- rep(0,nrow(d))
    }else
    {
      # cross distance matrix, need to compute all entries
      u <- expand.grid(1:length(other_diagrams),1:length(diagrams))
      R <- lapply(X = 1:nrow(u),FUN = function(X){
        
        return(list(other_diagrams[[u[X,1]]],diagrams[[u[X,2]]]))
        
      })
      
      # remove diagrams and other_diagrams to preserve memory
      rm(list = c("diagrams","other_diagrams"))
      
      # store distance calculations in matrix
      d[as.matrix(u)] <- foreach::`%dopar%`(foreach::foreach(r = R,.combine = cbind,.packages = c("TDApplied")),ex = {TDApplied::diagram_distance(D1 = r[[1]],D2 = r[[2]],dim = dim,distance = "fisher",sigma = sigma,rho = rho)})
      
    }
    
  }, warning = function(w){warning(w)},
  error = function(e){stop(e)},
  finally = {
    # close cluster
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
  })
  
  return(d)
  
}
diags <- list(race_PH$diag,gender_PH$diag,orientation_PH$diag,religion_PH$diag)
D0 <- parallel_approx_distance_matrix(diags,dim = 0,sigma = 0.01)
D1 <- parallel_approx_distance_matrix(diags,dim = 1,sigma = 0.01)
D2 <- parallel_approx_distance_matrix(diags,dim = 2,sigma = 0.01)
emb0 <- diagram_mds(diags,D = D0,dim = 0,k = 2,distance = "fisher",sigma = 0.01,rho = 1e-3)
emb1 <- diagram_mds(diags,D = D1,dim = 1,k = 2,distance = "fisher",sigma = 0.01,rho = 1e-3)
emb2 <- diagram_mds(diags,D = D2,dim = 2,k = 2,distance = "fisher",sigma = 0.01,rho = 1e-3)
labs <- c("Race","Gender","Orientation","Religion")
plot(emb0[,1],emb0[,2],xlim = c(-0.4,0.4),ylim = c(-0.4,0.4),xlab = "Embedding dim 1",ylab = "Embedding dim 2",main = "H0",pch = 19)
text(emb0[,1],emb0[,2],labels = labs,pos = 3)
plot(emb1[,1],emb1[,2],xlim = c(-0.2,0.2),ylim = c(-0.2,0.2),xlab = "Embedding dim 1",ylab = "Embedding dim 2",main = "H1",pch = 19)
text(emb1[,1],emb1[,2],labels = labs,pos = 3)
plot(emb2[,1],emb2[,2],xlim = c(-0.15,0.15),ylim = c(-0.15,0.15),xlab = "Embedding dim 1",ylab = "Embedding dim 2",main = "H2",pch = 19)
text(emb2[,1],emb2[,2],labels = labs,pos = 3)

# highly fragmented spaces except for the religion space
# let's plot a rips graph to see what's going on..
birth <- religion_PH$subsetted_diag$birth[[nrow(religion_PH$subsetted_diag)]]
death <- religion_PH$subsetted_diag$death[[nrow(religion_PH$subsetted_diag)]]
diag_with_reps <- TDA::ripsDiag(as.matrix(dist(religion_df)),maxdimension = 1,dist = "arbitrary",maxscale = 1.2,location = T,library = "dionysus")
rep_ind <- 553
rep <- unique(as.vector(diag_with_reps$cycleLocation[[rep_ind]]))
g <- rips_graphs(race_df,eps = c(birth))
cols <- rep("blue",nrow(religion_df))
cols[rep] <- "red"
plot_rips_graph(g,eps = birth,plot_isolated_vertices = F,vertex_labels = F,component_of = rep[[1]],cols = cols)




