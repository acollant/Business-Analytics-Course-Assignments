rm(list = ls(all.names = TRUE)) 
gc()

library(lubridate)
library(tictoc)
library(tidyverse)
library(gtrendsR)
library(markdown)
library(ggplot2)
library(rlang)
library(dplyr)
library(httr)
library(lpSolveAPI)

##=========================================================================================
##  Function to read all the rds files and load the content into a tibble data structure
##=========================================================================================
import_single_rds <- function ( file_path_, file_name_  ){
  
  full_path_ <- file.path ( file_path_ , file_name_ ) 
  
  file_ = readRDS(file = full_path_)
  
  return( file_ )
}




model_coverage_dom_set <- function(adjacency_matrix) {
  ##### Paramètres ######
  # P1
  mInter <- adjacency_matrix
  # P2
  MaxBoite <- 1
  # P3
  MinBoite <- 0
  
  ##### Variable ######
  # V1
  vX <- CVXR::Int(rows = 1, col = ncol(adjacency_matrix))
  
  ##### Objectif ######
  # O1
  objective <- CVXR::Minimize(sum(vX))
  
  ##### Contraintes ######
  constraints <- list(
    # C1
    vX >= MinBoite,
    # C2
    vX <= MaxBoite,
    # C3
    vX %*% mInter >= 1 # Chaque intersection doit être visible par au moins une boîte. 
  )
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- list(variables = vX,
                 problem = problem)
  
}  


is_sol_feasible <- function(x) {
  
  is_feasible = TRUE
  vX <- x %*% mInter
  
  for( i in 1:length(x) ){
    
    # C1 et C2
    if(x[i] < MinBoite || x[i] > MaxBoite){
      is_feasible = FALSE
      break
    }
    
    # C3
    else if(vX[i] < 1){
      is_feasible = FALSE
      break
    }
  }
  
  is_feasible
}

is_sol_x_better_than_y <- function(x,y) {
  sol_x_better = FALSE
  if (!is_sol_feasible(x)){
    sol_x_better = FALSE
  } else if (!is_sol_feasible(y)){
    sol_x_better = TRUE
  } else if (val_obj(x) < val_obj(y)) {
    sol_x_better = TRUE
  }
  sol_x_better
}

val_obj <- function(x){
  # objetif: minimiser sum(x)
  sum(x)
}

random_search <- function(domains,
                          max_iterations = 100,
                          is_sol_x_better_than_y,
                          verbose = FALSE) {
  
  best_incumbent <- rep(1, length(domains))
  best_incumbent_iteration <- 0
  curr_iteration <- 1
  # if ( verbose ) {
  #   cat('Une nouvelle solution candidate est trouvée à litération ',
  #   curr_iteration,
  #   ':',
  #   sum(best_incumbent),
  #   ' boîtes noires\n')
  # }
  candidate <- rep(1, length(domains))
  
  for (curr_iteration in 1:max_iterations) {
    
    for (i in 1:length(domains)) {
      candidate[i] <- sample(domains[[i]], 1)
    }
    
    if(is_sol_x_better_than_y(candidate, best_incumbent)) {
      
      best_incumbent <- candidate
      best_incumbent_iteration <- curr_iteration
      
      if ( verbose ) {
        
        cat('Une nouvelle solution candidate est trouvée à litération ',
            curr_iteration,
            ':',
            sum(best_incumbent),
            ' boîtes noires\n')
      }
    }
    
    curr_iteration <- curr_iteration + 1;
  }
  
  best_status <- is_sol_feasible(best_incumbent)
  
  result <- list(best_sol = sum(best_incumbent),
                 iteration = best_incumbent_iteration,
                 status = best_status)
  
  return(result)
}

# Une fonction de recherche aléatoire avec une probabilité favorisant 0 (minimise le nombre de boîte)
random_search_variant <- function(domains,
                                  max_iterations = 1000,
                                  is_sol_x_better_than_y,
                                  verbose = FALSE) {
  
  best_incumbent <- rep(1, length(domains))
  best_incumbent_iteration <- 0
  curr_iteration <- 1
  # if ( verbose ) {
  #   cat('Une nouvelle solution candidate est trouvée à litération ',
  #   curr_iteration,
  #   ':',
  #   sum(best_incumbent),
  #   ' boîtes noires\n')
  # }
  candidate <- rep(1, length(domains))
  
  for (curr_iteration in 1:max_iterations) {
    
    for (i in 1:length(domains)) {
      
      # #####################
      # Ici, si mean = 0.6, on obtient des solutions à la recherche sur large, car notre solution favorise l'obtention de 1 pour chaque case.
      # Par contre on obtient de moins bonnes solutions pour huge4 dans le prochain numéro.
      # À mean = 0.4, on répond à l'attente du prof, car on favorise un plus petit nombre de boîte et on obtient de meilleurs valeurs pour Huge4 (numéro 2.5)
      # #####################
      
      
      candidate[i] <- rnorm( 1,
                             mean = 0.4,
                             sd = 0.2) %>%
        round( digits = 0 )
    }
    
    if(is_sol_x_better_than_y(candidate, best_incumbent)) {
      
      best_incumbent <- candidate
      best_incumbent_iteration <- curr_iteration
      
      if ( verbose ) {
        cat('Une nouvelle solution candidate est trouvée à litération ',
            curr_iteration,
            ':',
            sum(best_incumbent),
            ' boîtes noires\n')
      }
    }
    
    curr_iteration <- curr_iteration + 1;
  }
  
  best_status <- is_sol_feasible(best_incumbent)
  
  result <- list(best_sol = sum(best_incumbent),
                 iteration = best_incumbent_iteration,
                 status = best_status)
  
  return(result)
}





experiment <- function(matrix_, instance_,  nb_runs_){
  
 
  
  nb_runs <- nb_runs_
  max_iterations = 100
  algo_list = c('LPSOLVE', 'aleatoire', 'variante')
  
  experiments.cols.size = 9
  
  experiments.table <- tibble(algo     = rep("",experiments.cols.size),
                              instance = rep("",experiments.cols.size),
                              no_run   = rep(0,experiments.cols.size),
                              time_sec = rep(0,experiments.cols.size),
                              obj_val  = rep(0,experiments.cols.size),
                              optimal  = rep(FALSE,experiments.cols.size),
                              realisable = rep(FALSE, experiments.cols.size))
  set.seed(45378)
  current_experiment <- 1
  mInter <- matrix_
  MaxBoite <- 1
  MinBoite <- 0
  
  for (algo in algo_list){
    for (curr_run in 1:nb_runs) {
      
      domains <- rep(list(0:1), nrow(matrix_))
      tic(quiet = TRUE) # Demarrer le timer
      
      if(algo == 'aleatoire'){
        return_value <- random_search(domains,
                                      max_iterations = max_iterations,
                                      is_sol_x_better_than_y =
                                        is_sol_x_better_than_y,
                                      verbose = FALSE)
        this_value = return_value$best_sol
        print('aleatoire')
      }
      
      if (instance_ == 'huge5') {
        opt_value = NA
        this_value = 0
      }  
      
      if(algo == 'LPSOLVE' && instance_ != 'huge5' ){
        opt_problem <- model_coverage_dom_set(matrix_)
        
        opt_result <- CVXR::psolve(opt_problem$problem,
                                   solver = 'LPSOLVE',
                                   verbose = FALSE)
        this_value = opt_result$value
        
        opt_value = NA
        if(opt_result$status == "optimal") opt_value = this_value
        print('LPSOLVE')
      }
      
      if (algo == 'variante'){
        return_value <- random_search_variant(domains,
                                              max_iterations = max_iterations,
                                              is_sol_x_better_than_y =
                                                is_sol_x_better_than_y,
                                              verbose = FALSE)
        
        this_value = return_value$best_sol
        print('variante')
      }
      
      chrono <- toc(quiet = TRUE) # arreter le timer
      opt_status <- (opt_value == this_value)
      status <- (0 < this_value && !is.na(opt_value))
      
      experiments.table$algo[current_experiment] <- algo
      experiments.table$instance[current_experiment] <- instance_#instances[rds]
      experiments.table$no_run[current_experiment] <- curr_run
      experiments.table$time_sec[current_experiment] <- (chrono$toc - chrono$tic)
      experiments.table$obj_val[current_experiment] <- this_value
      experiments.table$optimal[current_experiment] <- opt_status
      experiments.table$realisable[current_experiment] <- status
      current_experiment <- current_experiment + 1
    }
  }
  
  results <- experiments.table
  return(results)

}


   
  experiments_rds_filename = 'data/experiments.rds'
  INPUT_DATA_FILE_PATH  =  "data"
  
  files = c('Q2_huge1_interproblem_adjacency_matrix.rds',
            'Q2_huge2_interproblem_adjacency_matrix.rds',
            'Q2_huge3_interproblem_adjacency_matrix.rds',
            'Q2_huge4_interproblem_adjacency_matrix.rds',
            'Q2_huge5_interproblem_adjacency_matrix.rds')
  
  instances = c('huge1','huge2','huge3','huge4','huge5')
  current_experiment_ = c(1,10,19,28,37)
  MaxBoite <- 1
  MinBoite <- 0
  nb_runs = 3
  
  experiments.cols.size = 45 # nb_runs = 3   * nb_files = 5 * nb_algorithms = 3
  experiments.table <- tibble(algo     = rep("",experiments.cols.size),
                              instance = rep("",experiments.cols.size),
                              no_run   = rep(0,experiments.cols.size),
                              time_sec = rep(0,experiments.cols.size),
                              obj_val  = rep(0,experiments.cols.size),
                              optimal  = rep(FALSE,experiments.cols.size),
                              realisable = rep(FALSE, experiments.cols.size))
  
  
  if(file.exists(experiments_rds_filename)) {
    experiments.table <- import_single_rds(INPUT_DATA_FILE_PATH,"experiments.rds")
    
  } else{
  
      for (rds in 1:length(files)) {
        matrix_ <- import_single_rds(INPUT_DATA_FILE_PATH,files[rds])
        mInter <- matrix_
        instance_ <- instances[rds]
        e <- experiment (matrix_ , instance_, nb_runs)
        current_experiment = current_experiment_[rds]
        for ( x in 1:nrow(e)){
          experiments.table$algo[current_experiment] <- e$algo[x]
          experiments.table$instance[current_experiment] <- e$instance[x]
          experiments.table$no_run[current_experiment] <-  e$no_run[x]
          experiments.table$time_sec[current_experiment] <- e$time_sec[x]
          experiments.table$obj_val[current_experiment] <- e$obj_val[x]
          experiments.table$optimal[current_experiment] <- e$optimal[x]
          experiments.table$realisable[current_experiment] <- e$realisable[x]
          current_experiment = current_experiment + 1
        }
        
        View(experiments.table)  
      } 
  
       saveRDS(experiments.table,experiments_rds_filename)

}



# 
# 
# #
# max_iterations = 100
# nb_runs = 3
# 
# algo_list = c('LPSOLVE', 'aleatoire', 'variante')
# 
# experiments <- tibble(algo = rep("",
#                                  nb_runs * length(files) * length(algo_list)),
#                       instance = rep("",
#                                      nb_runs * length(files) * length(algo_list)),
#                       no_run = rep(0,
#                                    nb_runs * length(files) * length(algo_list)),
#                       time_sec = rep(0,
#                                      nb_runs * length(files) * length(algo_list)),
#                       obj_val = rep(0,
#                                     nb_runs * length(files) * length(algo_list)),
#                       optimal = rep(FALSE,
#                                     nb_runs * length(files) * length(algo_list)),
#                       realisable = rep(FALSE,
#                                        nb_runs * length(files) * length(algo_list)))
# 
# set.seed(45378)
# current_experiment <- 1
# 
# 
# #   # P1
# 
#   matrice <- import_single_rds(INPUT_DATA_FILE_PATH,INPUT_FILE_RDS_)
#   mInter <- matrice
#   
#   # P2
#   MaxBoite <- 1
#   # P3
#   MinBoite <- 0
#   for (algo in algo_list){
# 
#     for (curr_run in 1:nb_runs) {
# 
#       domains <- rep(list(0:1), nrow(matrice))
# 
#       tic(quiet = TRUE) # Demarrer le timer
# 
#       if(algo == 'aleatoire'){
#         return_value <- random_search(domains,
#                                       max_iterations = max_iterations,
#                                       is_sol_x_better_than_y =
#                                         is_sol_x_better_than_y,
#                                       verbose = FALSE)
#         this_value = return_value$best_sol
#         print("1")
#       }
# 
#       if(algo == 'LPSOLVE'){
#         opt_problem <- model_coverage_dom_set(matrice)
#         
#         # tryCatch(CVXR::psolve(opt_problem$problem,
#         #                       solver = 'LPSOLVE',
#         #                       verbose = FALSE), error = function(e) geterrmessage())
#         # 
#         opt_result <- CVXR::psolve(opt_problem$problem,
#                                    solver = 'LPSOLVE',
#                                    verbose = FALSE)
# 
#         this_value = opt_result$value
#         
#        
# 
#         opt_value = NA
#         if(opt_result$status == "optimal") opt_value = this_value
#         print("2")
#       }
#       if (algo == 'variante'){
#         return_value <- random_search_variant(domains,
#                                               max_iterations = max_iterations,
#                                               is_sol_x_better_than_y =
#                                                 is_sol_x_better_than_y,
#                                               verbose = FALSE)
#         this_value = return_value$best_sol
# 
#         print("3")
#       }
# 
#       chrono <- toc(quiet = TRUE) # arreter le timer
# 
#       opt_status <- (opt_value == this_value)
# 
#       status <- (0 < this_value && !is.na(opt_value))
# 
#       experiments$algo[current_experiment] <- algo
#       experiments$instance[current_experiment] <- instances[rds]
#       experiments$no_run[current_experiment] <- curr_run
#       experiments$time_sec[current_experiment] <- (chrono$toc - chrono$tic)
#       experiments$obj_val[current_experiment] <- this_value
#       experiments$optimal[current_experiment] <- opt_status
#       experiments$realisable[current_experiment] <- status
#       current_experiment <- current_experiment + 1
# 
#       # print(current_experiment)
#       # print(algo)
#     }
#   }
# 
# # #saveRDS(experiments, experiments_rds_filename)
# # saveRDS(experiments,experiments_rds_filename)
