## Read in swiss species data and environmental layers. Use switzerland as the study domain
## Read in fitted model
## Overlay butterflies with their interacting plant genera
##Jussi Makinen 7-11-2024
## Elena Quintero 7-9-2026 - path changes and matrix rows/cols refer for species

Overlay_analysis = function(papilio_test, model) {
  
  if (papilio_test) {
    output_dir = '/scratch/project_2011416/Output/Papilio'
  } else {
    output_dir = '/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/range_overlap/output'
  }
  
  #read in mean predicted probability
  if (model == 'spat') {
    pred_current_realized = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_current_non_spat_no_hpc.csv')
  }
  pred_current = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_current_non_spat_no_hpc.csv')
  pred_future_26 = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_future_26_non_spat_no_hpc.csv')
  pred_future_85 = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_future_85_non_spat_no_hpc.csv')
  
  #read in interaction matrix
  if (papilio_test) {
    species_matrix = read.csv('/scratch/project_2011416/Data/papilio/interaction_matrix_papilio.csv')
  } else {
    species_matrix = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/interaction/data/Swiss_web_gen_strict.csv') |> 
      column_to_rownames("X") |> t() |> as.data.frame() |> rownames_to_column("X")
  }
  
  #compute overlap between butterflies and their host plants
  #current realized realized
  if (model == 'spat') {
    range_butterfly = {}
    range_hosts = {}
    overlap_host = {}
    overlap_all = {}
    plant_names = colnames(species_matrix)
    butterfly_names = species_matrix$X
    butterfly_names = butterfly_names[butterfly_names %in% colnames(pred_current_realized)]
    
    for (i in 1:length(butterfly_names)) {
      host_names = plant_names[species_matrix[i,]==1]
      pred_current_host = as.matrix(pred_current_realized[,colnames(pred_current_realized) %in% host_names])
      pred_current_all = as.matrix(pred_current_realized[,colnames(pred_current_realized) %in% plant_names])
      if (ncol(pred_current_host) >= 2) {
        pred_current_host = apply(pred_current_host, 1, function(x) 1-prod(1-x))
      }
      pred_current_all = round(apply(pred_current_all, 1, function(x) 1-prod(1-x)),5)
      
      range_butterfly[i] = sum(pred_current_realized[,colnames(pred_current_realized)==butterfly_names[i]])
      range_hosts[i] = sum(pred_current_host)
      overlap_host[i] = sum(pred_current_realized[,colnames(pred_current_realized)==butterfly_names[i]] * pred_current_host)
      overlap_all[i] = sum(pred_current_realized[,colnames(pred_current_realized)==butterfly_names[i]] * pred_current_all)
    }
    
    #collect results
    range_current_df = data.frame(Species = butterfly_names,
                                  Range = range_butterfly,
                                  Range_hosts = range_hosts,
                                  Overlap_hosts = overlap_host,
                                  Overlap_any_plant = overlap_all)
    
    
    write.csv(range_current_df, file = paste0(output_dir, '/range_overlap_current_realized_', model, '.csv'), quote = F)
  }
  
  
  #current
  range_butterfly = {}
  range_hosts = {}
  overlap_host = {}
  overlap_all = {}
  plant_names = colnames(species_matrix)
  butterfly_names = species_matrix$X
  butterfly_names = butterfly_names[butterfly_names %in% colnames(pred_current)]
  
  for (i in 1:length(butterfly_names)) {
    # host_names = plant_names[species_matrix[i,]==1]
    # changing the way of getting host_names for the butterflies because now the order and no. of species in the 
    #  interaction matrix is not the same as in the pred_current data frame
    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 
    pred_current_host = as.matrix(pred_current[,colnames(pred_current) %in% host_names])
    pred_current_all = as.matrix(pred_current[,colnames(pred_current) %in% plant_names])
    if (ncol(pred_current_host) >= 2) {
      pred_current_host = apply(pred_current_host, 1, function(x) 1-prod(1-x))
    }
    pred_current_all = round(apply(pred_current_all, 1, function(x) 1-prod(1-x)),5)
    
    range_butterfly[i] = sum(pred_current[,colnames(pred_current)==butterfly_names[i]])
    range_hosts[i] = sum(pred_current_host)
    overlap_host[i] = sum(pred_current[,colnames(pred_current)==butterfly_names[i]] * pred_current_host)
    overlap_all[i] = sum(pred_current[,colnames(pred_current)==butterfly_names[i]] * pred_current_all)
  }
  
  #collect results
  range_current_df = data.frame(Species = butterfly_names,
                                Range = range_butterfly,
                                Range_hosts = range_hosts,
                                Overlap_hosts = overlap_host,
                                Overlap_any_plant = overlap_all)
  
  
  write.csv(range_current_df, file = paste0(output_dir, '/range_overlap_current_', model, '.csv'), quote = F)
  
  
  #future 26
  range_butterfly = {}
  range_hosts = {}
  overlap_host = {}
  overlap_all = {}
  
  for (i in 1:length(butterfly_names)) {
    # host_names = plant_names[species_matrix[i,]==1]
    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 
    pred_future_host = as.matrix(pred_future_26[,colnames(pred_future_26) %in% host_names])
    pred_future_all = as.matrix(pred_future_26[,colnames(pred_future_26) %in% plant_names])
    if (ncol(pred_future_host) >= 2) {
      pred_future_host = apply(pred_future_host, 1, function(x) 1-prod(1-x))
    }
    pred_future_all = apply(pred_future_all, 1, function(x) 1-prod(1-x))
    
    range_butterfly[i] = sum(pred_future_26[,colnames(pred_future_26)==butterfly_names[i]])
    range_hosts[i] = sum(pred_future_host)
    overlap_host[i] = sum(pred_future_26[,colnames(pred_future_26)==butterfly_names[i]] * pred_future_host)
    overlap_all[i] = sum(pred_future_26[,colnames(pred_future_26)==butterfly_names[i]] * pred_future_all)
  }
  
  #collect results
  range_future_df = data.frame(Species = butterfly_names,
                                Range = range_butterfly,
                                Range_hosts = range_hosts,
                                Overlap_hosts = overlap_host,
                                Overlap_any_plant = overlap_all)
  
  write.csv(range_future_df, file = paste0(output_dir, '/range_overlap_future_26_', model, '.csv'), quote = F)
  
  
  #future 85
  range_butterfly = {}
  range_hosts = {}
  overlap_host = {}
  overlap_all = {}
  
  for (i in 1:length(butterfly_names)) {
    # host_names = plant_names[species_matrix[i,]==1]
    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 
    pred_future_host = as.matrix(pred_future_85[,colnames(pred_future_85) %in% host_names])
    pred_future_all = as.matrix(pred_future_85[,colnames(pred_future_85) %in% plant_names])
    if (ncol(pred_future_host) >= 2) {
      pred_future_host = apply(pred_future_host, 1, function(x) 1-prod(1-x))
    }
    pred_future_all = apply(pred_future_all, 1, function(x) 1-prod(1-x))
    
    range_butterfly[i] = sum(pred_future_85[,colnames(pred_future_85)==butterfly_names[i]])
    range_hosts[i] = sum(pred_future_host)
    overlap_host[i] = sum(pred_future_85[,colnames(pred_future_85)==butterfly_names[i]] * pred_future_host)
    overlap_all[i] = sum(pred_future_85[,colnames(pred_future_85)==butterfly_names[i]] * pred_future_all)
  }
  
  #collect results
  range_future_df = data.frame(Species = butterfly_names,
                               Range = range_butterfly,
                               Range_hosts = range_hosts,
                               Overlap_hosts = overlap_host,
                               Overlap_any_plant = overlap_all)
  
  write.csv(range_future_df, file = paste0(output_dir, '/range_overlap_future_85_', model, '.csv'), quote = F)
  
}

Overlay_analysis(papilio_test = F, model = "non_spat")
