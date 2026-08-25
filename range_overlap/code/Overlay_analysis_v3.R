## Read in swiss species data and environmental layers. Use switzerland as the study domain
## Read in fitted model
## Overlay butterflies with their interacting plant genera
##Jussi Makinen 7-11-2024
## small add-on Elena Quintero 14-08-2026

Overlay_analysis = function(papilio_test, model) {
  
  if (papilio_test) {
    output_dir = '/scratch/project_2011416/Output/Papilio'
  } else {
    output_dir = '/scratch/project_2011416/Output/All_species'
  }
  
  #read in mean predicted probability
  if (model == 'spat') {
    pred_current_realized = read.csv(paste0(output_dir, '/pred_all_mean_current_', model, '.csv'))
  }
  pred_current = read.csv(paste0(output_dir, '/pred_all_mean_current_', model, '.csv'))
  pred_future_26 = read.csv(paste0(output_dir, '/pred_all_mean_future_26_', model, '.csv'))
  pred_future_85 = read.csv(paste0(output_dir, '/pred_all_mean_future_85_', model, '.csv'))
  
  #read in interaction matrix
  if (papilio_test) {
    species_matrix = read.csv('/scratch/project_2011416/Data/papilio/interaction_matrix_papilio.csv')
  } else {
    species_matrix = read.csv('/scratch/project_2011416/Data/interaction_matrix.csv')
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
    host_names = plant_names[species_matrix[i,]==1]
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
    host_names = plant_names[species_matrix[i,]==1]
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
    host_names = plant_names[species_matrix[i,]==1]
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
