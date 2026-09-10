## Read in swiss species data and environmental layers. Use switzerland as the study domain
## Read in fitted model
## Overlay butterflies with their interacting plant genera
## Jussi Makinen 7-11-2024
## Elena Quintero 7-9-2026 - path changes and matrix rows/cols refer for species
## Elena - update to calculate range for each host (not in general)
library(dplyr)
library(tibble)

  output_dir = '/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/range_overlap/output'
  
  #read in mean predicted probability
  pred_current = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_current_non_spat_no_hpc.csv')
  pred_future_26 = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_future_26_non_spat_no_hpc.csv')
  pred_future_85 = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/hmsc/output/pred_all_mean_future_85_non_spat_no_hpc.csv')
  
  #read in interaction matrix
  species_matrix = read.csv('/Users/elqu2194/Documents/SLU_EvoNets/swiss_butterflies_hmsc/interaction/data/Swiss_web_gen_strict.csv') |> 
      column_to_rownames("X") |> t() |> as.data.frame() |> rownames_to_column("X")
  
  #current
  range_current_df = data.frame()
  
  plant_names = colnames(species_matrix)
  butterfly_names = species_matrix$X
  butterfly_names = butterfly_names[butterfly_names %in% colnames(pred_current)]

  for (i in 1:length(butterfly_names)) {

    df = data.frame()

    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 

    for (j in 1:length(host_names)) {

      df = rbind(df, cbind(
            Species = butterfly_names[i],
            Host = host_names[j],
            Species_range = sum(pred_current[,colnames(pred_current)==butterfly_names[i]]),
            Host_range = sum(pred_current[,colnames(pred_current)==host_names[j]]),
            Overlap_host = sum(pred_current[,colnames(pred_current)==butterfly_names[i]] * 
                                pred_current[,colnames(pred_current)==host_names[j]])))
    }

    range_current_df = rbind(range_current_df, df)
  }
  
  write.csv(range_current_df, file = paste0(output_dir, '/range_overlap_hosts_current.csv'), quote = F)
  
  
  #future 26
  
  range_fut26_df = data.frame()
  
  plant_names = colnames(species_matrix)
  butterfly_names = species_matrix$X
  butterfly_names = butterfly_names[butterfly_names %in% colnames(pred_current)]

  for (i in 1:length(butterfly_names)) {

    df = data.frame()

    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 

    for (j in 1:length(host_names)) {

      df = rbind(df, cbind(
            Species = butterfly_names[i],
            Host = host_names[j],
            Species_range = sum(pred_current[,colnames(pred_current)==butterfly_names[i]]),
            Host_range = sum(pred_current[,colnames(pred_current)==host_names[j]]),
            Overlap_host = sum(pred_current[,colnames(pred_current)==butterfly_names[i]] * 
                                pred_current[,colnames(pred_current)==host_names[j]])))
    }

    range_fut26_df = rbind(range_fut26_df, df)
  }
  
  write.csv(range_fut26_df, file = paste0(output_dir, '/range_overlap_hosts_future_26.csv'), quote = F)
  
  
  #future 85

  range_fut85_df = data.frame()
  
  plant_names = colnames(species_matrix)
  butterfly_names = species_matrix$X
  butterfly_names = butterfly_names[butterfly_names %in% colnames(pred_current)]

  for (i in 1:length(butterfly_names)) {

    df = data.frame()

    host_names = plant_names[species_matrix |> filter(X == butterfly_names[i])==1] 

    for (j in 1:length(host_names)) {

      df = rbind(df, cbind(
            Species = butterfly_names[i],
            Host = host_names[j],
            Species_range = sum(pred_current[,colnames(pred_current)==butterfly_names[i]]),
            Host_range = sum(pred_current[,colnames(pred_current)==host_names[j]]),
            Overlap_host = sum(pred_current[,colnames(pred_current)==butterfly_names[i]] * 
                                pred_current[,colnames(pred_current)==host_names[j]])))
    }

    range_fut85_df = rbind(range_fut85_df, df)
  }
  
  write.csv(range_fut85_df, file = paste0(output_dir, '/range_overlap_hosts_future_85.csv'), quote = F)
  
