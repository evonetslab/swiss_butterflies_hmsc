GetData = function(papilio_test, plot_maps, epsg_temp, fact_pred_temp) {
  
  if (papilio_test) {
    data_suffix = 'papilio_test'
    plot_dir = 'output/papilio'
  } else {
    data_suffix = 'all'
    plot_dir = 'output/all'
  }
  
  if (file.exists(paste0('/scratch/project_2011416/Data/Compiled_data/Data_', data_suffix, '_v6.RData'))) {
    
    load(paste0('/scratch/project_2011416/Data/Compiled_data/Data_', data_suffix, '_v6.RData'))
    
  } else {
    
    ####BUTTERFLIES####
    #Read species interaction matrix (species on the rows and host genera on the columns)
    if (papilio_test) {
      species_matrix = read.csv('/scratch/project_2011416/Data/papilio/interaction_matrix_papilio.csv')
    } else {
      species_matrix = read.csv('/scratch/project_2011416/Data/interaction_matrix.csv')
    }
    
    #Butterfly study species
    butterfly_sp = species_matrix$X
    
    #Read butterfly occurrence matrix 
    # butterfly_occ = read.csv('/scratch/project_2011416/Data/infofauna_1x1km/infofauna_1x1km.csv', sep = ';', header = T)
    butterfly_occ = read.csv('/scratch/project_2011416/Data/Species/infofauna_1X1km.csv', header = T)
    
    #Filter temporally 2000-2022
    butterfly_occ = butterfly_occ[butterfly_occ$JAHR>=2000,]
    
    #Filter data by study species
    # butterfly_occ$scientificName = apply(matrix(butterfly_occ$TAXON), 1, function(x) paste0(unlist(str_split(x, ' '))[1], '_', unlist(str_split(x, ' '))[2]))
    # t$scientificName = apply(matrix(t$TAXON), 1, function(x) paste0(unlist(str_split(x, ' '))[1], '_', unlist(str_split(x, ' '))[2]))
    butterfly_occ$scientificName = str_replace(butterfly_occ$TAXON, ' ', '_')
    
    #add a column for occurrence
    butterfly_occ$occ = 1
    
    #filter species observations with spatial precision
    ind = which(butterfly_occ$RADIUS == '250 - 1000' & !is.na(butterfly_occ$RAUM))
    butterfly_occ = butterfly_occ[ind,]
    
    #turn into a data frame with species as columns
    butterfly_table = cast(butterfly_occ, KOORDX + KOORDY ~ scientificName, value = 'occ', sum)
    butterfly_table = butterfly_table[,c(1:2, which(colnames(butterfly_table) %in% butterfly_sp))]
    colnames(butterfly_table)[1:2] = c('x', 'y')
    
    #number of visits at each site
    survey_effort = matrix(NA, nrow = nrow(butterfly_table))
    for (i in 1:nrow(butterfly_table)) {
      print(i)
      temp = butterfly_occ[butterfly_occ$KOORDX==butterfly_table$x[i] &
                                         butterfly_occ$KOORDY==butterfly_table$y[i],]
      
      survey_effort[i] = nrow(unique(select(temp, c(JAHR, MONAT, TAG))))
    }
    
    butterfly_table$survey_effort_butterfly = survey_effort
    
    ####PLANTS####
    ##Repeat the same procedure to the plant data
    #Plant study genera
    plant_gen = colnames(species_matrix)[-1]
    
    #Get plant data
    # plant_occ = read.table('/scratch/project_2011416/Data/infoflora_1x1km/infoflora_1x1km.tab', sep = '\t', header = T)
    plant_occ = read.table('/scratch/project_2011416/Data/Species/infoflora_1x1km.tab', sep = '\t', header = T)
    
    #Filter temporally 2000-2025
    plant_occ = plant_occ[plant_occ$first_observation>=2000,]
    
    #add a column for occurrence
    plant_occ$occ = 1
    
    #get species genus from names
    plant_occ$genus = unlist(purrr::map(plant_occ$taxon, function(x) str_split_1(x, ' ')[1]))
    
    #turn into a data frame with species as columns
    plant_table = cast(plant_occ, x + y ~ genus, value = 'occ', sum)
    
    #keep only species from the genus
    plant_table = plant_table[,c(1:2, which(colnames(plant_table) %in% plant_gen))]
    
    #include a survey effort variable
    plant_table$survey_effort_plant = 0
    
    
    ####MERGE butterfly and plant tables####
    # sp_table = union_all(butterfly_table, plant_table)
    sp_table = bind_rows(butterfly_table, plant_table)
    
    #group by X and Y values
    sp_table = sp_table %>%
      group_by(x,y) %>%
      summarise(across(colnames(sp_table)[-c(1:2)], function(x) ifelse(all(is.na(x)), NA, sum(x, na.rm = T))))
    
    sp_table = data.frame(sp_table)
    
    #pick survey effort from the table
    survey_effort_butterfly = sp_table$survey_effort_butterfly
    survey_effort_plant = sp_table$survey_effort_plant
    sp_table[,'survey_effort_butterfly'] = NULL
    sp_table[,'survey_effort_plant'] = NULL
    
    #Remove cells with any plant has an na-observation
    #This needs to be updated
    # ind = which(apply(sp_table[,c(1:2,(ncol(butterfly_table)+1):ncol(sp_table))], 1, function(x) any(is.na(x))))
    # sp_table = sp_table[-ind,]
    # survey_effort_butterfly = survey_effort_butterfly[-ind]
    # survey_effort_plant = survey_effort_plant[-ind]

    ind = which(!is.na(survey_effort_butterfly) & !is.na(survey_effort_plant))
    sp_table = sp_table[ind,]
    survey_effort_butterfly = survey_effort_butterfly[ind]
    survey_effort_plant = survey_effort_plant[ind]
    
    #Get environmental covariates for the sampling locations
    #Transform sampling coordinates - the current coordinate reference system is CH1903+ / LV95 -- Swiss CH1903+ / LV95
    #EPSG:2056
    points_df = sp_table[,1:2]
    points_poly = st_as_sf(points_df, coords = c('x', 'y'), crs = epsg_temp)
    # plane_projection = 'epsg:3035'
    # points_poly = st_transform(points_poly, plane_projection)
    
    #read world data
    world = ne_countries(returnclass = "sf", scale = 'medium')
    world_plane = st_transform(world, epsg_temp)
    swiss_poly = world_plane$geometry[world_plane$sov_a3=='CHE']
    #Remove points from outside of Switzerland
    temp = st_within(points_poly, swiss_poly, sparse = F)
    sp_table = sp_table[temp,]
    survey_effort_butterfly = survey_effort_butterfly[temp]
    survey_effort_plant = survey_effort_plant[temp]
    points_poly = points_poly[temp,]
    
    #set domain corresponding to Switzerland
    domain = st_union(st_buffer(swiss_poly, dist = 100000))
    domain_wgs = st_transform(domain, crs = 'epsg:4326')
    
    # #get extent of the domain
    # extent_domain = st_bbox(domain_wgs)
    #read water body
    r = rast('/scratch/project_2011416/Data/Environment/Water_body/consensus_full_class_12.tif')
    #crop water body by extent
    r = crop(r, domain_wgs)
    #reproject water body
    domain = terra::project(r, paste0('epsg:', epsg_temp))
    
    #change land areas to 1 and water=100 areas 0
    domain[domain < 100] = 1
    domain[domain == 100] = 0
    
    #write domain raster
    writeRaster(domain, filename = '/scratch/project_2011416/Data/Environment/Domain/Domain.tif', overwrite = T)
    
    #get the resolution
    raster_res = res(domain)
    
    #read re-scaled rasters or re-scale original rasters
    if (file.exists('/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled.tif')) {
      raster_cov_t = list()
      raster_cov_t[[1]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled.tif')
      raster_cov_t[[2]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled_future_26.tif')
      raster_cov_t[[3]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled_future_85.tif')
      raster_cov_t[[4]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled.tif')
      raster_cov_t[[5]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled_future_26.tif')
      raster_cov_t[[6]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled_future_85.tif')
      raster_cov_t[[7]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled.tif')
      raster_cov_t[[8]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled_future_26.tif')
      raster_cov_t[[9]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled_future_85.tif')
      raster_cov_t[[10]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled.tif')
      raster_cov_t[[11]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled_future_26.tif')
      raster_cov_t[[12]] = rast('/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled_future_85.tif')
      raster_cov_t[[13]] = rast('/scratch/project_2011416/Data/Environment/Elevation/elevation_1KMmd_GMTEDmd_scaled.tif')
      raster_cov_t[[14]] = rast('/scratch/project_2011416/Data/Environment/Cloud_cover/MODCF_intraannualSD_scaled.tif')
      raster_cov_t[[15]] = rast('/scratch/project_2011416/Data/Environment/TRI/tri_1KMmd_GMTEDmd_scaled.tif')
      raster_cov_t[[16]] = rast('/scratch/project_2011416/Data/Environment/Reclassified-CLC-2018/LANDUSECLC1KM-roughedges-2018_scaled.tif')
      raster_cov_t[[17]] = rast('/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP2.6/rcp26yr2063_scaled.tif')
      raster_cov_t[[18]] = rast('/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP8.5/RCP85yr2063_scaled.tif')
      water = rast('/scratch/project_2011416/Data/Environment/Water_body/consensus_full_class_12_scaled.tif')
    } else {
      #get environmental covariates
      bioclim = apply(as.matrix(list.files('/scratch/project_2011416/Data/Environment/Bioclim')[-13]), 1, function(x) rast(paste0('/scratch/project_2011416/Data/Environment/Bioclim/', x)))
      #reorder bioclim layers
      bioclim = list(bioclim[[1]], bioclim[[2]], bioclim[[3]],
                     bioclim[[10]], bioclim[[11]], bioclim[[12]],
                     bioclim[[4]], bioclim[[5]], bioclim[[6]],
                     bioclim[[7]], bioclim[[8]], bioclim[[9]])
      
      elevation = rast('/scratch/project_2011416/Data/Environment/Elevation/elevation_1KMmd_GMTEDmd.tif')
      cloud = rast('/scratch/project_2011416/Data/Environment/Cloud_cover/MODCF_intraannualSD.tif')
      tri = rast('/scratch/project_2011416/Data/Environment/TRI/tri_1KMmd_GMTEDmd.tif')
      water = rast('/scratch/project_2011416/Data/Environment/Water_body/consensus_full_class_12.tif')
      land_cover = rast('/scratch/project_2011416/Data/Environment/Reclassified-CLC-2018/LANDUSECLC1KM-roughedges-2018.tif')
      land_cover_future_26 = rast('/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP2.6/rcp26yr2063.tif')
      land_cover_future_85 = rast('/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP8.5/RCP85yr2063.tif')
      
      #reproject land cover
      land_cover = project(land_cover, 'epsg:4326', method = 'mode')
      land_cover_future_26 = project(land_cover_future_26, 'epsg:4326', method = 'mode')
      land_cover_future_85 = project(land_cover_future_85, 'epsg:4326', method = 'mode')
      
      #crop covariates with the geographic extent
      land_cover = crop(land_cover, domain_wgs)
      land_cover_future_26 = crop(land_cover_future_26, domain_wgs)
      land_cover_future_85 = crop(land_cover_future_85, domain_wgs)
      
      #reproject covariates
      land_cover = project(land_cover, domain, method = 'mode')
      land_cover_future_26 = project(land_cover_future_26, domain, method = 'mode')
      land_cover_future_85 = project(land_cover_future_85, domain, method = 'mode')
      
      #set the same resolution for all rasters
      land_cover = resample(land_cover, domain, method = 'mode')
      land_cover_future_26 = resample(land_cover_future_26, domain, method = 'mode')
      land_cover_future_85 = resample(land_cover_future_85, domain, method = 'mode')
      
      #combine all raster layers into one stack
      raster_cov = c(bioclim, cloud, elevation, tri, water)  
      
      #crop covariates with the geographic extent
      raster_cov = lapply(raster_cov, function(x) crop(x, domain_wgs))
      
      #reproject covariates
      raster_cov = lapply(raster_cov, function(x) project(x, domain))
      
      #set the same resolution for all rasters
      raster_cov_t = lapply(raster_cov, function(x) resample(x, domain))
      
      writeRaster(raster_cov_t[[1]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[2]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled_future_26.tif', overwrite = T)
      writeRaster(raster_cov_t[[3]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio1_scaled_future_85.tif', overwrite = T)
      writeRaster(raster_cov_t[[4]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[5]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled_future_26.tif', overwrite = T)
      writeRaster(raster_cov_t[[6]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio2_scaled_future_85.tif', overwrite = T)
      writeRaster(raster_cov_t[[7]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[8]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled_future_26.tif', overwrite = T)
      writeRaster(raster_cov_t[[9]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio12_scaled_future_85.tif', overwrite = T)
      writeRaster(raster_cov_t[[10]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[11]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled_future_26.tif', overwrite = T)
      writeRaster(raster_cov_t[[12]], file = '/scratch/project_2011416/Data/Environment/Bioclim/Bio15_scaled_future_85.tif', overwrite = T)
      writeRaster(raster_cov_t[[13]], file = '/scratch/project_2011416/Data/Environment/Cloud_cover/MODCF_intraannualSD_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[14]], file = '/scratch/project_2011416/Data/Environment/Elevation/elevation_1KMmd_GMTEDmd_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[15]], file = '/scratch/project_2011416/Data/Environment/TRI/tri_1KMmd_GMTEDmd_scaled.tif', overwrite = T)
      writeRaster(land_cover, file = '/scratch/project_2011416/Data/Environment/Reclassified-CLC-2018/LANDUSECLC1KM-roughedges-2018_scaled.tif', overwrite = T)
      writeRaster(land_cover_future_26, file = '/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP2.6/rcp26yr2063_scaled.tif', overwrite = T)
      writeRaster(land_cover_future_85, file = '/scratch/project_2011416/Data/Environment/Projection-Land-use-1km-2015-2100/RCP8.5/RCP85yr2063_scaled.tif', overwrite = T)
      writeRaster(raster_cov_t[[16]], file = '/scratch/project_2011416/Data/Environment/Water_body/consensus_full_class_12_scaled.tif', overwrite = T)
      water = raster_cov_t[[16]]
      raster_cov_t[[16]] = land_cover
      raster_cov_t[[17]] = land_cover_future_26
      raster_cov_t[[18]] = land_cover_future_85
    }
    
    #set all cells where any covariate has NA or water = 100
    #get coordinates as a data frame
    train_coordinates = do.call(rbind, st_geometry(points_poly)) %>% 
      as_tibble() %>% setNames(c("x","y")) %>% data.frame()
    
    #turn list into one raster stack
    raster_cov_t = rast(raster_cov_t)
    train_covariates = terra::extract(raster_cov_t, train_coordinates, ID = F)
    train_covariates = as.matrix(train_covariates)
    
    #turn list into a matrix
    train_covariates[is.nan(train_covariates)] = NA
    
    #find rows where any covariate has NA value
    ind_na_covariates = apply(train_covariates, 1, anyNA)
    #find rows which have water body 100
    ind_water = terra::extract(water, train_coordinates, ID = F)[[1]]==100
    #take only points which have true in either ind:s
    ind_na_covariates = (ind_na_covariates + ind_water)>0
    
    #split the data to current and future bioclimatic and land cover layers
    bioclim_future = train_covariates[,c(2:3,5:6,8:9,11:12)]
    land_cover_future = train_covariates[,17:18]
    
    train_covariates = train_covariates[,-c(2:3,5:6,8:9,11:12,17:18)]
    
    colnames(train_covariates) = c('bio1', 'bio2', 'bio12', 'bio15', 'cloud', 'elevation', 'tri', 'land_cover')
    
    #standardize the covariates
    cov_mean = apply(train_covariates[,-ncol(train_covariates)],2,mean)
    cov_sd = apply(train_covariates[,-ncol(train_covariates)],2,sd)
    train_covariates_st = cbind(apply(train_covariates[,-ncol(train_covariates)], 2, function(x) (x-mean(x))/sd(x)),
                                train_covariates[,ncol(train_covariates)])
    colnames(train_covariates_st) = c('bio1', 'bio2', 'bio12', 'bio15', 'cloud', 'elevation', 'tri', 'land_cover')
    
    #create a data frame with coordinates and covariate values
    train_covariates_st = cbind(train_coordinates, train_covariates_st)
    train_covariates_df = data.frame(train_covariates_st)
    
    #remove NA points from environmental and species tables
    train_covariates_df = train_covariates_df[!ind_na_covariates,]
    sp_table = sp_table[!ind_na_covariates,]
    survey_effort_butterfly = survey_effort_butterfly[!ind_na_covariates]
    survey_effort_plant = survey_effort_plant[!ind_na_covariates]
    train_coordinates = train_coordinates[!ind_na_covariates,]
    
    #remove coordinates from species and covariate matrices
    sp_table = sp_table[,-c(1:2)]
    train_covariates_df = train_covariates_df[,-c(1,2)]
    
    #store data into a list
    PA_training_data = list(train_coordinates, sp_table, train_covariates_df, cov_mean, cov_sd, survey_effort_butterfly, survey_effort_plant)
    names(PA_training_data) = c('coordinates', 'species', 'covariates', 'cov_mean', 'cov_sd', 'survey_effort_butterfly', 'survey_effort_plant')
    
    ## Create data in 1km resolution for predictive maps
    # coarsen resolution of the predictive maps
    if (fact_pred_temp > 1) {
      raster_cov_coarse_env = aggregate(raster_cov_t[[1:15]], fact = fact_pred_temp, fun = 'median')
      raster_cov_coarse_lc = aggregate(raster_cov_t[[16:18]], fact = fact_pred_temp, fun = 'mode')
      raster_cov_coarse = c(raster_cov_coarse_env, raster_cov_coarse_lc)
      
      #aggregate domain
      domain_coarse = aggregate(domain, fact = fact_pred_temp, fun = 'mean')
      values(domain_coarse) = ifelse(values(domain_coarse)>=.5, 1, NA)
      writeRaster(domain_coarse, filename = '/scratch/project_2011416/Data/Environment/Domain/Domain_pred.tif', overwrite = T)
      
      #crop coarse raster with coarse domain
      raster_cov_coarse = lapply(raster_cov_coarse, function(x) crop(x, domain_coarse))
    } else {
      raster_cov_coarse = raster_cov_t
      domain_coarse = domain
      writeRaster(domain_coarse, filename = '/scratch/project_2011416/Data/Environment/Domain/Domain_pred.tif', overwrite = T)
    }

    
    #transform rasters into a data frame and get cell coordinates
    cov_df = as.matrix(raster_cov_coarse, na.rm = F)
    cov_df[is.nan(cov_df)] = NA
    cov_na = apply(cov_df,1,anyNA)
    domain_na = is.na(values(domain_coarse))
    ind_na = (cov_na + domain_na) > 0
    cov_matrix = cov_df[!ind_na,]
    # cov_matrix = cbind(cov_matrix, cov_matrix^2)
    #create matrices for current and future conditions
    #keep only bioclim of current or future conditions
    cov_matrix_current = cov_matrix[,-c(2:3,5:6,8:9,11:12,17:18)]
    cov_matrix_future_26 = cov_matrix[,-c(1,3,4,6,7,9,10,12,16,18)]
    cov_matrix_future_85 = cov_matrix[,-c(1:2,4:5,7:8,10:11,16:17)]
    
    colnames(cov_matrix_current) = colnames(train_covariates)
    colnames(cov_matrix_future_26) = colnames(train_covariates)
    colnames(cov_matrix_future_85) = colnames(train_covariates)
    
    #get coordinates of the cells
    xy_matrix = crds(domain_coarse, df = F, na.rm = F)[!ind_na,]
    colnames(xy_matrix) = c('X', 'Y')
    
    #scale covariate values
    for (i in 1:(ncol(cov_matrix_current)-1)) {
      cov_matrix_current[,i] = (cov_matrix_current[,i]-cov_mean[i])/cov_sd[i]
      cov_matrix_future_26[,i] = (cov_matrix_future_26[,i]-cov_mean[i])/cov_sd[i]
      cov_matrix_future_85[,i] = (cov_matrix_future_85[,i]-cov_mean[i])/cov_sd[i]
    }
    
    pred_data = list(coordinates = xy_matrix, covariates_current = cov_matrix_current, covariates_future_26 = cov_matrix_future_26, covariates_future_85 = cov_matrix_future_85, ind_na = ind_na)
    
    # #for INLA    
    # hull = inla.nonconvex.hull(as.matrix(PA_training_data$coordinates)/1000, convex = -0.01, resolution = 500)
    # mesh = inla.mesh.2d(boundary = hull, cutoff = 4, offset = c(20, 40), max.edge = c(15,30))
    
    data = list(PA_training_data, pred_data)
    names(data) = c('PA_training_data', 'pred_data')
    
    save(data, file = paste0('/scratch/project_2011416/Data/Compiled_data/Data_', data_suffix, '_v6.RData'))
    
  }
  return(data)
  
  if (plot_maps) {
    
    ##Plot presence points
    # swissfill = map("world", "switzerland", fill=TRUE, plot=FALSE)
    # IDs = sapply(strsplit(swissfill$names, ":"), function(x) x[1])
    # swiss.poly = map2SpatialPolygons(swissfill, IDs = IDs, proj4string = CRS('+init=epsg:4326'))
    # swiss.poly = spTransform(swiss.poly, CRS(paste0('+init=',data$proj_raster)))
    # 
    # 
    #turn species matrix into a spatial points data frame
    #Read species interaction matrix (species on the rows and host genera on the columns)
    if (papilio_test) {
      species_matrix = read.csv('/scratch/project_2011416/Data/papilio/interaction_matrix_papilio.csv')
    } else {
      species_matrix = read.csv('/scratch/project_2011416/Data/interaction_matrix.csv')
    }
    
    #study species/genus
    butterfly_sp = species_matrix$X
    plant_gen = colnames(species_matrix)[-1]
    
    #butterflies
    butterflies_coord = NULL
    butterflies_name = NULL
    prevalence = NULL
    
    for (i in butterfly_sp) {
      row_ind = data$PA_training_data$species[,colnames(data$PA_training_data$species) %in% i]==1 &
        !is.na(data$PA_training_data$species[,colnames(data$PA_training_data$species) %in% i])
      
      butterflies_coord = rbind(butterflies_coord, data$PA_training_data$coordinates[row_ind,])
      butterflies_name = c(butterflies_name, rep(i,sum(row_ind)))
      
      prevalence = c(prevalence, sum(row_ind))
    }
    
    temp_data = data.frame(lon = butterflies_coord[,1], lat = butterflies_coord[,2], species = butterflies_name)
    ann_text = data.frame(lon = 2550000, lat = 1300000, prev = prevalence, species = butterfly_sp)
    
    g1 = ggplot(data = temp_data) + 
      geom_sf(data = swiss_poly) + 
      geom_point(aes(x=lon, y=lat)) +
      facet_wrap(vars(species)) +
      labs(x = 'Longitude', y = 'Latitude', col = 'Grid Observation') +
      ggtitle('Butterflies') +
      theme_light() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_text(data = ann_text, mapping = aes(x = lon, y = lat, label = prev))
    
    pdf(file = paste0(plot_dir, '/Butterflies_occ_v5.pdf'), width = 10, height = 10)
    g1
    dev.off()
    
    
    #plants
    plants_coord = NULL
    plants_name = NULL
    prevalence = NULL
    for (i in plant_gen) {
      row_ind = data$PA_training_data$species[,colnames(data$PA_training_data$species) %in% i]==1 &
        !is.na(data$PA_training_data$species[,colnames(data$PA_training_data$species) %in% i])
      
      plants_coord = rbind(plants_coord, data$PA_training_data$coordinates[row_ind,])
      plants_name = c(plants_name, rep(i,sum(row_ind)))
      
      prevalence = c(prevalence, sum(row_ind))
    }
    
    temp_data = data.frame(lon = plants_coord[,1], lat = plants_coord[,2], species = plants_name)
    ann_text = data.frame(lon = 2550000, lat = 1300000, prev = prevalence, species = plant_gen)
    
    g2 = ggplot(data = temp_data) + 
      geom_sf(data = swiss_poly) + 
      geom_point(aes(x=lon, y=lat)) +
      facet_wrap(vars(species)) +
      labs(x = 'Longitude', y = 'Latitude', col = 'Grid Observation') +
      ggtitle('Plants') +
      theme_light() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_text(data = ann_text, mapping = aes(x = lon, y = lat, label = prev))
    
    pdf(file = paste0(plot_dir, '/Plants_occ_v5.pdf'), width = 10, height = 10)
    g2
    dev.off()
    
  }
}


#how many cells have plants and butterflies surveys - both should give the same result
# sum(apply(sp_table[,3:8], 1, function(x) !all(is.na(x))) & apply(sp_table[,9:22], 1, function(x) !all(is.na(x))))
# sum(apply(matrix(sp_table[,9]), 1, function(x) !is.na(x)) & apply(matrix(sp_table[22]), 1, function(x) !is.na(x)))
