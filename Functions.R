############## Gets the max drainage value (from most to least by x) to get their streamflow data ##########################
max_drainage <- function(drainage, x) { # x is the number of places from the max
  drain <- drainage[order(drainage$Drainage, na.last = FALSE) , ] # orders the dataframe by drainage (increasing top to bottom) and puts NAs at top
  row <- drain[length(drain$Drainage) - x, ] # Grabs the desired row (first max, second, third, etc...)
  return(row)
}

############# Tests the quality of the Streamflow data that is supplied to it ############################
streamflow_quality <- function(our_station) {
  start_date <- "2000-01-01"
  end_date <- "2019-12-31"
  if (is.na(our_station$Station[1])) {
    streamflow <- NULL
    return(streamflow)
  }
  streamflow <- try(readNWISdv(our_station$Station[1], parameterCd = '00060', statCd = '00003', startDate = start_date, endDate = end_date))
  if (nrow(streamflow) != 0) {
    return(streamflow)
  } else if (nrow(streamflow) == 0) {
    streamflow <- NULL
    return(streamflow)
  }
}


########################### Read all Streamflow Data ###########################
read_streamflow <- function(streamflow_files) {
    all_streamflow <- read.csv(streamflow_files[1], header = FALSE) # Insert the first column of data to allow for consecutive merges
    parse_1 <- str_split_1(streamflow_files[1], '_') # Parse the filename into its pieces
    all_streamflow <- all_streamflow %>%
      dplyr::select(-1) %>%
      rename(Date = V2,
             !!parse_1[4] := V3) # Names the column with the data as "HUC"
    for (i in 2:(length(streamflow_files) - 120)) {
      local <- read.csv(streamflow_files[i], header = FALSE)
      parsed <- str_split_1(streamflow_files[i], '_') # Parse the file name
      local <- local %>%
        dplyr::select(-1) %>% # Remove the column that counts the observations
        rename(Date = V2,
               !!parsed[4] := V3) # Names the column with the data as "ID_HUC"
      as.Date(local$Date, format = "%Y%m%d")
      all_streamflow <- merge(all_streamflow, local, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE) # Join the dataframes by date, keeping all observations
      rm(local)
    }
    all_streamflow <- all_streamflow %>%
      dplyr::select(-c(22,23))
    for (i in 153:length(streamflow_files)) {
      local <- read.csv(streamflow_files[i], header = FALSE)
      parsed <- str_split_1(streamflow_files[i], '_') # Parse the file name
      local <- local %>%
        dplyr::select(-1) %>% # Remove the column that counts the observations
        rename(Date = V2,
               !!parsed[4] := V3) # Names the column with the data as "ID_HUC"
      as.Date(local$Date, format = '%Y%m%d')
      all_streamflow <- merge(all_streamflow, local, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE) # Join the dataframes by date, keeping all observations
      rm(local)
    }
    return(all_streamflow)
} 


################# Make Dataframe of SNOTEL Station Quality #####################

snotel_quality <- function(snotel_files) {
  snotel_quality <- data.frame(matrix(nrow = 0, ncol = 11))
  colnames(snotel_quality) <- c('Filename', 'Id', 'Lat', 'Lon', 'HUC', 'SWE_Count', 'All_Count', 'Total_Observations', 'SWE_Percent', 'All_Percent', 'SWE_Quality')
  
  for (i in snotel_files) {
    local <- data.frame(matrix(nrow = 1, ncol = 11))
    colnames(local) <- c('Filename', 'Id', 'Lat', 'Lon', 'HUC', 'SWE_Count', 'All_Count', 'Total_Observations', 'SWE_Percent', 'All_Percent', 'SWE_Quality')
    snotel <- read.csv(i, header = FALSE)
    parsed <- str_split_1(i, '_')
    local$Filename <- i
    local$Id <- parsed[1]
    local$Lat <- as.numeric(parsed[2])
    local$Lon <- as.numeric(parsed[3])
    local$HUC <- parsed[4]
    local$SWE_Count <- nrow(snotel[!is.na(snotel$V4), ]) # number of non-NA precip observations
    local$All_Count <- nrow(snotel[!is.na(snotel$V4) & !is.na(snotel$V3) & !is.na(snotel$V5) & !is.na(snotel$V6) , ]) # number of days when precip, Airtemp, precip, and VSM_2 are valid
    local$Total_Observations <- nrow(snotel)
    local$SWE_Percent <- round((local$SWE_Count/local$Total_Observations)*100 , digits = 1)
    local$All_Percent <- round((local$All_Count/local$Total_Observations)*100 , digits = 1)
    if (local$SWE_Percent < 20) {
      local$SWE_Quality <- 'Poor'
    } else if (local$SWE_Percent >= 20 & local$SWE_Percent < 50) {
      local$SWE_Quality <- 'Decent'
    } else if (local$SWE_Percent >= 50 & local$SWE_Percent < 70) {
      local$SWE_Quality <- 'Good'
    } else if (local$SWE_Percent >= 70) {
      local$SWE_Quality <- 'Great'
    }
    snotel_quality <- rbind(snotel_quality, local)
  }
  return(snotel_quality)
}


##################### Clip The Dataframe to the correct size ###################
clip_data <- function(df) {
  df$Year <- year(ymd(df$Date))
  df$Month <- month(ymd(df$Date))
  df$Day <- day(ymd(df$Date))
  df <- df %>%
    relocate(c('Year', 'Month', 'Day')) %>%
    relocate('Date') %>%
    filter(Date >= '2000-10-01') %>%
    filter(Date < '2019-10-01')
  return(df)
}

################# Creates and Writes a Streamflow File for the station with the max drainage area in each HUC8 ###################
HUC_streamflow <- function(unique_huc, HUC_no_data) {
  # Initialize our USGS station metadata dataframe
  USGS_metadata <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(USGS_metadata) <- c("Station", "Latitude", "Longitude", "Elevation", "HUC", "Drainage")
  
  for (i in unique_huc) { ## HUC8 14060006 does not return any data for 2000-2019? (contains station 461 in East Willow Utah)
    # Collect USGS Data of the stations in that HUC - looking for streamflow (parm_cd: 00060_00003 (mean daily discharge, cubic ft/second))
    huc_code <- paste0("", i, "")
    start_date <- "2000-01-01"
    end_date <- "2019-12-31"
    usgs <- whatNWISdata(huc = huc_code)
    
    # Remove all groundwater stations (ID length is > 8) and make a list of unique sites in this huc
    usgs <- usgs %>% 
      filter(nchar(site_no) == 8)
    usgs_stations <- unique(usgs$site_no)
    if(nrow(usgs) == 0) {
      print(huc_code)
      rm(usgs)
      next
    }
   
    # Collect the drainage area for each unique station in the HUC8
    drainage <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(drainage) <- c("Station", "Latitude", "Longitude", "Elevation", "HUC", "Drainage")
    for (j in usgs_stations) {
      data <- readNWISsite(j)
      if (nrow(data) > 1) {
        data <- data %>%
          filter(data$agency_cd == 'USGS')
      }
      if (nchar(j) != 8) {
        print(j)
        next
      }
      local <- data.frame(matrix(nrow = 1, ncol = 6))
      colnames(local) <- c("Station", "Latitude", "Longitude", "Elevation", "HUC", "Drainage")
      data$drain_area_va -> local$Drainage
      data$dec_lat_va -> local$Latitude
      data$dec_long_va -> local$Longitude
      data$site_no -> local$Station
      data$alt_va -> local$Elevation
      data$huc_cd -> local$HUC
      drainage <- rbind(drainage, local)
    }
    
    # Collect average daily streamflow (00060_00003) for the station with the max drainage area
    # This for loop will try to collect streamflow from the max drainage area in each HUC, but if there is none, it will go to the next largest drainage area
    for (x in 0:length(drainage)) {
      our_station <- max_drainage(drainage, x) # Should return a dataframe with a single row 
      streamflow <- streamflow_quality(our_station) # if there are multiple stations returned, this line keeps the station with the best daily data
      if (is.null(streamflow)) {
        rm(our_station)
        next # If streamflow_quality did not find any streamflow data for that station, go to the next highest drainage area station. 
      } else { break } # If streamflow was found, end the for loop 
    } 
    
    # If the entire HUC did not have any streamflow data, go to next HUC
    if (is_empty(streamflow)) {
      HUC_no_data <- append(HUC_no_data, huc_code)
      next
    } 
    
    # Bind our station (with the selected max drainage) to our USGS site metadata file
    USGS_metadata <- rbind(USGS_metadata, our_station)
    
    max_station <- streamflow$site_no[1]
    max_latitude <- our_station$Latitude
    max_longitude <- our_station$Longitude
    streamflow <- streamflow %>%
      dplyr::select(-c(1,2,5))
    names(streamflow) <- NULL
    # Write a csv file containing average daily average streamflow for this max-drainage station
    write_path <- "/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/.shortcut-targets-by-id/1GX8mj4pqw1SoB9JTetKChQp2Y0mgEAss/Ethan_Heidtman_Summer_2023/SWE/USGS_Streamflow/"
    file_name <- paste(max_station, max_latitude, max_longitude, huc_code, '.csv', sep = '_')
    write.csv(streamflow, file.path(write_path, file_name))
    # Clear the workspace for next loops
    rm(usgs, streamflow, max_station, drainage, max_longitude, max_latitude, usgs_stations)

    return(USGS_metadata)
  }     
  
}

############################# Normalizes a set of Data #############################################
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


############# Creates and Writes a Dataframe for Each Snotel Station ####################
station_files <- function(metadata, date, AirTemp, SWE, Precip, VSM_2, VSM_12, VSM_32) { # The arguments are metadata, date, and each of the created datasets
  for (i in 1:length(metadata$StationName)) {
    # Create temporary dataframe for each of the metadata with length of the time period
    station <- metadata$StationName[i]
    id <- metadata$StationId[i]
    state <- metadata$State[i]
    elevation <- metadata$Elevation[i]
    lat <- metadata$Latitude[i]
    lon <- metadata$Longitude[i]
    HUC8 <- metadata$HUC8[i]
    file <- date
    # Bind the actual data to the dataframe
    air_local <- AirTemp %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    swe_local <- SWE %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    precip_local <- Precip %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    vsm_2_local <- VSM_2 %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    vsm_12_local <- VSM_12 %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    vsm_32_local <- VSM_32 %>%
      dplyr::select(contains(paste(' ', as.character(id))))
    file <- cbind(file, air_local, swe_local, precip_local, vsm_2_local, vsm_12_local, vsm_32_local)
    if (ncol(file) > 7) { # Some columns were copied multiple times 
      file <- file %>%
        dplyr::select(-c(8:12))
    }
    names(file) <- NULL # remove the names of the colummns
    # Save the dataframe with the appropriate name 
    write_path <- "/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/.shortcut-targets-by-id/1GX8mj4pqw1SoB9JTetKChQp2Y0mgEAss/Ethan_Heidtman_Summer_2023/SWE/SWE_Station_Files/"
    file_name <- paste(id, lat, lon, HUC8, '.csv', sep = '_')
    write.csv(file, file.path(write_path, file_name))
    # Clear all the stuff to start again
    rm(station, id, state, elevation, lat, lon, HUC8, file_name) 
  }
}


########### Creates Column of Water Year when used with mutate ##################
wtr_yr <- function(water_year_date, start_month=10) {
  # Convert dates into POSIXlt
  date.posix <-  as.POSIXlt(water_year_date)
  # Year offset
  offset = ifelse(date.posix$mon >= start_month - 1, 1, 0)
  # Water year
  adj.year = date.posix$year + 1900 + offset
  # Return the water year
  adj.year
}

################ Creates a dataframe of peak SWE -> 5% swe by water year ################
peak_perc_swe <- function(a, b, c, d, e) {
  for (i in 2:length(a)) {
    local <- data.frame(matrix(nrow = 1, ncol = 4))
    colnames(local) <- c('Station', 'Date_Peak', 'Date_5', 'Duration')
    maxes <- a %>%
      filter(Water_Year == c)
    percs <- b %>%
      filter(Water_Year == c)
    max <- as.numeric(maxes[colnames(maxes) == colnames(a)[i]])
    five <- as.numeric(percs[colnames(percs) == colnames(a)[i]])
    if (is.na(max)) { # If there is no peak (i.e. all NA e observations)
      # max_date <- NA
      # perc_date <- NA
      local$Station <- as.numeric(colnames(a[i]))
      local$Date_Peak <- NA
      local$Date_5 <- NA
      local$Duration <- NA
      d <- rbind(d, local)
      rm(local)
      next
    }
    station <- e %>%
      dplyr::select(c(2, colnames(a)[i]))
    max_dates <- station$Date[station[2] == max]
    five_dates <- station$Date[station[2] <= five]
    if (length(max_dates) == 0) {
      local$Station <- as.numeric(colnames(a[i]))
      local$Date_Peak <- NA
      local$Date_5 <- NA
      local$Duration <- NA
      d <- rbind(d, local)
      rm(local)
      next
    }
    max_date <- as.Date(min(max_dates))
    if (length(five_dates) == 0) { # If no dates are less than the peak
      local$Station <- as.numeric(colnames(a[i]))
      local$Date_Peak <- max_date
      local$Date_5 <- as.Date(paste(c, '-09-30', sep = ''))
      local$Duration <- as.numeric(local$Date_5 - max_date)
      d <- rbind(d, local)
      rm(local)
      next
    }
    
    five_dates <- five_dates[five_dates >= max_date]
    if (length(five_dates) == 0) { # If the snowpack doesn't reach 5%, just use the end of the water c
      local$Date_5 <- as.Date(paste(c, '-09-30', sep = ''))
      local$Date_Peak <- max_date
      local$Duration <- as.numeric(local$Date_5 - local$Date_Peak)
      local$Station <- as.numeric(colnames(a[i]))
      d <- rbind(d, local)
      rm(local)
      next
    }
    
    five_date <- as.Date(min(five_dates))
    
    local$Station <- as.numeric(colnames(a[i]))
    local$Date_Peak <- max_date
    local$Date_5 <- five_date
    local$Duration <- as.numeric(five_date - max_date)
    d <- rbind(d, local)
    
    rm(local, maxes, percs, max, five, max_dates, five_dates, max_date, five_date)
    
  }  
  return(d)
}


###########  Average the Streamflow sites by HUC2 by Year ######################
mean_streamflow_func <- function(mean_streamflow) {
  mean_HUC2 <- mean_streamflow
  colnames(mean_HUC2)[2:251] <- substr(colnames(mean_HUC2[2:251]), 1, 2) # Makes the columns HUC2s by keeping only the first two digits of the HUC8

  # Sort by HUC2, average across all of them
  nine <- as.data.frame(mean_HUC2[,1])
  nine <- cbind(nine, mean_HUC2[names(mean_HUC2) == '09'])
  colnames(nine) <- c('Water_Year', 'HUC09')
  
  ten <- mean_HUC2[,1]
  ten <- cbind(ten, mean_HUC2[names(mean_HUC2) == '10'])
  ten <- as.data.frame(rowMeans(ten[, 2:length(ten)], na.rm = TRUE))
  colnames(ten) <- 'HUC10'
  
  eleven <- mean_HUC2[,1]
  eleven <- cbind(eleven, mean_HUC2[names(mean_HUC2) == '11'])
  eleven <- as.data.frame(rowMeans(eleven[, 2:length(eleven)], na.rm = TRUE))
  colnames(eleven) <- 'HUC11'
  
  # There is no HUC2 == 12
  # twelve <- mean_HUC2[,1]
  # twelve <- cbind(twelve, mean_HUC2[names(mean_HUC2) == '12'])
  # twelve <- as.data.frame(rowMeans(twelve[, 2:length(twelve)], na.rm = TRUE))
  # colnames(twelve) <- 'HUC12'
  
  thirteen <- mean_HUC2[,1]
  thirteen <- cbind(thirteen, mean_HUC2[names(mean_HUC2) == '13'])
  thirteen <- as.data.frame(rowMeans(thirteen[, 2:length(thirteen)], na.rm = TRUE))
  colnames(thirteen) <- 'HUC13'
  
  fourteen <- mean_HUC2[,1]
  fourteen <- cbind(fourteen, mean_HUC2[names(mean_HUC2) == '14'])
  fourteen <- as.data.frame(rowMeans(fourteen[, 2:length(fourteen)], na.rm = TRUE))
  colnames(fourteen) <- 'HUC14'
  
  fifteen <- mean_HUC2[,1]
  fifteen <- cbind(fifteen, mean_HUC2[names(mean_HUC2) == '15'])
  fifteen <- as.data.frame(rowMeans(fifteen[, 2:length(fifteen)], na.rm = TRUE))
  colnames(fifteen) <- 'HUC15'
  
  sixteen <- mean_HUC2[,1]
  sixteen <- cbind(sixteen, mean_HUC2[names(mean_HUC2) == '16'])
  sixteen <- as.data.frame(rowMeans(sixteen[, 2:length(sixteen)], na.rm = TRUE))
  colnames(sixteen) <- 'HUC16'
  
  seventeen <- mean_HUC2[,1]
  seventeen <- cbind(seventeen, mean_HUC2[names(mean_HUC2) == '17'])
  seventeen <- as.data.frame(rowMeans(seventeen[, 2:length(seventeen)], na.rm = TRUE))
  colnames(seventeen) <- 'HUC17'
  
  eighteen <- mean_HUC2[,1]
  eighteen <- cbind(eighteen, mean_HUC2[names(mean_HUC2) == '18'])
  eighteen <- as.data.frame(rowMeans(eighteen[, 2:length(eighteen)], na.rm = TRUE))
  colnames(eighteen) <- 'HUC18'
  
  rm(mean_HUC2)
  mean_HUC2 <- bind_cols(nine, ten, eleven, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen)
  rm(nine, ten, eleven, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen)
  
  # Round all values to 2 decimal places
  mean_HUC2 <- mean_HUC2 %>%
    mutate_if(is.numeric, round, digits = 2)
  mean_HUC2$Water_Year <- as.factor(mean_HUC2$Water_Year) # Make factor variable for plotting purposes
  return(mean_HUC2)
  
}

###### Other expand ######
other_expand <- function(a, b, xcolname) {
  hucs <- as.character(unique(a$HUC8))
  for (i in hucs) {
    if (i %in% colnames(b)) {
      df <- b[names(b) == i]
    df <- as.data.frame(t(df))
    colnames(df) <- c(paste(xcolname, '01', sep = '_'), paste(xcolname, '02', sep = '_'), paste(xcolname, '03', sep = '_'), paste(xcolname, '04', sep = '_'),
                      paste(xcolname, '05', sep = '_'), paste(xcolname, '06', sep = '_'), paste(xcolname, '07', sep = '_'), paste(xcolname, '08', sep = '_'),
                      paste(xcolname, '09', sep = '_'), paste(xcolname, '10', sep = '_'), paste(xcolname, '11', sep = '_'), paste(xcolname, '12', sep = '_'), 
                      paste(xcolname, '13', sep = '_'), paste(xcolname, '14', sep = '_'), paste(xcolname, '15', sep = '_'), paste(xcolname, '16', sep = '_'), 
                      paste(xcolname, '17', sep = '_'), paste(xcolname, '18', sep = '_'), paste(xcolname, '19', sep = '_'))
    rows <- which(a$HUC8 == i)
    for (x in rows) {
      a[x, 9:27] <- df
    }
    } else { next }
  }
  return(a)
}

#################### Create Expanded Metadata with Peak SWE and Streamflow by Water Year for each HUC8 #########################
expand_metadata <- function(metadata, x, y, z, xcolname, ycolname, zcolname) {
  
  if (missing(z)) {
    # Initialize new dataframe to add to our metadata
    new <- data.frame(matrix(nrow = 653, ncol = 38))
    colnames(new) <- c(paste(xcolname, '01', sep = '_'), paste(xcolname, '02', sep = '_'), paste(xcolname, '03', sep = '_'), paste(xcolname, '04', sep = '_'),
                       paste(xcolname, '05', sep = '_'), paste(xcolname, '06', sep = '_'), paste(xcolname, '07', sep = '_'), paste(xcolname, '08', sep = '_'),
                       paste(xcolname, '09', sep = '_'), paste(xcolname, '10', sep = '_'), paste(xcolname, '11', sep = '_'), paste(xcolname, '12', sep = '_'), 
                       paste(xcolname, '13', sep = '_'), paste(xcolname, '14', sep = '_'), paste(xcolname, '15', sep = '_'), paste(xcolname, '16', sep = '_'), 
                       paste(xcolname, '17', sep = '_'), paste(xcolname, '18', sep = '_'), paste(xcolname, '19', sep = '_'),
                       paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                       paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                       paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                       paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                       paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'))
    expanded <- cbind(metadata, new) # Add all the new cells to insert the data
    expanded$StationId <- as.character(expanded$StationId) # Make the Ids characters for joining
    # Create list of HUC8s that have no data
    bad_huc8 <- vector(mode = 'list', length = 0)
    
    # Add our first set of data to the metadata
    f <- unique(metadata$HUC8) # List of unique HUC8s
    for (i in f) {
      df <- expanded %>% # Get the SNOTEL stations that are in this huc 8
        filter(HUC8 == i)
      peak_year <- x[names(x) %in% df$StationId]
      names <- colnames(peak_year)
      peak_year <- as.data.frame(t(peak_year)) # Transpose the dataframe
      rownames(peak_year) <- names
      colnames(peak_year) <- c(paste(xcolname, '01', sep = '_'), paste(xcolname, '02', sep = '_'), paste(xcolname, '03', sep = '_'), paste(xcolname, '04', sep = '_'),
                               paste(xcolname, '05', sep = '_'), paste(xcolname, '06', sep = '_'), paste(xcolname, '07', sep = '_'), paste(xcolname, '08', sep = '_'),
                               paste(xcolname, '09', sep = '_'), paste(xcolname, '10', sep = '_'), paste(xcolname, '11', sep = '_'), paste(xcolname, '12', sep = '_'), 
                               paste(xcolname, '13', sep = '_'), paste(xcolname, '14', sep = '_'), paste(xcolname, '15', sep = '_'), paste(xcolname, '16', sep = '_'), 
                               paste(xcolname, '17', sep = '_'), paste(xcolname, '18', sep = '_'), paste(xcolname, '19', sep = '_'))
      peak_year$StationId <- rownames(peak_year) # Create a new column of the station ID
      peak_year <- peak_year %>%
        relocate(StationId)
      for (j in peak_year$StationId) {
        row <- which(expanded$StationId == j) # The row in the expanded data frame to add data to
        row2 <- which(peak_year$StationId == j) # The row in the peak_year data frame
        expanded[row, 9:27] <- peak_year[row2, 2:length(peak_year)]
      }
      rm(peak_year, row, row2, names, df)
    }
    
    #### Add the Next variable too: 2 options (if colnames are HUCs or StationIds)
    if (nchar(colnames(y[2])) == 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- y[names(y) == huc8] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 28:46] <- NA
          bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                         paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                         paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                         paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                         paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'))
        q$HUC8 <- rownames(q) # Create a new column of the HUC8
        q <- q %>%
          relocate(HUC8)
        expanded[row, 28:46] <- q[, 2:length(q)] # Insert the second variable's data into the correct cells of the expanded dataframe
        rm(row, huc8, q)
      } 
    }
    if (nchar(colnames(y[2])) < 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        #huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- y[names(y) == k] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 28:46] <- NA
          #bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                         paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                         paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                         paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                         paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'))
        q$StationId <- rownames(q) # Create a new column of the StationId
        q <- q %>%
          relocate(StationId)
        expanded[row, 28:46] <- q[, 2:length(q)] # Insert the summed streamflow data into the correct cells of the expanded dataframe
        rm(row, q)
      } 
    }
    return(expanded)
  }
  if (!missing(z)) {
    new <- data.frame(matrix(nrow = 653, ncol = 57))
    colnames(new) <- c(paste(xcolname, '01', sep = '_'), paste(xcolname, '02', sep = '_'), paste(xcolname, '03', sep = '_'), paste(xcolname, '04', sep = '_'),
                       paste(xcolname, '05', sep = '_'), paste(xcolname, '06', sep = '_'), paste(xcolname, '07', sep = '_'), paste(xcolname, '08', sep = '_'),
                       paste(xcolname, '09', sep = '_'), paste(xcolname, '10', sep = '_'), paste(xcolname, '11', sep = '_'), paste(xcolname, '12', sep = '_'), 
                       paste(xcolname, '13', sep = '_'), paste(xcolname, '14', sep = '_'), paste(xcolname, '15', sep = '_'), paste(xcolname, '16', sep = '_'), 
                       paste(xcolname, '17', sep = '_'), paste(xcolname, '18', sep = '_'), paste(xcolname, '19', sep = '_'),
                       paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                       paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                       paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                       paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                       paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'), 
                       paste(zcolname, '01', sep = '_'), paste(zcolname, '02', sep = '_'), paste(zcolname, '03', sep = '_'), paste(zcolname, '04', sep = '_'),
                       paste(zcolname, '05', sep = '_'), paste(zcolname, '06', sep = '_'), paste(zcolname, '07', sep = '_'), paste(zcolname, '08', sep = '_'),
                       paste(zcolname, '09', sep = '_'), paste(zcolname, '10', sep = '_'), paste(zcolname, '11', sep = '_'), paste(zcolname, '12', sep = '_'), 
                       paste(zcolname, '13', sep = '_'), paste(zcolname, '14', sep = '_'), paste(zcolname, '15', sep = '_'), paste(zcolname, '16', sep = '_'), 
                       paste(zcolname, '17', sep = '_'), paste(zcolname, '18', sep = '_'), paste(zcolname, '19', sep = '_'))
  
    expanded <- cbind(metadata, new) # Add all the new cells to insert the data
    expanded$StationId <- as.character(expanded$StationId) # Make the Ids characters for joining
    
    # Create list of HUC8s that have no data
    bad_huc8 <- vector(mode = 'list', length = 0)
    
    # Add our first set of data to the metadata
    f <- unique(metadata$HUC8) # List of unique HUC8s
    for (i in f) {
      df <- expanded %>% # Get the SNOTEL stations that are in this huc 8
        filter(HUC8 == i)
      peak_year <- x[names(x) %in% df$StationId]
      names <- colnames(peak_year)
      
      peak_year <- as.data.frame(t(peak_year)) # Transpose the dataframe
      rownames(peak_year) <- names
      colnames(peak_year) <- c(paste(xcolname, '01', sep = '_'), paste(xcolname, '02', sep = '_'), paste(xcolname, '03', sep = '_'), paste(xcolname, '04', sep = '_'),
                               paste(xcolname, '05', sep = '_'), paste(xcolname, '06', sep = '_'), paste(xcolname, '07', sep = '_'), paste(xcolname, '08', sep = '_'),
                               paste(xcolname, '09', sep = '_'), paste(xcolname, '10', sep = '_'), paste(xcolname, '11', sep = '_'), paste(xcolname, '12', sep = '_'), 
                               paste(xcolname, '13', sep = '_'), paste(xcolname, '14', sep = '_'), paste(xcolname, '15', sep = '_'), paste(xcolname, '16', sep = '_'), 
                               paste(xcolname, '17', sep = '_'), paste(xcolname, '18', sep = '_'), paste(xcolname, '19', sep = '_'))
      
      peak_year$StationId <- rownames(peak_year) # Create a new column of the station ID
      peak_year <- peak_year %>%
        relocate(StationId)
      for (j in peak_year$StationId) {
        row <- which(expanded$StationId == j) # The row in the expanded data frame to add data to
        row2 <- which(peak_year$StationId == j) # The row in the peak_year data frame
        expanded[row, 9:27] <- peak_year[row2, 2:length(peak_year)]
      }
      rm(peak_year, row, row2, names, df)
    }
    
    #### Add the Next variable too: 2 options (if colnames are HUCs or StationIds)
    if (nchar(colnames(y[2])) == 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- y[names(y) == huc8] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 28:46] <- NA
          bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                         paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                         paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                         paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                         paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'))
        q$HUC8 <- rownames(q) # Create a new column of the HUC8
        q <- q %>%
          relocate(HUC8)
        expanded[row, 28:46] <- q[, 2:length(q)] # Insert the second variable's data into the correct cells of the expanded dataframe
        rm(row, huc8, q)
      } 
    }
    
    if (nchar(colnames(y[2])) < 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        #huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- y[names(y) == k] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 28:46] <- NA
          #bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(ycolname, '01', sep = '_'), paste(ycolname, '02', sep = '_'), paste(ycolname, '03', sep = '_'), paste(ycolname, '04', sep = '_'),
                         paste(ycolname, '05', sep = '_'), paste(ycolname, '06', sep = '_'), paste(ycolname, '07', sep = '_'), paste(ycolname, '08', sep = '_'),
                         paste(ycolname, '09', sep = '_'), paste(ycolname, '10', sep = '_'), paste(ycolname, '11', sep = '_'), paste(ycolname, '12', sep = '_'), 
                         paste(ycolname, '13', sep = '_'), paste(ycolname, '14', sep = '_'), paste(ycolname, '15', sep = '_'), paste(ycolname, '16', sep = '_'), 
                         paste(ycolname, '17', sep = '_'), paste(ycolname, '18', sep = '_'), paste(ycolname, '19', sep = '_'))
        q$StationId <- rownames(q) # Create a new column of the StationId
        q <- q %>%
          relocate(StationId)
        expanded[row, 28:46] <- q[, 2:length(q)] # Insert the summed streamflow data into the correct cells of the expanded dataframe
        rm(row, q)
      } 
    }
    
    # Add the third variable's data to the metadata dataframe
    if (nchar(colnames(z[2])) == 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- z[names(z) == huc8] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 47:66] <- NA
          bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(zcolname, '01', sep = '_'), paste(zcolname, '02', sep = '_'), paste(zcolname, '03', sep = '_'), paste(zcolname, '04', sep = '_'),
                         paste(zcolname, '05', sep = '_'), paste(zcolname, '06', sep = '_'), paste(zcolname, '07', sep = '_'), paste(zcolname, '08', sep = '_'),
                         paste(zcolname, '09', sep = '_'), paste(zcolname, '10', sep = '_'), paste(zcolname, '11', sep = '_'), paste(zcolname, '12', sep = '_'), 
                         paste(zcolname, '13', sep = '_'), paste(zcolname, '14', sep = '_'), paste(zcolname, '15', sep = '_'), paste(zcolname, '16', sep = '_'), 
                         paste(zcolname, '17', sep = '_'), paste(zcolname, '18', sep = '_'), paste(zcolname, '19', sep = '_'))
        q$HUC8 <- rownames(q) # Create a new column of the HUC8
        q <- q %>%
          relocate(HUC8)
        expanded[row, 47:66] <- q[, 2:length(q)] # Insert the second variable's data into the correct cells of the expanded dataframe
        rm(row, huc8, q)
      } 
    }
    
    if (nchar(colnames(z[2])) < 8) {
      for (k in expanded$StationId) {
        row <- which(expanded$StationId == k)
        #huc8 <- as.character(expanded$HUC8[row]) # Get the HUC8 for the SNOTEL station in question
        q <- z[names(z) == k] # Grab only the streamflow data for this HUC8
        if (ncol(q) == 0) {
          expanded[row, 47:66] <- NA
          #bad_huc8 <- append(bad_huc8, huc8)
          rm(row, huc8)
          next
        }
        q <- as.data.frame(t(q)) # Transpose the data
        colnames(q) <- c(paste(zcolname, '01', sep = '_'), paste(zcolname, '02', sep = '_'), paste(zcolname, '03', sep = '_'), paste(zcolname, '04', sep = '_'),
                         paste(zcolname, '05', sep = '_'), paste(zcolname, '06', sep = '_'), paste(zcolname, '07', sep = '_'), paste(zcolname, '08', sep = '_'),
                         paste(zcolname, '09', sep = '_'), paste(zcolname, '10', sep = '_'), paste(zcolname, '11', sep = '_'), paste(zcolname, '12', sep = '_'), 
                         paste(zcolname, '13', sep = '_'), paste(zcolname, '14', sep = '_'), paste(zcolname, '15', sep = '_'), paste(zcolname, '16', sep = '_'), 
                         paste(zcolname, '17', sep = '_'), paste(zcolname, '18', sep = '_'), paste(zcolname, '19', sep = '_'))
        q$StationId <- rownames(q) # Create a new column of the StationId
        q <- q %>%
          relocate(StationId)
        expanded[row, 47:66] <- q[, 2:length(q)] # Insert the third variable's data data into the correct cells of the expanded dataframe
        rm(row, q)
      } 
    }
    
    rm(new, f, i, j, k)
    return(expanded)
  }
}

################ Average the SNOTEL Sites By HUC2 By Water Year ######################
mean_swe_HUC2 <- function(great_swe, great_metadata) {
  # First Create a Mean SWE df
  mean_swe <- great_swe %>%
    dplyr::select(-c(3,4,5)) %>%
    group_by(Water_Year) %>%
    summarise(across(c('1000' : '990'), ~mean(.x, na.rm = TRUE))) %>%
    mutate_if(is.numeric, round, digits = 2)
  mean_swe <- mean_swe %>%
    mutate_all(~ifelse(is.nan(.), NA, .))
  # Second group the stations by HUC2, average across rows, and then bind again. 
  nine <- great_metadata[great_metadata$HUC2 == '09' ,]
  nine <- nine$StationId
  ten <- great_metadata[great_metadata$HUC2 == '10' ,]
  ten <- ten$StationId
  eleven <- great_metadata[great_metadata$HUC2 == '11' ,]
  eleven <- eleven$StationId
  thirteen <- great_metadata[great_metadata$HUC2 == '13' ,]
  thirteen <- thirteen$StationId
  fourteen <- great_metadata[great_metadata$HUC2 == '13' ,]
  fourteen <- fourteen$StationId
  fifteen <- great_metadata[great_metadata$HUC2 == '15' ,]
  fifteen <- fifteen$StationId
  sixteen <- great_metadata[great_metadata$HUC2 == '16' ,]
  sixteen <- sixteen$StationId
  seventeen <- great_metadata[great_metadata$HUC2 == '17' ,]
  seventeen <- seventeen$StationId
  eighteen <- great_metadata[great_metadata$HUC2 == '18' ,]
  eighteen <- eighteen$StationId
  
  swe_09 <- mean_swe[names(mean_swe) %in% nine]
  swe_10 <- mean_swe[names(mean_swe) %in% ten]
  swe_11 <- mean_swe[names(mean_swe) %in% eleven]
  swe_13 <- mean_swe[names(mean_swe) %in% thirteen]
  swe_14 <- mean_swe[names(mean_swe) %in% fourteen]
  swe_15 <- mean_swe[names(mean_swe) %in% fifteen]
  swe_16 <- mean_swe[names(mean_swe) %in% sixteen]
  swe_17 <- mean_swe[names(mean_swe) %in% seventeen]
  swe_18 <- mean_swe[names(mean_swe) %in% eighteen]
  
  # Average Across Rows in Each
  swe_09 <- as.data.frame(rowMeans(swe_09[, 1:length(swe_09)], na.rm = TRUE))
  colnames(swe_09) <- 'HUC09'
  swe_10 <- as.data.frame(rowMeans(swe_10[, 1:length(swe_10)], na.rm = TRUE))
  colnames(swe_10) <- 'HUC10'
  swe_11 <- as.data.frame(rowMeans(swe_11[, 1:length(swe_11)], na.rm = TRUE))
  colnames(swe_11) <- 'HUC11'
  swe_13 <- as.data.frame(rowMeans(swe_13[, 1:length(swe_13)], na.rm = TRUE))
  colnames(swe_13) <- 'HUC13'
  swe_14 <- as.data.frame(rowMeans(swe_14[, 1:length(swe_14)], na.rm = TRUE))
  colnames(swe_14) <- 'HUC14'
  swe_15 <- as.data.frame(rowMeans(swe_15[, 1:length(swe_15)], na.rm = TRUE))
  colnames(swe_15) <- 'HUC15'
  swe_16 <- as.data.frame(rowMeans(swe_16[, 1:length(swe_16)], na.rm = TRUE))
  colnames(swe_16) <- 'HUC16'
  swe_17 <- as.data.frame(rowMeans(swe_17[, 1:length(swe_17)], na.rm = TRUE))
  colnames(swe_17) <- 'HUC17'
  swe_18 <- as.data.frame(rowMeans(swe_18[, 1:length(swe_18)], na.rm = TRUE))
  colnames(swe_18) <- 'HUC18'
  
  swe_HUC2 <- bind_cols(as.data.frame(mean_swe[, 1]), swe_09, swe_10, swe_11, swe_13, swe_14, swe_15, swe_16, swe_17, swe_18)
  swe_HUC2 <- swe_HUC2 %>%
    mutate_if(is.numeric, round, digits = 2)
  rm(nine, ten, eleven, thirteen, fourteen, fifteen, sixteen, seventeen, eighteen)
  rm(swe_09, swe_10, swe_11, swe_13, swe_14, swe_15, swe_16, swe_17, swe_18)
  return(swe_HUC2)
}

####################### Gather the Model Metadata ################################
model_metadata <- function(noah_files, metadata) {
  # Create dataframe of Model metadata (bring in equivalent elevation, Snotel ID, and HUC based off lat/lon)
  model_metadata <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(model_metadata) <- c('Latitude', 'Longitude', 'Snotel_Station', 'Elevation', 'HUC8')
  
  # Run through each file name
  for (i in noah_files) {
    local <- data.frame(matrix(nrow = 1, ncol = 5))
    colnames(local) <- c('Latitude', 'Longitude', 'Snotel_Station', 'Elevation', 'HUC8')
    parsed <- str_split_1(i, '_') # Parse the file name
    parsed2 <- str_split_1(parsed[3], '.txt') # Parse the longitude from the '.txt'
    
    if (parsed[2] %in% metadata$Latitude) {
      snotel <- metadata[metadata$Latitude == as.numeric(parsed[2]), ] # Grab the snotel metadata
    local$Latitude <- parsed[2]
    local$Longitude <- parsed2[1]
    local$Snotel_Station <- snotel$StationId
    local$Elevation <- snotel$Elevation
    local$HUC8 <- snotel$HUC8
    
    model_metadata <- rbind(model_metadata, local) # Bind to the overall df
    rm(local, parsed, parsed2, snotel)
    }
  
  }
  rm(snotel, local, parsed, parsed2)
  return(model_metadata)
}


################### Creates Model Data frames from the text files #################
model_data <- function(noah_files, model_meta, param) {
  model <- data.frame(matrix(nrow = 6939, ncol = 0))
  # Run through model file
  for (i in noah_files) {
    data <- read.table(i, sep = ' ')
    data <- data %>%
      mutate(Date = make_date(V1, V2, V3)) %>% # Create a date column
      relocate(Date) %>%
      filter(Date > '2000-09-30') %>% # Remove unnecessary dates
      filter(Date < '2019-10-01') %>%
      mutate(Water_Year = wtr_yr(Date)) %>% # Calls the water year function
      relocate('Water_Year')
    # Parse the name of the file
    parsed <- str_split_1(i, '_')
    if (parsed[2] %in% model_meta$Latitude) {
      station <- model_meta[model_meta$Latitude == as.numeric(parsed[2]), ] # Get corresponding row from metadata (based on the latitude from the file name)
      # Get the Snotel Station Id
      id <- as.character(station$Snotel_Station)
      # Select the desired parameter (param) from the read-in data
      local <- data %>%
        dplyr::select(all_of(param)) %>%
        rename(!!id := param)
      # Bind the resulting one row to the overall dataframe
      model <- cbind(model, local)
      # Create a data frame of dates to add at the end of the loop
      date_cols <- as.data.frame(data[,1:5])
      colnames(date_cols) <- c('Water_Year', 'Date', 'Year', 'Month', 'Day')
      # Remove stuff for the next loop
      rm(data, parsed, station, id, local)
    }
  }
  model <- cbind(date_cols, model) # Bind the date columns to the front
  return(model)
}


#################### Adds Dataframes by Matching Columns #######################
# Sums matching (named) numeric columns of 2 data frames and returns a summed data frame
sum_df <- function(df1, df2, colskip) {
  # Make sure both df's have same # of variables
  
  names_df1 <- names(df1) # Create list of columns names for df1
  names_df1 <- names_df1[-(1:colskip)] # Remove the specified number of beginning columns from that list
  names_df2 <- names(df2)
  names_df2 <- names_df2[-(1:colskip)]
  
  # Create blank dataframe to return
  summed_df <- as.data.frame(df1[, 1:colskip])
  
  for (i in names_df1) {
    local <- as.data.frame(df1[, names(df1) == i] + df2[, names(df2) == i]) # Add the columns from each data frame with matching column-names
    local <- local %>%
      rename(!!i := names(local)[1])
    summed_df <- cbind(summed_df, local) # Bind to the summed dataframe
    rm(local)
  }
  return(summed_df)
}

############# Formula to Calculate Area of Grid Cell on Sphere #################
cell_area <- function(radius, lon, lat, size) {
  # Radius in meters, lon/lat in , size of grid box in degrees
  longitude <- lon * (pi/180) # Convert decimal degrees to radians, keeping the sign convention
  latitude <- lat * (pi/180) 
  grid <- size * (pi/180) # Convert the size of the grid box to radians
  
  lon1 <- longitude - (grid / 2) # Get the lower longitude bound of the grid box
  lon2 <- longitude + (grid / 2) # Get the upper longitude bound of the grid box
  lat1 <- latitude - (grid / 2)
  lat2 <- latitude + (grid / 2)
  
  area = (sin(lat2) - sin(lat1)) * abs(lon2 - lon1) * radius^2
  return(area)
}


################ Create a Dataframe of Summed Precipitation ####################

sum_precipitation <- function(a, b, c) {
  sum <- as.data.frame(c$Water_Year)
  colnames(sum) <- c('Water_Year')
  sum$Water_Year <- as.factor(sum$Water_Year)
  for (i in 6:length(a)) {
    col <- a %>% # Select the water year, date, and ith column 
      dplyr::select(c(1,2, colnames(a)[i]))
    col1 <- col %>% # Filter for the first water year, and for dates before the peak swe date from the b dataframe
      filter(Water_Year == 2001) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 4])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col2 <- col %>%
      filter(Water_Year == 2002) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 5])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col3 <- col %>%
      filter(Water_Year == 2003) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 6])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col4 <- col %>%
      filter(Water_Year == 2004) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 7])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col5 <- col %>%
      filter(Water_Year == 2005) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 8])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col6 <- col %>%
      filter(Water_Year == 2006) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 9])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col7 <- col %>%
      filter(Water_Year == 2007) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 10])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col8 <- col %>%
      filter(Water_Year == 2008) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 11])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col9 <- col %>%
      filter(Water_Year == 2009) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 12])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col10 <- col %>%
      filter(Water_Year == 2010) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 13])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col11 <- col %>%
      filter(Water_Year == 2011) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 14])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col12 <- col %>%
      filter(Water_Year == 2012) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 15])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col13 <- col %>%
      filter(Water_Year == 2013) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 16])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col14 <- col %>%
      filter(Water_Year == 2014) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 17])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col15 <- col %>%
      filter(Water_Year == 2015) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 18])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col16 <- col %>%
      filter(Water_Year == 2016) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 19])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col17 <- col %>%
      filter(Water_Year == 2017) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 20])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col18 <- col %>%
      filter(Water_Year == 2018) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 21])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    col19 <- col %>%
      filter(Water_Year == 2019) %>%
      filter(Date <= as.Date(b[b$Station == as.character(colnames(a)[i]), 22])) %>%
      summarise(across(c(colnames(a)[i]), ~sum(.x, na.rm = TRUE)))
    
    cols <- rbind(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19)
    rm(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19)
    sum <- cbind(sum, cols)
  }
  return(sum)
}

#################### Average Temperature for Each Snotel Station ###############
mean_temperature <- function(all_temp, all_peaks) {
  mean <- as.data.frame(sum_streamflow$Water_Year)
  colnames(mean) <- c('Water_Year')
  for (i in 6:length(all_temp)) {
    col <- all_temp %>%
      dplyr::select(c(1,2, colnames(all_temp)[i]))
    col1 <- col %>%
      filter(Water_Year == 2001) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 4])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col2 <- col %>%
      filter(Water_Year == 2002) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 5])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col3 <- col %>%
      filter(Water_Year == 2003) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 6])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col4 <- col %>%
      filter(Water_Year == 2004) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 7])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col5 <- col %>%
      filter(Water_Year == 2005) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 8])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col6 <- col %>%
      filter(Water_Year == 2006) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 9])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col7 <- col %>%
      filter(Water_Year == 2007) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 10])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col8 <- col %>%
      filter(Water_Year == 2008) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 11])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col9 <- col %>%
      filter(Water_Year == 2009) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 12])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col10 <- col %>%
      filter(Water_Year == 2010) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 13])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col11 <- col %>%
      filter(Water_Year == 2011) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 14])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col12 <- col %>%
      filter(Water_Year == 2012) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 15])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col13 <- col %>%
      filter(Water_Year == 2013) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 16])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col14 <- col %>%
      filter(Water_Year == 2014) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 17])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col15 <- col %>%
      filter(Water_Year == 2015) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 18])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col16 <- col %>%
      filter(Water_Year == 2016) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 19])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col17 <- col %>%
      filter(Water_Year == 2017) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 20])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col18 <- col %>%
      filter(Water_Year == 2018) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 21])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    col19 <- col %>%
      filter(Water_Year == 2019) %>%
      filter(Date <= as.Date(all_peaks[all_peaks$Station == as.character(colnames(all_temp)[i]), 22])) %>%
      summarise(across(c(colnames(all_temp)[i]), ~mean(.x, na.rm = TRUE)))
    
    cols <- rbind(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19)
    rm(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19)
    mean <- cbind(mean, cols)
  }
  mean <- mean %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) # Remove NaNs that come from there being no data to take the average of before the peak date
  mean$Water_Year <- as.factor(mean$Water_Year)
  return(mean)
}


############## Returns a Contingency Table based on Two Dataframes #############

contingency_table <- function(obs, model, metadata, model_meta, huc8) {
  ####### Observations
  if (huc8 %in% metadata$HUC8) {   
      # Filter for the selected stations
      stations <- metadata %>%
        filter(HUC8 == as.character(huc8))
      
      # Create list of Unique ids
      ids <- as.character(stations$StationId)
      
      # Get the three smallest years for OBSERVATIONS
      # Depends on the type of dataframe supplied (Hucs or ids as column names)
      if (nchar(colnames(obs)[2]) == 8) { # If the column names are huc8s
        obs_new <- data.frame(obs[, names(obs) == huc8])
        if (is_empty(obs_new)) {
          contingency <- data.frame(matrix(ncol = 5, nrow = 1))
          colnames(contingency) <- c('HUC8', 'Hits', 'FAs', 'Misses', 'CNs') # Where FA = False Alarm and CN = Correct Negative
          contingency$HUC8 <- as.factor(huc8)
          contingency$Hits <- NA # Hit: obs has the year, so does the model
          contingency$FAs <- NA # False Alarm: Year in the model that isnt in the obs
          contingency$Misses <- NA # Miss: Obs has the year, the model doesnt
          contingency$CNs <- NA # Correct Negative: Number of years that arent in either
          return(contingency)
        }
        obs_new <- cbind(obs$Water_Year, obs_new)
        obs_new <- obs_new %>%
          rename(Water_Year = `obs$Water_Year`) %>%
          arrange(pick(2))
        # %>% # Sort the column by ascending
        # filter(!is.na) # Get rid of NA observations
        ordered <- as.character(obs_new$Water_Year) # Save the ordered years to a character vector
        least_obs <- ordered[1:3] # Pick the three smallest years
        
      }
      if (nchar(colnames(obs)[2]) < 5) { # If the column names are ids
        obs_new <- data.frame(obs[, names(obs) %in% ids])
        colnames(obs_new) <- ids
        obs_new <- cbind(obs$Water_Year, obs_new)
        obs_new <- obs_new %>%
          rename(Water_Year = `obs$Water_Year`) 
        if (length(obs_new) > 2) {
          obs_new <- obs_new %>%
            mutate(HUC8 = rowMeans(.[, 2:ncol(obs_new)], na.rm = TRUE)) %>%
            rename(!!huc8 := HUC8)
          obs_new <- obs_new %>%
            dplyr::select(c(1, ncol(obs_new))) %>%
            arrange(pick(2))
        } else {
          obs_new <- obs_new %>%
            arrange(pick(2)) %>%
            rename(!!huc8 := colnames(obs_new)[2])
        }
        ordered <- as.character(obs_new$Water_Year) # Save the ordered years to a character vector
        least_obs <- ordered[1:3] # Pick the three smallest years
      }
      rm(stations)
  } else { 
    contingency <- data.frame(matrix(ncol = 5, nrow = 1))
    colnames(contingency) <- c('HUC8', 'Hits', 'FAs', 'Misses', 'CNs') # Where FA = False Alarm and CN = Correct Negative
    test <- as.integer(least_model %in% least_obs) # Returns 1 = True, 0 = False
    contingency$HUC8 <- as.factor(huc8)
    contingency$Hits <- NA # Hit: obs has the year, so does the model
    contingency$FAs <- NA # False Alarm: Year in the model that isnt in the obs
    contingency$Misses <- NA # Miss: Obs has the year, the model doesnt
    contingency$CNs <- NA # Correct Negative: Number of years that arent in either
    return(contingency)
  }
  
  ######## Model
  if (huc8 %in% model_meta$HUC8) {
      # Filter for the selected stations
      stations <- model_meta %>%
        filter(HUC8 == as.character(huc8))
      # Create list of Unique ids
      ids <- as.character(stations$Snotel_Station)
      # Get the three smallest years for MODEL
      if (nchar(colnames(model)[2]) == 8) { # If the column names are huc8s
        model_new <- model[, names(model) == huc8]
        if (is_empty(model_new)) {
          contingency <- data.frame(matrix(ncol = 5, nrow = 1))
          colnames(contingency) <- c('HUC8', 'Hits', 'FAs', 'Misses', 'CNs') # Where FA = False Alarm and CN = Correct Negative
          contingency$HUC8 <- as.factor(huc8)
          contingency$Hits <- NA # Hit: obs has the year, so does the model
          contingency$FAs <- NA # False Alarm: Year in the model that isnt in the obs
          contingency$Misses <- NA # Miss: Obs has the year, the model doesnt
          contingency$CNs <- NA # Correct Negative: Number of years that arent in either
          return(contingency)
        }
        model_new <- cbind(model$Water_Year, model_new)
        model_new <- model_new %>%
          rename(Water_Year = `model$Water_Year`) %>%
          arrange(pick(2))
        ordered <- as.character(model_new$Water_Year) # Save the ordered years to a character vector
        least_model <- ordered[1:3] # Pick the three smallest years
        
      }
      if (nchar(colnames(model)[2]) < 5) { # If the column names are ids
        model_new <- data.frame(model[, names(model) %in% ids])
        colnames(model_new) <- ids
        model_new <- cbind(model$Water_Year, model_new)
        model_new <- model_new %>%
          rename(Water_Year = `model$Water_Year`) 
        if (length(model_new) > 2) {
          model_new <- model_new %>%
            mutate(HUC8 = rowMeans(.[, 2:ncol(model_new)], na.rm = TRUE)) %>%
            rename(!!huc8 := HUC8)
          model_new <- model_new %>%
            dplyr::select(c(1, ncol(model_new))) %>%
            arrange(pick(2))
        } else {
          model_new <- model_new %>%
            arrange(pick(2)) %>%
            rename(!!huc8 := colnames(model_new)[2])
        }
        ordered <- as.character(model_new$Water_Year) # Save the ordered years to a character vector
        least_model <- ordered[1:3] # Pick the three smallest years
      }
      rm(stations)
  } else { 
    contingency <- data.frame(matrix(ncol = 5, nrow = 1))
    colnames(contingency) <- c('HUC8', 'Hits', 'FAs', 'Misses', 'CNs') # Where FA = False Alarm and CN = Correct Negative
    test <- as.integer(least_model %in% least_obs) # Returns 1 = True, 0 = False
    contingency$HUC8 <- as.factor(huc8)
    contingency$Hits <- NA # Hit: obs has the year, so does the model
    contingency$FAs <- NA # False Alarm: Year in the model that isnt in the obs
    contingency$Misses <- NA # Miss: Obs has the year, the model doesnt
    contingency$CNs <- NA # Correct Negative: Number of years that arent in either
    return(contingency)
  }
  
  # We have the list for obs and model: now make our contingency table
  contingency <- data.frame(matrix(ncol = 5, nrow = 1))
  colnames(contingency) <- c('HUC8', 'Hits', 'FAs', 'Misses', 'CNs') # Where FA = False Alarm and CN = Correct Negative
  
  test <- as.integer(least_model %in% least_obs) # Returns 1 = True, 0 = False
  contingency$HUC8 <- as.factor(huc8)
  
  if (sum(test == 1) == 0) { # three 0s
    contingency$Hits <- 0 # Hit: obs has the year, so does the model
    contingency$FAs <- 3 # False Alarm: Year in the model that isnt in the obs
    contingency$Misses <- 3 # Miss: Obs has the year, the model doesnt
    contingency$CNs <- 13 # Correct Negative: Number of years that arent in either
  }
  if (sum(test == 1) == 1) { # 1 one, 2 Os
    contingency$Hits <- 1
    contingency$FAs <- 2
    contingency$Misses <- 2
    contingency$CNs <- 14
  }
  if (sum(test == 1) == 2) { # 2 ones, 1 zero
    contingency$Hits <- 2
    contingency$FAs <- 1
    contingency$Misses <- 1
    contingency$CNs <- 15
  }
  if (sum(test == 1) == 3) { # 3 ones, 0 zeros
    contingency$Hits <- 3
    contingency$FAs <- 0
    contingency$Misses <- 0
    contingency$CNs <- 16
  }
  
  return(contingency)
}




