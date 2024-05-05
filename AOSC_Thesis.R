############################ AOSC Senior Thesis Code ##############################
# This script is written entirely by Ethan Heidtman, starting in the Fall of 2023

### Table of Contents ###
# Functions and Packages
# Importing Metadata
# Importing SNOTEL Data (SWE, Air Temperature)
# Importing PDO Index
# Calculating Pertinent Hydrograph Characteristic Time Series
# Linear Inverse Modeling


############################# Functions and Packages ###########################

          path1 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/All_SWE_Data/'
          # Where the Streamflow Data Is Stored
          path2 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/Streamflow_Data/'
          # Where the Air Temperature Data is Stored
          path3 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/AirTemp/'
          # Where the Precipitation Data is Stored
          path4 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/Precip/'
          # Where any figures will be saved to
          path5 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/Figures/'
          # Where the PDO index is stored
          path6 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/My Drive/Senior Thesis/Indices/'
          # Where the HUC2 Shapefiles are located
          path7 <- '/Users/ethanheidtman/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/.shortcut-targets-by-id/1GX8mj4pqw1SoB9JTetKChQp2Y0mgEAss/Ethan_Heidtman_Summer_2023/HUC2_Shapefiles'
          # Names of the shapefiles
          shapes_path <- c('WBD_09_HU2_Shape/Shape/', 'WBD_10_HU2_Shape/Shape/', 'WBD_11_HU2_Shape/Shape/', 'WBD_13_HU2_Shape/Shape/', 'WBD_14_HU2_Shape/Shape/',
                                    'WBD_15_HU2_Shape/Shape/', 'WBD_16_HU2_Shape/Shape/', 'WBD_17_HU2_Shape/Shape/', 'WBD_18_HU2_Shape/Shape/')
          
          setwd(path1)
          
          # Source our functions 
          source('~/Library/CloudStorage/GoogleDrive-eheidtma@terpmail.umd.edu/.shortcut-targets-by-id/1GX8mj4pqw1SoB9JTetKChQp2Y0mgEAss/Ethan_Heidtman_Summer_2023/SWE/Functions.R')
          
          # Load lots of packages 
          library(easypackages)
          libraries('tidyverse', 'sf', 'ggplot2', 'dplyr', 'lubridate', 'dataRetrieval', 'data.table', 'stringr', 'usmap', 'broom', 'ggthemes', 'viridis',
                    'maptools', 'maps', 'mapdata', 'anytime', 'reshape2', 'scales', 'ggpubr', 'ggbreak', 'gganimate', 'raster', 'rgdal', 'ncdf4',
                    'ggridges', 'grwat', 'WaveletComp', 'multispatialCCM', 'rEDM', 'trend', 'FlowScreen', 'zoo', 'changepoint', 'mcp', 'segmented',
                    'strucchange', 'ggpmisc', 'moments', 'Kendall', 'limSolve', 'LIM', 'ggtext', 'RColorBrewer', 'cowplot', 'mblm', 'htmlTable')
############################ Metadata and MISC. ################################

          # Map of the US 
          us <- map_data('state')
          us <- subset(us, region %in% c("california", "oregon", 'nevada', 'colorado', 'washington', 'wyoming', 'idaho', 'new mexico', 'montana', 'texas',
                                         'arizona', 'north dakota', 'south dakota', 'nebraska', 'kansas', 'oklahoma'))
          
          
          # Read in the Full SNOTEL SWE Metadata: 897 Stations
          swe_meta <- read.csv('SNOTEL.csv', skip = 916)
          swe_meta$Elevation <- swe_meta$Elevation / 3.281 # Convert from feet to meters
          swe_meta$HUC8[465] <- '09040001' # Add an extra 0 to the front of this HUC8
          swe_meta <- swe_meta %>%
            mutate_at('Elevation', ~ round(., 2)) %>% # Round the elevation values to 2 digits
            dplyr::select(-c('EndDate', 'StartDate')) %>% # Get rid of these two columns
            mutate(HUC2 = substr(HUC8, 1, 2)) %>% # Create a HUC2 column
            relocate(HUC2, .after = HUC8)
          
          # Gather just the Alaskan Stations so they can be removed later
          alaska <- swe_meta[swe_meta$StateName == 'ALASKA', ] # Alaskan Stations 
          
          for (i in 1:length(shapes_path)) {
            new_path <- paste(path7, shapes_path[i], sep = '/')
            setwd(new_path)
            if (i == 1) { huc9 <- st_read('WBDHU2.shp') }
            else if (i == 2) { huc10 <- st_read('WBDHU2.shp') }
            else if (i == 3) { huc11 <- st_read('WBDHU2.shp') }
            else if (i == 4) { huc13 <- st_read('WBDHU2.shp') }
            else if (i == 5) { huc14 <- st_read('WBDHU2.shp') }
            else if (i == 6) { huc15 <- st_read('WBDHU2.shp') }
            else if (i == 7) { huc16 <- st_read('WBDHU2.shp') }
            else if (i == 8) { huc17 <- st_read('WBDHU2.shp') }
            else if (i == 9) { huc18 <- st_read('WBDHU2.shp') }
            rm(new_path)
          }
          rm(i)
          labels <- data.frame(longitude = c(-97, -103, -99, -106, -109, -112, -116, -119, -119.5),
                               latitude = c(48, 44, 36, 34, 39.5, 33.5, 40, 47, 36),
                               label = c('09', '10', '11', '13', '14', '15', '16', '17', '18'))
############################# SNOTEL Data ######################################

          ####
              # Read all the SWE Data
              setwd(path1)
              swe <- read_snotel(path1)
              colnames(swe)[2:length(swe)] <- gsub("[^0-9-]{3,4}", "", colnames(swe)[2:length(swe)]) # Make column names just the station ID
              colnames(swe)[2:length(swe)] <- str_extract(colnames(swe)[2:length(swe)], "[[:digit:]]+(?!.)")
              
              # Remove station 1317 from SWE so all columns match
              swe <- swe %>%
                dplyr::select(-c('1317'))
              
              # Tidy the swe dataset
              swe$Year <- year(ymd(swe$Date))
              swe$Month <- month(ymd(swe$Date))
              swe$Day <- day(ymd(swe$Date))
              swe <- swe %>%
                relocate(c('Year', 'Month', 'Day')) %>%
                relocate('Date') %>%
                mutate(Water_Year = wtr_yr(Date)) %>% # Calls the water year function
                relocate('Water_Year') %>% 
                mutate(across(c('301' : '1033'), ~ .x * 25.4)) %>% # Inches to millimeters
                mutate_if(is.numeric, round, digits = 1) # Round to 2 sig figs
              
              swe <- swe %>%
                dplyr::select(c(1:5, colnames(swe[colnames(swe)[6:length(swe)] %in% swe_meta$StationId])))
              
              # Trim the data to our time frame
              swe <- swe %>%
                filter(Date > '1979-09-30') %>%
                filter(Date < '2022-10-01')
              
              # Remove the Alaskan stations from SWE and from the metadata
              swe <- swe %>%
                dplyr::select(-c(colnames(swe[colnames(swe) %in% alaska$StationId]))) 
          ####
          
          #### Gather the Air Temperature Data
              setwd(path3)
              
              one <- read.csv('80to90.csv', skip = 952)
              two <- read.csv('91to01.csv', skip = 952)
              three <- read.csv('02to11.csv', skip = 952)
              four <- read.csv('12to22.csv', skip = 952)
              temp <- rbind(one, two, three, four)
              rm(one, two, three, four)
              
              colnames(temp)[2:length(temp)] <- gsub("[^0-9-]{3,4}", "", colnames(temp)[2:length(temp)]) # Make column names just the station ID
              colnames(temp)[2:length(temp)] <- str_extract(colnames(temp)[2:length(temp)], "[[:digit:]]+(?=s)")
              
              temp <- temp %>% # Remove alaskan stations and make sure temp and swe have the same exact stations
                dplyr::select(-c(colnames(temp[colnames(temp) %in% alaska$StationId]))) 
              temp <- temp %>%
                dplyr::select(-c(setdiff(colnames(temp)[2:length(temp)], swe_meta$StationId)))
              
              # Tidy the temperature data and change its form a bit
              temp$Year <- year(ymd(temp$Date))
              temp$Month <- month(ymd(temp$Date))
              temp$Day <- day(ymd(temp$Date))
              temp <- temp %>%
                relocate(c('Year', 'Month', 'Day')) %>%
                relocate('Date') %>%
                mutate(Water_Year = wtr_yr(Date)) %>% # Calls the water year function
                relocate('Water_Year') %>% 
                mutate(across(c('301' : '1033'), ~ ((.x - 32)*(5/9)))) %>% # Fahrenheit to celsius
                mutate_if(is.numeric, round, digits = 2) %>% # Round to 2 sig figs
                filter(Date > '1979-09-30') %>%
                filter(Date < '2022-10-01')
              
              

###################### Pacific Decadal Oscillation Index #######################
          # https://www.ncei.noaa.gov/access/monitoring/pdo/
            # https://psl.noaa.gov/pdo/
          # https://psl.noaa.gov/enso/dashboard.html
          # https://bookdown.org/rdpeng/exdata/the-base-plotting-system-1.html
          
          setwd(path6)
          
          # Read in the PDO index
              pdo <- read.csv('PDO_ERSST_V5.csv')
              colnames(pdo)[2] <- 'PDO'
              pdo$Date <- as.Date(pdo$Date)
              pdo <- pdo %>%
                filter(Date > '1979-09-01') %>%
                filter(Date < '2022-10-01') %>%
                mutate(Water_Year = wtr_yr(Date)) %>% # Calls the water year function
                relocate('Water_Year')
              pdo$Water_Year <- as.factor(pdo$Water_Year)
              # pdo$Year <- year(ymd(pdo$Date))
              pdo$Month <- month(ymd(pdo$Date))
              pdo$Day <- day(ymd(pdo$Date))
              pdo <- pdo %>%
                relocate(c('Month', 'Day'), .after = Date)
              
              mean_pdo <- pdo %>%
                group_by(Water_Year) %>%
                summarise_if(is.numeric, mean) %>%
                mutate_if(is.numeric, round, digits = 2)

          ## Plot the PDO over the course of the time series
              
              plot <- ggplot(pdo, aes(x = Date, y = PDO)) + 
                geom_line() + 
                labs(x = 'Date', y = 'PDO Index', title = 'NCEI ERSST V5 Monthly PDO Index: 1980-2022', 
                     caption = '**Figure 4** The raw NCEI ERSST V5 PDO index, where positive values are warm phase, and negative values are the cold phase. 43 years of<br>monthly data are shown. Note that the PDO has consistently been decreasing since the mid 2010s.') +
                scale_x_date(breaks = '2 years', date_labels = '%Y') + 
                theme_clean() + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      legend.key.size = unit(1, 'cm'),
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 12, vjust = 0.9, hjust = 1, angle = 30),
                      axis.text.y = element_text(size = 14))+ 
                theme(aspect.ratio = 2/10) 
              # ggsave(filename = 'pdo_line.png', width = 13, height = 5, plot, path = path5)
                
             

############# Calculating Time Series of Characteristics #######################

          # 1. Filter SNOTEL Stations that have no more than 1000 missing observations for SWE and Air Temperature
              na <- swe %>%
                summarise(across(c('301' : '1033'),  ~ sum(is.na(.x)))) %>%
                dplyr::select_if(~ any(. < 1000))
              swe <- swe %>%
                dplyr::select(c(1:5), colnames(na))
              swe_meta <- swe_meta %>%
                filter(StationId %in% colnames(swe))
                
              na <- temp %>%
                summarise(across(c('301' : '1033'),  ~ sum(is.na(.x)))) %>%
                dplyr::select_if(~ any(. < 1000))
              temp <- temp %>%
                dplyr::select(c(1:5, colnames(na)))
              swe_meta <- swe_meta %>%
                filter(StationId %in% colnames(temp))
              swe <- swe %>%
                dplyr::select(c(1:5, which(colnames(swe) %in% swe_meta$StationId)))
              temp <- temp %>%
                dplyr::select(which(colnames(temp) %in% colnames(swe)))
          
              # Remove NA dataframe
              rm(na, alaska)
              
              # Remove station 613 in HUC09
              swe_meta <- swe_meta %>%
                filter(!StationId == 613)
              swe <- swe %>%
                dplyr::select(-c('613'))
              temp <- temp %>%
                dplyr::select(-c('613'))
              
          # 2. MAGNITUDE OF PEAK SWE
              peak_swe <- swe %>%
                group_by(Water_Year) %>%
                summarise(across(c('303' : '877'), ~max(., na.rm = TRUE))) # Collect the max swe for each column for each water year
              peak_swe[peak_swe == '-Inf'] <- NA
              
          # 3. DATE OF PEAK SWE
              date_peak_swe <- as.data.frame(peak_swe$Water_Year) # Initialize empty data frame
              colnames(date_peak_swe) <- 'Water_Year'
              for (i in 2:length(peak_swe)) { # For each SNOTEL station
                df <- swe %>%
                  dplyr::select(c(1,2, colnames(swe)[i + 4])) # Select SNOTEL station i
                local <- data.frame(matrix(ncol = 1, nrow = 0))
                colnames(local) <- colnames(peak_swe)[i]
                for (x in date_peak_swe$Water_Year) { # For each year 
                  data <- df %>%
                    filter(Water_Year == x) # Filter for water year x
                  index <- which.max(data[, 3]) # Find the index that is the max value of SWE
                  if (is_empty(index)) { # if there is no max index (i.e., this year is all NA obs for SWE)
                    date <- NA
                    local <- rbind(local, date)
                    next
                  }
                  if (max(data[, 3], na.rm = TRUE) == 0) { # If the max SWE is 0 for this year
                    date <- NA
                    local <- rbind(local, date)
                    next
                  }
                  date <- data$Date[index] # Collect the date at that index
                  local <- rbind(local, as.character(date)) # bind to local data frame
                  rm(data, date, index)
                }
                colnames(local) <- colnames(peak_swe)[i]
                local[, 1] <- as.Date(local[, 1]) # Make class date
                date_peak_swe <- cbind(date_peak_swe, as.data.frame(local)) # bind to product dataframe
                rm(df, local)
              }
              
              # Change from date to # of days from start of water year
              date_peak_swe <- date_peak_swe %>%
                mutate_if(is.Date, as.character)
              for (i in 2:length(date_peak_swe)) {
                for (x in 1:length(date_peak_swe$Water_Year)) {
                  # Compute number of days that have elapsed since the start of the water year for each observation
                  ref <- as.Date(paste0(as.numeric(date_peak_swe$Water_Year[x]) - 1, '-10-01'))
                  diff <- as.integer(difftime(date_peak_swe[x, i], ref))
                  date_peak_swe[x, i] <- diff
                  rm(ref, dif)
                }
              }

              date_peak_swe <- date_peak_swe %>%
                mutate_if(is.character, as.numeric)
              
          # 4. Date of 5% of Peak SWE
              date_5_swe <- as.data.frame(peak_swe$Water_Year) # Initialize empty data frame
              colnames(date_5_swe) <- 'Water_Year'
              for (i in 2:length(peak_swe)) {
                data <- swe %>%
                  dplyr::select(c(1, 2, colnames(peak_swe)[i])) # Select the proper column of SWE data
                local <- data.frame(matrix(ncol = 1, nrow = 0))
                colnames(local) <- colnames(peak_swe)[i]
                for (x in 1:length(date_5_swe$Water_Year)) {
                  df <- data %>%
                    filter(Water_Year == peak_swe$Water_Year[x]) # Filter for the right water year
                  peak = as.numeric(peak_swe[x,i]) # Get the peak value
                  if (peak == 0 | is.na(peak)) { # If the peak is 0 or there is no peak
                    local <- rbind(local, NA) # 5% date is NA
                    colnames(local) <- colnames(peak_swe)[i]
                    next
                  }
                  end = peak*.05 # Compute the 5% value
                   
                  df <- df %>%
                    filter(df[, 3] <= end) # Filter for SWE values that are less than the 5% value
                  df <- df %>% # Filter for Dates that are after the date of peak SWE
                    filter(Date > as.Date(paste0(peak_swe$Water_Year[x] - 1, '-10-01')) + as.numeric(date_peak_swe[x, i]))
                  
                  if (nrow(df) == 0) { # If there are no dates after the date of peak SWE
                    date <- as.Date(paste0(peak_swe$Water_Year[x], '-09-30'))  # Assign the end of the water year as the date of 5%
                    local <- rbind(local, as.Date(date))
                    next
                  }
                    
                  date <- df$Date[1]
                  local <- rbind(local, date) # Bind the date of 5%
                  rm(peak, end, date)
                }
                colnames(local) <- colnames(peak_swe)[i]
                date_5_swe <- cbind(date_5_swe, local) # Bind the results for this station to the overall df
                rm(df, local, data)
              }
              date_5_swe[1,108] <- NA
              date_5_swe[4,285] <- NA
              
              # Change from date to # of days from start of water year
              date_5_swe <- date_5_swe %>%
                mutate_if(is.Date, as.character)
              for (i in 2:length(date_5_swe)) {
                for (x in 1:length(date_5_swe$Water_Year)) {
                  # Compute number of days that have elapsed since the start of the water year for each observation
                  ref <- as.Date(paste0(as.numeric(date_5_swe$Water_Year[x]) - 1, '-10-01'))
                  diff <- as.integer(difftime(date_5_swe[x, i], ref)) 
                  date_5_swe[x, i] <- diff
                  rm(ref, dif)
                }
              }
              
              date_5_swe <- date_5_swe %>%
                mutate_if(is.character, as.numeric)
              
          # 5. Melt Time 
              melt_time <- as.data.frame(peak_swe$Water_Year)
              colnames(melt_time) <- 'Water_Year'
              for (i in 2:length(date_peak_swe)) { # For each SNOTEL station 
                  df1 <- as.data.frame(date_peak_swe[,i]) # Gather the date of peak SWE
                  df2 <- as.data.frame(date_5_swe[, i]) # Gather the date of 5% SWE
                  df1 <- cbind(df1, df2) # Bind the two together so they can be subtracted
                  df1[, 1] <- as.numeric(df1[, 1]) # Make both numeric so they can be subtracted
                  df1[, 2] <- as.numeric(df1[, 2])
                  df1 <- df1 %>%
                    mutate(diff = .[[2]] - .[[1]]) %>% # Subtract the columns
                    rename(!!colnames(date_peak_swe)[i] := diff) # Rename to the correct SNOTEL
                melt_time <- cbind(melt_time, df1[3])
                rm(df1, df2)
              }
              
############################################## Linear Inverse Model Stuff ##############################
              
              # Create a matrix of our PDO indices by year and month
              pdo_matrix <- pdo %>%
                pivot_wider(-c(Date, Day), names_from = Month, values_from = PDO) %>%
                mutate_if(is.factor, as.numeric) %>%
                dplyr::select(-1) %>%
                as.matrix()
              rownames(pdo_matrix) <- peak_swe$Water_Year
              
              # Convert the data to z-scores
              # Use the mean and standard deviation of the whole sample
              z_peak_swe <- peak_swe %>%
                mutate(across(c('303' : '877'), ~ ((. - mean(as.matrix(peak_swe[,2:314]), na.rm = TRUE)))/sd(as.matrix(peak_swe[,2:314]), na.rm = TRUE))) %>%
                mutate_if(is.numeric, round, digits = 2)
              
              z_date_peak_swe <- date_peak_swe %>%
                mutate(across(c('303' : '877'), ~ ((. - mean(as.matrix(date_peak_swe[,2:314]), na.rm = TRUE)))/sd(as.matrix(date_peak_swe[,2:314]), na.rm = TRUE))) %>%
                mutate_if(is.numeric, round, digits = 2)
              
              z_date_5_swe <- date_5_swe %>%
                mutate(across(c('303' : '877'), ~ ((. - mean(as.matrix(date_5_swe[,2:314]), na.rm = TRUE)))/sd(as.matrix(date_5_swe[,2:314]), na.rm = TRUE))) %>%
                mutate_if(is.numeric, round, digits = 2)
              
              z_melt_time <- melt_time %>%
                mutate(across(c('303' : '877'), ~ ((. - mean(as.matrix(melt_time[,2:314]), na.rm = TRUE)))/sd(as.matrix(melt_time[,2:314]), na.rm = TRUE))) %>%
                mutate_if(is.numeric, round, digits = 2)
              
            
              # Time Period 1 1980-1989
                  # Peak SWE
                  lim1_peakSWE <- as.data.frame(z_peak_swe$Water_Year[1:10]) # Select the rows corresponding to the time frame
                  colnames(lim1_peakSWE) <- 'Water_Year'
                  for (i in 2:length(z_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_peak_swe, startyear = 1980, endyear = 1989, station = i) # Compute the linear model
                    lim1_peakSWE <- merge(lim1_peakSWE, result, by = 'Water_Year', all = TRUE) # Join the linear model to the main df
                  }
                  lim1_peakSWE <- lim1_peakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>% # Create a Decade column
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('One')) %>% # Create a time period column
                    relocate(Period, .after = Decade)
                 
                  # Date Peak SWE
                  lim1_datePeakSWE <- as.data.frame(z_date_peak_swe$Water_Year[1:10])
                  colnames(lim1_datePeakSWE) <- 'Water_Year'
                  for (i in 2:length(z_date_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_peak_swe, startyear = 1980, endyear = 1989, station = i)
                    lim1_datePeakSWE <- merge(lim1_datePeakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim1_datePeakSWE <- lim1_datePeakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('One')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date 5 SWE
                  lim1_date5swe <- as.data.frame(z_date_5_swe$Water_Year[1:10])
                  colnames(lim1_date5swe) <- 'Water_Year'
                  for (i in 2:length(z_date_5_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_5_swe, startyear = 1980, endyear = 1989, station = i)
                    lim1_date5swe <- merge(lim1_date5swe, result, by = 'Water_Year', all = TRUE)
                  }
                  lim1_date5swe <- lim1_date5swe %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('One')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Melt Time
                  lim1_meltTime <- as.data.frame(z_melt_time$Water_Year[1:10])
                  colnames(lim1_meltTime) <- 'Water_Year'
                  for (i in 2:length(z_melt_time)) {
                    result <- LIM(E = pdo_matrix, F = z_melt_time, startyear = 1980, endyear = 1989, station = i)
                    lim1_meltTime <- merge(lim1_meltTime, result, by = 'Water_Year', all = TRUE)
                  }
                  lim1_meltTime <- lim1_meltTime %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('One')) %>%
                    relocate(Period, .after = Decade)
                  
          
              # Time Period 2 1990-1999
                  # Peak SWE
                  lim2_peakSWE <- as.data.frame(z_peak_swe$Water_Year[11:20])
                  colnames(lim2_peakSWE) <- 'Water_Year'
                  for (i in 2:length(z_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_peak_swe, startyear = 1990, endyear = 1999, station = i)
                    lim2_peakSWE <- merge(lim2_peakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim2_peakSWE <- lim2_peakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Two')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date Peak SWE
                  lim2_datePeakSWE <- as.data.frame(z_date_peak_swe$Water_Year[11:20])
                  colnames(lim2_datePeakSWE) <- 'Water_Year'
                  for (i in 2:length(z_date_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_peak_swe, startyear = 1990, endyear = 1999, station = i)
                    lim2_datePeakSWE <- merge(lim2_datePeakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim2_datePeakSWE <- lim2_datePeakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Two')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date 5 SWE
                  lim2_date5swe <- as.data.frame(z_date_5_swe$Water_Year[11:20])
                  colnames(lim2_date5swe) <- 'Water_Year'
                  for (i in 2:length(z_date_5_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_5_swe, startyear = 1990, endyear = 1999, station = i)
                    lim2_date5swe <- merge(lim2_date5swe, result, by = 'Water_Year', all = TRUE)
                  }
                  lim2_date5swe <- lim2_date5swe %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Two')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Melt Time
                  lim2_meltTime <- as.data.frame(z_melt_time$Water_Year[11:20])
                  colnames(lim2_meltTime) <- 'Water_Year'
                  for (i in 2:length(z_melt_time)) {
                    result <- LIM(E = pdo_matrix, F = z_melt_time, startyear = 1990, endyear = 1999, station = i)
                    lim2_meltTime <- merge(lim2_meltTime, result, by = 'Water_Year', all = TRUE)
                  }
                  lim2_meltTime <- lim2_meltTime %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Two')) %>%
                    relocate(Period, .after = Decade)
              
              # Time Period 3 2000-2009
                  # Peak SWE
                  lim3_peakSWE <- as.data.frame(z_peak_swe$Water_Year[21:30])
                  colnames(lim3_peakSWE) <- 'Water_Year'
                  for (i in 2:length(z_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_peak_swe, startyear = 2000, endyear = 2009, station = i)
                    lim3_peakSWE <- merge(lim3_peakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim3_peakSWE <- lim3_peakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Three')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date Peak SWE
                  lim3_datePeakSWE <- as.data.frame(z_date_peak_swe$Water_Year[21:30])
                  colnames(lim3_datePeakSWE) <- 'Water_Year'
                  for (i in 2:length(z_date_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_peak_swe, startyear = 2000, endyear = 2009, station = i)
                    lim3_datePeakSWE <- merge(lim3_datePeakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim3_datePeakSWE <- lim3_datePeakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Three')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date 5 SWE
                  lim3_date5swe <- as.data.frame(z_date_5_swe$Water_Year[21:30])
                  colnames(lim3_date5swe) <- 'Water_Year'
                  for (i in 2:length(z_date_5_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_5_swe, startyear = 2000, endyear = 2009, station = i)
                    lim3_date5swe <- merge(lim3_date5swe, result, by = 'Water_Year', all = TRUE)
                  }
                  lim3_date5swe <- lim3_date5swe %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Three')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Melt Time
                  lim3_meltTime <- as.data.frame(z_melt_time$Water_Year[21:30])
                  colnames(lim3_meltTime) <- 'Water_Year'
                  for (i in 2:length(z_melt_time)) {
                    result <- LIM(E = pdo_matrix, F = z_melt_time, startyear = 2000, endyear = 2009, station = i)
                    lim3_meltTime <- merge(lim3_meltTime, result, by = 'Water_Year', all = TRUE)
                  }
                  lim3_meltTime <- lim3_meltTime %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Three')) %>%
                    relocate(Period, .after = Decade)
              
              # Time Period 4 2010-2022
                  # Peak SWE
                  lim4_peakSWE <- as.data.frame(z_peak_swe$Water_Year[31:43])
                  colnames(lim4_peakSWE) <- 'Water_Year'
                  for (i in 2:length(z_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_peak_swe, startyear = 2010, endyear = 2022, station = i)
                    lim4_peakSWE <- merge(lim4_peakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim4_peakSWE <- lim4_peakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Four')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date Peak SWE
                  lim4_datePeakSWE <- as.data.frame(z_date_peak_swe$Water_Year[31:43])
                  colnames(lim4_datePeakSWE) <- 'Water_Year'
                  for (i in 2:length(z_date_peak_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_peak_swe, startyear = 2010, endyear = 2022, station = i)
                    lim4_datePeakSWE <- merge(lim4_datePeakSWE, result, by = 'Water_Year', all = TRUE)
                  }
                  lim4_datePeakSWE <- lim4_datePeakSWE %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Four')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Date 5 SWE
                  lim4_date5swe <- as.data.frame(z_date_5_swe$Water_Year[31:43])
                  colnames(lim4_date5swe) <- 'Water_Year'
                  for (i in 2:length(z_date_5_swe)) {
                    result <- LIM(E = pdo_matrix, F = z_date_5_swe, startyear = 2010, endyear = 2022, station = i)
                    lim4_date5swe <- merge(lim4_date5swe, result, by = 'Water_Year', all = TRUE)
                  }
                  lim4_date5swe <- lim4_date5swe %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Four')) %>%
                    relocate(Period, .after = Decade)
                  
                  # Melt Time
                  lim4_meltTime <- as.data.frame(z_melt_time$Water_Year[31:43])
                  colnames(lim4_meltTime) <- 'Water_Year'
                  for (i in 2:length(z_melt_time)) {
                    result <- LIM(E = pdo_matrix, F = z_melt_time, startyear = 2010, endyear = 2022, station = i)
                    lim4_meltTime <- merge(lim4_meltTime, result, by = 'Water_Year', all = TRUE)
                  }
                  lim4_meltTime <- lim4_meltTime %>%
                    mutate(Decade = Water_Year - Water_Year %% 10) %>%
                    relocate(Decade, .after = Water_Year) %>%
                    mutate(Period = as.factor('Four')) %>%
                    relocate(Period, .after = Decade)
                  
              # Join the Results together into one for each time frame and pivot them
              lim_peakSWE <- rbind(lim1_peakSWE, lim2_peakSWE, lim3_peakSWE, lim4_peakSWE)
              lim_datePeakSWE <- rbind(lim1_datePeakSWE, lim2_datePeakSWE, lim3_datePeakSWE, lim4_datePeakSWE)
              lim_date5SWE <- rbind(lim1_date5swe, lim2_date5swe, lim3_date5swe, lim4_date5swe)
              lim_meltTime <- rbind(lim1_meltTime, lim2_meltTime, lim3_meltTime, lim4_meltTime)
              rm(lim1_peakSWE, lim1_datePeakSWE, lim1_date5swe, lim1_meltTime, 
                 lim2_peakSWE, lim2_datePeakSWE, lim2_date5swe, lim2_meltTime, 
                 lim3_peakSWE, lim3_datePeakSWE, lim3_date5swe, lim3_meltTime, 
                 lim4_peakSWE, lim4_datePeakSWE, lim4_date5swe, lim4_meltTime)
              
              lim_peakSWE <- lim_peakSWE %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade, Period), names_to = 'Station', values_to = 'LIMPeakSWE')
              
              lim_datePeakSWE <- lim_datePeakSWE %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade, Period), names_to = 'Station', values_to = 'LIMDatePeakSWE')
              
              lim_date5SWE <- lim_date5SWE %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade, Period), names_to = 'Station', values_to = 'LIMDate5SWE')
              
              lim_meltTime <- lim_meltTime %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade, Period), names_to = 'Station', values_to = 'LIMMeltTime')
                  
   
              # Create a master dataframe for plotting of all LIM results
              # Want a long column for each time period of zPeak SWE, zDate Peak SWE, zDate5SWE, zMelttime
              
              # Pivot the Zscore Dataframes
              z_peak_swe <- z_peak_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              z_peak_swe$Decade[41:43] <- 2010
              z_peak_swe <- z_peak_swe %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'PeakSWE')
              
              z_date_peak_swe <- z_date_peak_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              z_date_peak_swe$Decade[41:43] <- 2010
              z_date_peak_swe <- z_date_peak_swe %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'DatePeakSWE')
              
              z_date_5_swe <- z_date_5_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              z_date_5_swe$Decade[41:43] <- 2010
              z_date_5_swe <- z_date_5_swe %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'Date5SWE')
              
              z_melt_time <- z_melt_time %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              z_melt_time$Decade[41:43] <- 2010
              z_melt_time <- z_melt_time %>%
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'MeltTime')
                  
              
              # Pivot the SNOTEL Metadata
              data <- data.frame(matrix(nrow = 13459, ncol = 6)) 
              colnames(data) <- c('Station', 'Longitude', 'Latitude', 'Elevation', 'HUC8', 'HUC2')
              for (i in 1:length(z_peak_swe$Station)) {
                station <- z_peak_swe$Station[i]
                index <- which(swe_meta$StationId == station) # index for the station in question 
                huc8 <- swe_meta$HUC8[index] # the huc that corresponds to that station
                huc2 <- swe_meta$HUC2[index] # the HUC2 that corresponds to that station 
                elev <- swe_meta$Elevation[index] # the elevation that corresponds to that station
                lon <- swe_meta$Longitude[index]
                lat <- swe_meta$Latitude[index]
                data[i, 1] <- station
                data[i, 2] <- lon
                data[i, 3] <- lat
                data[i, 4] <- elev
                data[i, 5] <- huc8 
                data[i, 6] <- huc2
                rm(station, index, huc8, huc2, elev, lon, lat)
              }
              
              # Pivot the Mean PDO data
              pdo_pivot <- data.frame(matrix(nrow = 13459, ncol = 1))
              colnames(pdo_pivot) <- c('MeanPDO')
              for (i in 1:length(z_peak_swe$Water_Year)) {
                year <- z_peak_swe$Water_Year[i]
                val <- mean_pdo$PDO[mean_pdo$Water_Year == year]
                pdo_pivot[i, 1] <- val
                rm(year, val)
              }

              # Create a Master Dataframe
              master <- cbind(z_peak_swe, data[, 2:6], pdo_pivot, z_date_peak_swe$DatePeakSWE, z_date_5_swe$Date5SWE, z_melt_time$MeltTime,
                              lim_peakSWE[, c(3,5)], lim_datePeakSWE$LIMDatePeakSWE, lim_date5SWE$LIMDate5SWE, lim_meltTime$LIMMeltTime)
              master <- master %>%
                relocate(PeakSWE, .after = MeanPDO) %>%
                relocate(Period, .after = Decade)
              colnames(master) <- c('Water_Year', 'Decade', 'Period', 'Station', 'Longitude', 'Latitude', 'Elevation', 'HUC8', 'HUC2', 'MeanPDO', 'PeakSWE', 
                                    'DatePeakSWE', 'Date5SWE', 'MeltTime', 'LIMPeakSWE', 'LIMDatePeakSWE', 'LIMDate5SWE', 'LIMMeltTime')
              master <- master %>%
                mutate(Water_Year = as.factor(Water_Year),
                       Station = as.factor(Station),
                       HUC8 = as.factor(HUC8),
                       HUC2 = as.factor(HUC2))
              
              rm(data, lim_peakSWE, lim_datePeakSWE, lim_date5SWE, lim_meltTime)
              
              

############################### Model vs Observed Plots ######################## 
               master1 <- master %>%
                filter(Period == 'One')
              master2 <- master %>%
                filter(Period == 'Two')
              master3 <- master %>%
                filter(Period == 'Three')
              master4 <- master %>%
                filter(Period == 'Four')
              
              p1 <- ggplot(master1, aes(x = LIMPeakSWE, y = PeakSWE, group = as.factor(Period), color = Water_Year)) +
                labs(x = 'Modeled Peak SWE Z-Score', y = 'Observed Peak SWE Z-Score', title = '1980s: Observed Peak SWE vs Modeled Peak SWE', color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5,3)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3,4)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 16),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p2 <- ggplot(master2, aes(x = LIMPeakSWE, y = PeakSWE, group = as.factor(Period), color = Water_Year)) +
                labs(x = 'Modeled Peak SWE Z-Score', y = 'Observed Peak SWE Z-Score', title = '1990s: Observed Peak SWE vs Modeled Peak SWE', color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-3,5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3,5)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 16),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p3 <- ggplot(master3, aes(x = LIMPeakSWE, y = PeakSWE, group = as.factor(Period), color = Water_Year)) +
                labs(x = 'Modeled Peak SWE Z-Score', y = 'Observed Peak SWE Z-Score', title = '2000s: Observed Peak SWE vs Modeled Peak SWE', color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5,5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3,3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 16),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) +
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p4 <- ggplot(master4, aes(x = LIMPeakSWE, y = PeakSWE, group = as.factor(Period), color = Water_Year)) +
                labs(x = 'Modeled Peak SWE Z-Score', y = 'Peak SWE Z-Score', title = '2010s: Observed Peak SWE vs Modeled Peak SWE', color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                                              "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 'grey')) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3,6)) +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 16),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = 'AUTO', rel_heights = c(1,1,1,1), label_size = 20)
              ggsave(filename = 'PeakModObs.png', width = 15, height = 15, p, path = path5)
              

              # Date of Peak SWE
              p1 <- ggplot(master1, aes(x = LIMDatePeakSWE, y = DatePeakSWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of Peak SWE Z-Score', y = 'Observed Date of Peak SWE Z-Score', title = '1980s: Observed Date of Peak SWE vs Modeled Date of Peak SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3,3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p2 <- ggplot(master2, aes(x = LIMDatePeakSWE, y = DatePeakSWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of Peak SWE Z-Score', y = 'Observed Date of Peak SWE Z-Score', title = '1990s: Observed Date of Peak SWE vs Modeled Date of Peak SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 0.5), limits = c(-3, 3)) +
                scale_y_continuous(breaks = seq(-5, 5, 1), limits = c(-3, 3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p3 <- ggplot(master3, aes(x = LIMDatePeakSWE, y = DatePeakSWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of Peak SWE Z-Score', y = 'Observed Date of Peak SWE Z-Score', title = '2000s: Observed Date of Peak SWE vs Modeled Date of Peak SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) +  
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p4 <- ggplot(master4, aes(x = LIMDatePeakSWE, y = DatePeakSWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of Peak SWE Z-Score', y = 'Observed Date of Peak SWE Z-Score', title = '2010s: Observed Date of Peak SWE vs Modeled Date of Peak SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                                              "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 'grey')) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-4, 3)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-4, 3)) +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = 'AUTO', rel_heights = c(1,1,1,1), label_size = 20)
              ggsave(filename = 'DatePeakModObs.png', width = 15, height = 15, p, path = path5)
              
              
              # Date of 5% SWE
              p1 <- ggplot(master1, aes(x = LIMDate5SWE, y = Date5SWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of 5% SWE Z-Score', y = 'Observed Date of 5% SWE Z-Score', title = '1980s: Observed Date of 5% SWE vs Modeled Date of 5% SWE', 
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p2 <- ggplot(master2,  aes(x = LIMDate5SWE, y = Date5SWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of 5% SWE Z-Score', y = 'Observed Date of 5% SWE Z-Score', title = '1990s: Observed Date of 5% SWE vs Modeled Date of 5% SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-4, 3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p3 <- ggplot(master3, aes(x = LIMDate5SWE, y = Date5SWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of 5% SWE Z-Score', y = 'Observed Date of 5% SWE Z-Score', title = '2000s: Observed Date of 5% SWE vs Modeled Date of 5% SWE',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 3)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p4 <- ggplot(master4,  aes(x = LIMDate5SWE, y = Date5SWE, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Date of 5% SWE Z-Score', y = 'Observed Date of 5% SWE Z-Score', title = '2010s: Observed Date of 5% SWE vs Modeled Date of 5% SWE', 
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                                              "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 'grey')) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 3)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 3)) +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = 'AUTO', rel_heights = c(1,1,1,1), label_size = 20)
              ggsave(filename = 'Date5ModObs.png', width = 15, height = 15, p, path = path5)
              
              
              # Melt Duration
              p1 <- ggplot(master1, aes(x = LIMMeltTime, y = MeltTime, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Melt Duration Z-Score', y = 'Observed Melt Duration Z-Score', title = '1980s: Observed Melt Duration vs Observed Melt Duration', 
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 4)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6, show.legend = FALSE) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p2 <- ggplot(master2,   aes(x = LIMMeltTime, y = MeltTime, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Melt Duration Z-Score', y = 'Observed Melt Duration Z-Score', title = '1990s: Observed Melt Duration vs Modeled Melt Duration',
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 5)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6, show.legend = FALSE) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p3 <- ggplot(master3,  aes(x = LIMMeltTime, y = MeltTime, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Melt Duration Z-Score', y = 'Observed Melt Duration Z-Score', title = '2000s: Observed Melt Duration vs Modeled Melt Duration', 
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-5, 5)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-4, 4)) +
                scale_color_brewer(palette = 'Paired') +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6, show.legend = FALSE) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p4 <- ggplot(master4, aes(x = LIMMeltTime, y = MeltTime, group = Period, color = Water_Year)) +
                labs(x = 'Modeled Melt Duration Z-Score', y = 'Observed Melt Duration Z-Score', title = '2010s: Observed Melt Duration vs Modeled Melt Duration', 
                     color = 'Water Year') +
                geom_point(na.rm = TRUE, alpha = 0.7) +
                scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                                              "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", 'grey')) +
                scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 3)) +
                scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-3, 3)) +
                geom_hline(aes(yintercept = 0), color = 'black') + 
                geom_vline(aes(xintercept = 0), color = 'black') +
                stat_poly_line(method = 'lm', se = FALSE, na.rm = TRUE, color = 'red', linetype = 'dashed', linewidth = 2) +
                stat_poly_eq(use_label(c("eq", "R2", 'P')), label.x = 'left', label.y = 'top', color = 'red', size = 6, show.legend = FALSE, na.rm = TRUE) +
                # stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = '<br>')), 
                #              label.x = 'left', label.y = 'top', color = 'red', size = 6, parse = TRUE) +
                geom_abline(intercept = 0, slope = 1) +
                theme_clean() + 
                theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
                      plot.title = element_text(size = 15),
                      strip.text = element_text(size = 12),
                      plot.caption = ggtext::element_markdown(size = 15, hjust = 0),
                      legend.title = element_text(size = 11),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      axis.title = element_text(size = 16),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14)) + 
                guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
              
              p <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = 'AUTO', rel_heights = c(1,1,1,1), label_size = 20)
              ggsave(filename = 'MeltTimeModObs.png', width = 15, height = 15, p, path = path5)
              
              
################################### Variance Reduction Mapping #################
          
              # Compute the Variance Reduction
                  peak_swe_stats <- master %>%
                    group_by(Station, Period) %>%
                    summarise(varRed = (1 - (var(PeakSWE - LIMPeakSWE))/(var(PeakSWE) - 0))* 100) %>%
                    mutate_if(is.numeric, round, digits = 2)
                  date_peak_swe_stats <- master %>%
                    group_by(Station, Period) %>%
                    summarise(varRed = (1 - (var(DatePeakSWE - LIMDatePeakSWE))/(var(DatePeakSWE) - 0))* 100) %>%
                    mutate_if(is.numeric, round, digits = 2)
                  date_5_swe_stats <- master %>%
                    group_by(Station, Period) %>%
                    summarise(varRed = (1 - (var(Date5SWE - LIMDate5SWE))/(var(Date5SWE) - 0))* 100) %>%
                    mutate_if(is.numeric, round, digits = 2)
                  melt_time_stats <- master %>%
                    group_by(Station, Period) %>%
                    summarise(varRed = (1 - (var(MeltTime - LIMMeltTime))/(var(MeltTime) - 0))* 100) %>%
                    mutate_if(is.numeric, round, digits = 2)
              
              
              # Pivot the Metadata
                  data <- data.frame(matrix(nrow = 1252, ncol = 5)) # 1252
                  colnames(data) <- c('Longitude', 'Latitude', 'Elevation', 'HUC8', 'HUC2')
                  for (i in 1:length(peak_swe_stats$Station)) {
                    station <- peak_swe_stats$Station[i]
                    index <- which(swe_meta$StationId == station) # index for the station in question 
                    huc8 <- swe_meta$HUC8[index] # the huc that corresponds to that station
                    huc2 <- swe_meta$HUC2[index] # the HUC2 that corresponds to that station 
                    elev <- swe_meta$Elevation[index] # the elevation that corresponds to that station
                    lon <- swe_meta$Longitude[index]
                    lat <- swe_meta$Latitude[index]
                    data[i, 1] <- lon
                    data[i, 2] <- lat
                    data[i, 3] <- elev
                    data[i, 4] <- huc8 
                    data[i, 5] <- huc2
                    rm(station, index, huc8, huc2, elev, lon, lat)
                  }
                  
              # Bind the statistics to the metadata
                  # Peak SWE
                  peak_swe_stats <- cbind(data, peak_swe_stats)
                  peak_swe_stats <- peak_swe_stats %>%
                    relocate(Station) %>%
                    mutate(ElevationBinned = findInterval(Elevation, c(0, 1000, 1500, 2000, 2500, 3000, 3500), all.inside = TRUE)) %>%
                    relocate(ElevationBinned, .after = Elevation)  %>% # Bin the Elevations
                    mutate(ElevationBinned = factor(ElevationBinned, labels = c('0-1000', '1000-1500', '1500-2000', '2000-2500', '2500-3000', '3000-3500'))) %>%
                    mutate(LatBinned = findInterval(Latitude, c(32, 35, 38, 41, 44, 47, 50), all.inside = TRUE)) %>%
                    relocate(LatBinned, .after = Latitude)  %>% # Bin the Latitude
                    mutate(LatBinned = factor(LatBinned, labels = c('32-35', '35-38', '38-41', '41-44', '44-47', '47-50'))) %>%
                    mutate(LonBinned = findInterval(Longitude, c(-124, -120, -115, -110, -105), all.inside = TRUE)) %>%
                    relocate(LonBinned, .after = Longitude)  %>% # Bin the Longitude
                    mutate(LonBinned = factor(LonBinned, labels = c('-124 to -120', '-120 to -115', '-115 to -110', '-110 to -105'))) 
                  #rm(data)

                  # Date of Peak SWE
                  date_peak_swe_stats <- cbind(data, date_peak_swe_stats)
                  date_peak_swe_stats <- date_peak_swe_stats %>%
                    relocate(Station) %>%
                    mutate(ElevationBinned = findInterval(Elevation, c(0,1000, 1500, 2000, 2500, 3000, 3500), all.inside = TRUE)) %>%
                    relocate(ElevationBinned, .after = Elevation)  %>% # Bin the Elevations
                    mutate(ElevationBinned = factor(ElevationBinned, labels = c('0-1000', '1000-1500', '1500-2000', '2000-2500', '2500-3000', '3000-3500'))) %>%
                    mutate(LatBinned = findInterval(Latitude, c(32, 35, 38, 41, 44, 47, 50), all.inside = TRUE)) %>%
                    relocate(LatBinned, .after = Latitude)  %>% # Bin the Latitude
                    mutate(LatBinned = factor(LatBinned, labels = c('32-35', '35-38', '38-41', '41-44', '44-47', '47-50'))) %>%
                    mutate(LonBinned = findInterval(Longitude, c(-124, -120, -115, -110, -105), all.inside = TRUE)) %>%
                    relocate(LonBinned, .after = Longitude)  %>% # Bin the Longitude
                    mutate(LonBinned = factor(LonBinned, labels = c('-124 to -120', '-120 to -115', '-115 to -110', '-110 to -105'))) 
                  #rm(data)
              
                  # Date of 5% Peak SWE
                  date_5_swe_stats <- cbind(data, date_5_swe_stats)
                  date_5_swe_stats <- date_5_swe_stats %>%
                    relocate(Station) %>%
                    mutate(ElevationBinned = findInterval(Elevation, c(0,1000, 1500, 2000, 2500, 3000, 3500), all.inside = TRUE)) %>%
                    relocate(ElevationBinned, .after = Elevation)  %>% # Bin the Elevations
                    mutate(ElevationBinned = factor(ElevationBinned, labels = c('0-1000', '1000-1500', '1500-2000', '2000-2500', '2500-3000', '3000-3500'))) %>%
                    mutate(LatBinned = findInterval(Latitude, c(32, 35, 38, 41, 44, 47, 50), all.inside = TRUE)) %>%
                    relocate(LatBinned, .after = Latitude)  %>% # Bin the Latitude
                    mutate(LatBinned = factor(LatBinned, labels = c('32-35', '35-38', '38-41', '41-44', '44-47', '47-50'))) %>%
                    mutate(LonBinned = findInterval(Longitude, c(-124, -120, -115, -110, -105), all.inside = TRUE)) %>%
                    relocate(LonBinned, .after = Longitude)  %>% # Bin the Longitude
                    mutate(LonBinned = factor(LonBinned, labels = c('-124 to -120', '-120 to -115', '-115 to -110', '-110 to -105'))) 
                  #rm(data)
              
                  # Melt Duration
                  melt_time_stats <- cbind(data, melt_time_stats)
                  melt_time_stats <- melt_time_stats %>%
                    relocate(Station) %>%
                    mutate(ElevationBinned = findInterval(Elevation, c(0,1000, 1500, 2000, 2500, 3000, 3500), all.inside = TRUE)) %>%
                    relocate(ElevationBinned, .after = Elevation)  %>% # Bin the Elevations
                    mutate(ElevationBinned = factor(ElevationBinned, labels = c('0-1000', '1000-1500', '1500-2000', '2000-2500', '2500-3000', '3000-3500'))) %>%
                    mutate(LatBinned = findInterval(Latitude, c(32, 35, 38, 41, 44, 47, 50), all.inside = TRUE)) %>%
                    relocate(LatBinned, .after = Latitude)  %>% # Bin the Latitude
                    mutate(LatBinned = factor(LatBinned, labels = c('32-35', '35-38', '38-41', '41-44', '44-47', '47-50'))) %>%
                    mutate(LonBinned = findInterval(Longitude, c(-124, -120, -115, -110, -105), all.inside = TRUE)) %>%
                    relocate(LonBinned, .after = Longitude)  %>% # Bin the Longitude
                    mutate(LonBinned = factor(LonBinned, labels = c('-124 to -120', '-120 to -115', '-115 to -110', '-110 to -105'))) 
                  rm(data)
              
             
              
              # Plot Variance Reduction
                  # Peak SWE
                  p <- ggplot(us, aes(long, lat, group = group)) +
                    geom_polygon(fill = 'white', color = 'grey50') +
                    geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                    coord_sf() +
                    geom_point(data = peak_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed, size = Elevation), shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                    geom_point(data = peak_swe_stats %>% filter(Period == 'Four' & varRed < 0), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 3) + 
                    labs(x = 'Longitude', y = 'Latitude', title = 'The Role of PDO in Peak SWE Magnitude: 2010-2022', 
                         caption = '**Figure #** Map of the variance reduction (%) computed for the role of the PDO in<br>Period Four (2010-2022). Greater percentages indicate more variance in the<br>peak SWE magnitude explained by the model. Grey dots represent locations<br>where the PDO had no significant input. Shaded regions are large river basins<br>(HUC2s) that are numbered here.') +
                    theme(plot.title = element_text(size = 20),
                          plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                          plot.subtitle = ggtext::element_markdown(size = 14),
                          legend.title = element_text(size = 14), 
                          legend.text = element_text(size = 11), 
                          legend.position = 'right', 
                          legend.key.size = unit(1, 'cm')) +
                    geom_label(data = labels, aes(x = longitude, y = latitude, label = label), inherit.aes = FALSE) +
                    scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
                  ggsave(filename = 'varRedPeakSWE4.png', width = 9, height = 9, p, path = path5)
                  
                  # Date of Peak SWE
                  p <- ggplot(us, aes(long, lat, group = group)) +
                    geom_polygon(fill = 'white', color = 'grey50') +
                    geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                    coord_sf() +
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed, size = Elevation), shape = 21, inherit.aes = FALSE, alpha = 0.8) + # size = 3
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & varRed < 0), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 3) + 
                    labs(x = 'Longitude', y = 'Latitude', title = 'The Role of PDO in Date of Peak SWE: 2010-2022', 
                         caption = '**Figure #** Map of the variance reduction (%) computed for the role of the PDO in<br>Period Four (2010-2022). Greater percentages indicate more variance in the<br>date of peak SWE explained by the model. Grey dots represent locations<br>where the PDO had no significant input. Shaded regions are large river basins<br>(HUC2s) that are numbered here.') +
                    theme(plot.title = element_text(size = 20),
                          plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                          plot.subtitle = ggtext::element_markdown(size = 14),
                          legend.title = element_text(size = 14), 
                          legend.text = element_text(size = 11), 
                          legend.position = 'right', 
                          legend.key.size = unit(1, 'cm')) +
                    geom_label(data = labels, aes(x = longitude, y = latitude, label = label), inherit.aes = FALSE) +
                    scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
                  ggsave(filename = 'varRedDatePeakSWE2.png', width = 9, height = 9, p, path = path5) # varRedDatePeakSWE4
                  
                  
                  # Date of 5% Peak SWE
                  p <- ggplot(us, aes(long, lat, group = group)) +
                    geom_polygon(fill = 'white', color = 'grey50') +
                    geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                    coord_sf() +
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed, size = Elevation), shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & varRed < 0), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 3) + 
                    labs(x = 'Longitude', y = 'Latitude', title = 'The Role of PDO in SWE Melt Date: 2010-2022', 
                         caption = '**Figure #** Map of the variance reduction (%) computed for the role of the PDO in<br>Period Four (2010-2022). Greater percentages indicate more variance in the<br>date of 5% peak SWE explained by the model. Grey dots represent locations<br>where the PDO had no significant input. Shaded regions are large river basins<br>(HUC2s) that are numbered here.') +
                    theme(plot.title = element_text(size = 20),
                          plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                          plot.subtitle = ggtext::element_markdown(size = 14),
                          legend.title = element_text(size = 14), 
                          legend.text = element_text(size = 11), 
                          legend.position = 'right', 
                          legend.key.size = unit(1, 'cm')) +
                    geom_label(data = labels, aes(x = longitude, y = latitude, label = label), inherit.aes = FALSE) +
                    scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
                  ggsave(filename = 'varRedDate5SWE4.png', width = 9, height = 9, p, path = path5)
                  
                  # Melt Duration
                  p <- ggplot(us, aes(long, lat, group = group)) +
                    geom_polygon(fill = 'white', color = 'grey50') +
                    geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                    geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                    coord_sf() +
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Two' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed, size = Elevation), shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                    geom_point(data = date_peak_swe_stats %>% filter(Period == 'Two' & varRed < 0), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 3) + 
                    labs(x = 'Longitude', y = 'Latitude', title = 'The Role of PDO in SWE Melt Duration: 1990-1999', 
                         caption = '**Figure #** Map of the variance reduction (%) computed for the role of the PDO in<br>Period Two (1990-1999). Greater percentages indicate more variance in the<br> SWE melt duration explained by the model. Grey dots represent locations<br>where the PDO had no significant input. Shaded regions are large river basins<br>(HUC2s) that are numbered here.') +
                    theme(plot.title = element_text(size = 20),
                          plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                          plot.subtitle = ggtext::element_markdown(size = 14),
                          legend.title = element_text(size = 14), 
                          legend.text = element_text(size = 11), 
                          legend.position = 'right', 
                          legend.key.size = unit(1, 'cm')) +
                    geom_label(data = labels, aes(x = longitude, y = latitude, label = label), inherit.aes = FALSE) +
                    scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
                  ggsave(filename = 'varRedMeltTime2.png', width = 9, height = 9, p, path = path5)
              
                
                
                
                  
                  
              
              
                
            
##################################### Create Panel Plots for Poster ################################
              # Period One
              p1 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1980-1989') +
                geom_point(data = peak_swe_stats %>% filter(Period == 'One' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = peak_swe_stats %>% filter(Period == 'One' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
             
              # Period 2
              p2 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1990-1999') +
                geom_point(data = peak_swe_stats %>% filter(Period == 'Two' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = peak_swe_stats %>% filter(Period == 'Two' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
    

              # Period 3
              p3 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2000-2009') +
                geom_point(data = peak_swe_stats %>% filter(Period == 'Three' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = peak_swe_stats %>% filter(Period == 'Three' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 4
              p4 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2010-2022') +
                geom_point(data = peak_swe_stats %>% filter(Period == 'Four' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = peak_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none' ) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
    
              legend <- get_legend(p4 + theme(legend.position = 'right',
                                              legend.key.height = unit(1.665, 'cm')))
              p <- plot_grid(p1, p2, p3, p4, legend, ncol = 5, nrow = 1, rel_heights = c(1,1,1,1, 3), rel_widths = c(1,1,1,1,0.38))
              ggsave(filename = 'VR_PeakSWE.png', width = 16, height = 4, p, path = path5)
              
              
              
              # Period One
              p1 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1980-1989') +
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'One' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'One' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              # Period 2
              p2 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1990-1999') +
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Two' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Two' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 3
              p3 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2000-2009') +
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Three' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Three' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 4
              p4 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2010-2022') +
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_peak_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none' ) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              legend <- get_legend(p4 + theme(legend.position = 'right',
                                              legend.key.height = unit(1.665, 'cm')))
              p <- plot_grid(p1, p2, p3, p4, legend, ncol = 5, nrow = 1, rel_heights = c(1,1,1,1, 3), rel_widths = c(1,1,1,1,0.38))
              ggsave(filename = 'VR_DatePeakSWE.png', width = 16, height = 4, p, path = path5)
              
              
              
              # Period One
              p1 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1980-1989') +
                geom_point(data = date_5_swe_stats %>% filter(Period == 'One' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_5_swe_stats %>% filter(Period == 'One' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              # Period 2
              p2 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1990-1999') +
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Two' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Two' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 3
              p3 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2000-2009') +
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Three' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Three' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 4
              p4 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2010-2022') +
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Four' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = date_5_swe_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none' ) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              legend <- get_legend(p4 + theme(legend.position = 'right',
                                              legend.key.height = unit(1.665, 'cm')))
              p <- plot_grid(p1, p2, p3, p4, legend, ncol = 5, nrow = 1, rel_heights = c(1,1,1,1, 3), rel_widths = c(1,1,1,1,0.38))
              ggsave(filename = 'VR_Date5SWE.png', width = 16, height = 4, p, path = path5)
              
              
              
              
              # Period One
              p1 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1980-1989') +
                geom_point(data = melt_time_stats %>% filter(Period == 'One' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = melt_time_stats %>% filter(Period == 'One' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              # Period 2
              p2 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '1990-1999') +
                geom_point(data = melt_time_stats %>% filter(Period == 'Two' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = melt_time_stats %>% filter(Period == 'Two' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 3
              p3 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2000-2009') +
                geom_point(data = melt_time_stats %>% filter(Period == 'Three' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = melt_time_stats %>% filter(Period == 'Three' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size = 16),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none', 
                      legend.key.size = unit(1, 'cm')) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              
              # Period 4
              p4 <- ggplot(us, aes(long, lat, group = group)) +
                geom_polygon(fill = 'white', color = 'grey50') +
                geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
                geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
                coord_sf() +
                labs(x = 'Longitude', y = 'Latitude', title = '2010-2022') +
                geom_point(data = melt_time_stats %>% filter(Period == 'Four' & (varRed < 0 | is.na(varRed))), mapping = aes(x = Longitude, y = Latitude), shape = 21, inherit.aes = FALSE, fill = 'grey', size = 1) + 
                geom_point(data = melt_time_stats %>% filter(Period == 'Four' & varRed > 0), mapping = aes(x = Longitude, y = Latitude, fill = varRed), size = 3, shape = 21, inherit.aes = FALSE, alpha = 0.8) + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      plot.subtitle = ggtext::element_markdown(size = 14),
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'none' ) +
                scale_fill_viridis('Variance\nReduction (%)', option = 'plasma')
              
              legend <- get_legend(p4 + theme(legend.position = 'right',
                                              legend.key.height = unit(1.665, 'cm')))
              p <- plot_grid(p1, p2, p3, p4, legend, ncol = 5, nrow = 1, rel_heights = c(1,1,1,1, 3), rel_widths = c(1,1,1,1,0.38))
              ggsave(filename = 'VR_MeltTime.png', width = 16, height = 4, p, path = path5)
             
              
##################################### Create Master dataframe for Plotting #########################
              # Pivot peak SWE magnitude
              peak_swe <- peak_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              peak_swe$Decade[41:43] <- 2010
              peak_swe <- peak_swe %>%  
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'PeakSWE')
              
              # Pivot Date of Peak SWE
              date_peak_swe <- date_peak_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              date_peak_swe$Decade[41:43] <- 2010
              date_peak_swe <- date_peak_swe %>%  
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'DatePeakSWE')
              
              # Pivot Date of 5% SWE
              date_5_swe <- date_5_swe %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              date_5_swe$Decade[41:43] <- 2010
              date_5_swe <- date_5_swe %>%  
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'Date5SWE')
              
              # Pivot Melt Duration
              melt_time <- melt_time %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year)
              melt_time$Decade[41:43] <- 2010
              melt_time <- melt_time %>%  
                mutate(Decade = as.factor(Decade)) %>%
                pivot_longer(-c(Water_Year, Decade), names_to = 'Station', values_to = 'MeltTime')
              
              # Pivot All Metadata
              data <- data.frame(matrix(nrow = 13459, ncol = 5))
              colnames(data) <- c('Longitude', 'Latitude', 'Elevation', 'HUC8', 'HUC2')
              for (i in 1:length(peak_swe$Station)) {
                station <- peak_swe$Station[i]
                index <- which(swe_meta$StationId == station) # index for the station in question 
                huc8 <- swe_meta$HUC8[index] # the huc that corresponds to that station
                huc2 <- swe_meta$HUC2[index] # the HUC2 that corresponds to that station 
                elev <- swe_meta$Elevation[index] # the elevation that corresponds to that station
                lon <- swe_meta$Longitude[index]
                lat <- swe_meta$Latitude[index]
                data[i, 1] <- lon
                data[i, 2] <- lat
                data[i, 3] <- elev
                data[i, 4] <- huc8 
                data[i, 5] <- huc2
                rm(station, index, huc8, huc2, elev, lon, lat)
              }
              
              # Pivot the PDO data
              df <- data.frame(matrix(nrow = 13459, ncol = 1))
              colnames(df) <- c('MeanPDO')
              for (i in 1:length(peak_swe$Water_Year)) {
                year <- peak_swe$Water_Year[i]
                index <- which(mean_pdo$Water_Year == year)
                val <- mean_pdo$PDO[index]
                df[i, 1] <- as.numeric(val)
              }
              
              # Join them all together
              master <- cbind(peak_swe, data)
              master <- master %>%
                relocate(c(5:9), .after = Station)
              master <- cbind(master, date_peak_swe$DatePeakSWE, date_5_swe$Date5SWE, melt_time$MeltTime, df$MeanPDO)
              colnames(master) <- c('Water_Year', 'Decade', 'Station', 'Longitude', 'Latitude', 'Elevation', 'HUC8', 'HUC2', 
                                    'PeakSWE', 'DatePeakSWE', 'Date5SWE', 'MeltTime', 'MeanPDO')
              master <- master %>%
                mutate(Water_Year = as.factor(Water_Year),
                       Station = as.factor(Station),
                       HUC8 = as.factor(HUC8),
                       HUC2 = as.factor(HUC2),
                       DatePeakSWE = as.numeric(DatePeakSWE), 
                       Date5SWE = as.numeric(Date5SWE)) 

########################################## Burgess Junction, Wyoming: SNOTEL 377 #################################
              
          # Burgess Junction SWE
              burg_swe <- swe %>%
                dplyr::select(c(1:5, '377')) %>%
                mutate(Decade = Water_Year - Water_Year %% 10) %>%
                relocate(Decade, .after = Water_Year) %>%
                filter(Water_Year == 2021) %>%
                mutate(Date = as.Date(Date))
              
              plot <- ggplot(burg_swe, aes(x = Date, y = `377`, group = Water_Year)) + 
                geom_area(na.rm = TRUE, fill = 'lightblue', color = 'black') + 
                scale_x_date(date_breaks = '1 month') + 
                labs(x = 'Date', y = 'Snow Water Equivalent (mm)', title = 'SWE Hydrograph at SNOTEL 377, Burgess Junction: Water Year 2021',
                     caption =  "**Figure 5** Snowpack hydrograph for water year 2021 at Burgess Junction, Wyoming. Light blue area is the raw snow water equivalent value<br>collected directly from the SNOTEL instrument. Red line is the date of peak SWE, April 27. Blue line is the date by which 5% of the<br>peak SWE value remains, June 6. Here, the melt duration is 36 days.") +
                theme_clean() + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 11),
                      legend.position = 'bottom',
                      legend.key.size = unit(1, 'cm'),
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 12, vjust = 0.9, hjust = 1, angle = 30),
                      axis.text.y = element_text(size = 14)) +
                geom_vline(aes(xintercept = as.Date('2021-06-02')), color = 'blue') + 
                geom_vline(aes(xintercept = as.Date('2021-04-27')), color = 'red')
              ggsave(filename = 'Burgess_swe_hydrograph.png', width = 13, height = 9, plot, path = path5)
              
              
          # Create df of all characteristics for this station
              burgess <- master %>%
                filter(Station == 377)
              burgess[burgess == 0] <- NA
        
          # Plots for Burgess Junction, Wyoming
              ggplot(burgess, aes(x = PeakSWE, y = DatePeakSWE, group = Station)) +
                geom_point(na.rm = TRUE) + 
                geom_smooth()
              
              
              # Duration bars and peak SWE line
              plot <- ggplot(burgess) + 
                geom_bar(mapping = aes(x = Water_Year, y = Duration), stat = 'identity', na.rm = TRUE, fill = 'pink') + 
                geom_line(mapping = aes(x = Water_Year, y = Mag), stat = 'identity', na.rm = TRUE, color = 'blue', linewidth = 2) +
                scale_y_continuous(name = 'Peak SWE (centimeters)', 
                                   sec.axis = sec_axis(~., name = 'Days from Peak SWE to 5% SWE')) + # Create a second axis
                geom_line(blah, mapping = aes(x = Water_Year, y = PDO)) +
                labs(x = 'Water Year', title = 'Number of Days from Peak SWE to 5% SWE, Peak SWE Magnitude: 1980-2022',
                     caption = 'Figure 7: SNOTEL 377, Burgess Junction, Wyoming. Pink bars are the number of days that elapsed each year from the date of the peak SWE <br>to the date of 5% SWE. Blue line is the magnitude, in centimeters, of peak SWE for each year. <br>Time series comprises each water year from 1980 through 2022.') + 
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'right', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle = 30, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14),
                      axis.line.y.right = element_line(color = "pink", size = 2), 
                      axis.ticks.y.right = element_line(color = "pink"),
                      axis.line.y.left = element_line(color = "blue", size = 2), 
                      axis.ticks.y.left = element_line(color = "blue"))
              ggsave(filename = 'Burgess_jctn_bars_and_line_aosc.png', width = 13, height = 8.5, plot, path = path5)
              
          
              # Duration from Peak SWE to 5% SWE vs the Magntiude of Peak SWE
              plot <- ggplot(burgess, aes(x = Mag, y = Duration)) + 
                stat_poly_line(se = FALSE) + 
                stat_poly_eq(use_label(c('eq', 'R2')), label.x = 'right') + 
                geom_point(aes(color = as.factor(Decade)), na.rm = TRUE, size = 3) + 
                labs(x = 'Peak Snow Water Equivalent (cm)', y = '# of Days from Peak SWE to 5% SWE',
                     title = 'Duration from Peak SWE to 5% SWE vs Magnitude of Peak SWE: SNOTEL 377', 
                     color = 'Decade',
                     caption = 'Figure 8: SNOTEL 377, Burgess Junction, Wyoming. Scatterplot of the duration between Peak SWE and 5% SWE versus the magnitude of Peak SWE. <br>Note that there is no noticeable significant correlation between the two variables.') + 
                theme_clean() +
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'bottom', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 12, angle = 30, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14))
              
              # Magnitude of Peak SWE over the time series
              plot <- ggplot(burgess, aes(x = Water_Year, y = Mag)) + 
               geom_line() + 
               geom_smooth() + 
               labs(x = 'Water Year', y = 'Magnitude of Peak SWE (cm)', title = 'Peak Snow Water Equivalent at Burgess Junction, Wyoming: 1980-2022',
                    caption = 'Figure 7: SNOTEL 377, Burgess Junction Wyoming. Line plot of the magnitude of peak SWE each year from 1980 through 2022. <br>A Loess-smoothing was applied at a 95% confidence interval. Note the gradual decreasing trend, and the two humps that <br>could indicate potential oscillations with the PDO.') + 
               theme_clean() +
               theme(plot.title = element_text(size = 20),
                   plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                   legend.title = element_text(size = 14), 
                   legend.text = element_text(size = 11), 
                   legend.position = 'bottom', 
                   legend.key.size = unit(1, 'cm'), 
                   axis.title = element_text(size = 18),
                   axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 1), 
                   axis.text.y = element_text(size = 14)) 
              ggsave(filename = 'Burgess_jctn_peakSWE_aosc.png', width = 13, height = 8.5, plot, path = path5)
              
              # Date of Peak SWE over the time Period
              plot <- ggplot(burgess, aes(x = Water_Year, y = Peak)) + 
                geom_line() + 
                geom_smooth() + 
                labs(x = 'Water Year', y = '# Of Days Since the Start of the Water Year', title = 'Date of Peak SWE at Burgess Junction, Wyoming: 1980-2022',
                     caption = 'Figure 8: SNOTEL 377, Burgess Junction Wyoming. Line plot of the date of peak SWE each year from 1980 through 2022. <br>A Loess-smoothing was applied at a 95% confidence interval. Note the gradual increasing trend, and the two humps that <br>could indicate potential oscillations with the PDO. The horizontal red lines represent April 18th and May 8th respectively.') + 
                theme_clean() +
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'bottom', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14)) +
                geom_hline(yintercept = 200, color = 'red') + 
                geom_hline(yintercept = 220, color = 'red')
              ggsave(filename = 'Burgess_jctn_date_peak_aosc.png', width = 13, height = 8.5, plot, path = path5)
              
              # Date of 5% SWE over the time Period
              plot <- ggplot(burgess, aes(x = Water_Year, y = Five)) + 
                geom_line() + 
                geom_smooth() + 
                labs(x = 'Water Year', y = '# Of Days Since the Start of the Water Year', title = 'SWE Melt Date at Burgess Junction, Wyoming: 1980-2022',
                     caption = 'Figure : SNOTEL 377, Burgess Junction Wyoming. Line plot of the SWE melt date each year from 1980 through 2022. <br>A Loess-smoothing was applied at a 95% confidence interval. Note the gradual increasing trend, and the two humps that <br>could indicate potential oscillations with the PDO. The horizontal red lines represent April 18th and May 8th respectively.') + 
                theme_clean() +
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'bottom', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14)) +
                geom_hline(yintercept = 240, color = 'red') 
              ggsave(filename = 'Burgess_jctn_date_peak_aosc.png', width = 13, height = 8.5, plot, path = path5)
              
              # Date of Peak SWE with Date of 5% SWE
              plot <- ggplot(burgess) + 
                geom_line(aes(x = Water_Year, y = Peak, color = 'blue'), na.rm = TRUE) +
                geom_line(aes(x = Water_Year, y = Five, color = 'red'), na.rm = TRUE) +
                labs(x = 'Water Year', y = '# Of Days Since the Start of the Water Year', title = 'Date of Peak SWE and SWE Melt Date at Burgess Junction, Wyoming: 1980-2022',
                     caption = 'Figure 9: SNOTEL 377, Burgess Junction Wyoming. Blue line represents the date of peak SWE, red line represents the melt date, <br>or the date by which the snowpack reaches 5% of its max value. Note the two horizontal black lines at 210 and 240 days <br>represent April 28 and May 29 respectively.') + 
                theme_clean() +
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'bottom', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14)) +
                geom_hline(yintercept = 240, color = 'black') + 
                geom_hline(yintercept = 210, color = 'black') + 
                scale_colour_manual(name = 'Color', 
                                    values =c('blue'='blue','red'='red'), labels = c('Date of Peak SWE','Date of 5% SWE'))
              ggsave(filename = 'Burgess_jctn_date_lines_aosc.png', width = 13, height = 8.5, plot, path = path5)
              
              # Date of Peak SWE vs Melt Date
              labs <- data.frame(x = c(198, 208, 218, 230, 230, 230),
                                 y = c(222, 222, 222, 231, 241, 251), 
                                 label = c('April 19', 'April 29', 'May 9', 'May 19', 'May 29', 'June 9'))
              plot <- ggplot(burgess, aes(x = Peak, y = Five)) +
                stat_poly_line(se = FALSE) +
                stat_poly_eq(use_label(c("eq", "R2")), label.x = 'right') +
                geom_point(aes(color = as.factor(Decade)), na.rm = TRUE, size = 3) + 
                labs(x = 'Date of Peak SWE (# of days since October 1)', y = 'Date of 5% SWE (# of Days since October 1)', 
                     title = 'Date of Peak SWE vs Melt Date, Burgess Junction, Wyoming: 1980-2022', 
                     color = 'Decade', 
                     caption = 'Figure 10: Plot of the date of peak SWE versus the melt date, AKA the date by which <5% of the peak SWE remains. The dots are colored <br>by the decade that they are a part of. Note the vertical black lines represent dates (see labels). All values are measured since the start of <br>the water year, October 1. A Loess linear regression model was fitted to all data points.') + 
                theme_clean() +
                theme(plot.title = element_text(size = 20),
                      plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                      legend.title = element_text(size = 14), 
                      legend.text = element_text(size = 11), 
                      legend.position = 'bottom', 
                      legend.key.size = unit(1, 'cm'), 
                      axis.title = element_text(size = 18),
                      axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 1), 
                      axis.text.y = element_text(size = 14)) + 
                geom_vline(xintercept = 200) + 
                geom_vline(xintercept = 210) + 
                geom_vline(xintercept = 220) + 
                geom_label(data = labs, aes(x = x, y = y, label = label)) 
              ggsave(filename = 'Burgess_jctn_date_vs_date_aosc.png', width = 13, height = 8.5, plot, path = path5)
                
                
            ggplot(pdo, aes(x = Date, y = PDO)) +
              geom_line()
              

                 

            
############################################### Plotting ####################################


          # Make map of filtered stations
          labels <- data.frame(longitude = c(-97, -103, -99, -106, -109, -112, -118, -119, -119.5),
                               latitude = c(48, 44, 36, 34, 39.5, 32, 40, 47, 36),
                               label = c('09', '10', '11', '13', '14', '15', '16', '17', '18'))
          # id_labels <- data.frame(longitude = c(-123, -107.9, -113.5, -121, -112),
          #                         latitude = c(41.24, 35.7, 41.9, 45.24, 48.17),
          #                         label = c('301', '1084', '1114','921', '917'))
          # Plot all Stations from all datasets
          station_map <- ggplot(us, aes(long, lat, group = group)) +
            geom_polygon(fill = 'white', color = 'grey50') +
            geom_sf(data = huc9, color  = '#E41A1C', fill = '#E41A1C', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc10, color  = '#377EB8', fill = '#377EB8', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc11, color  = '#4DAF4A', fill = '#4DAF4A', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc13, color  = '#984EA3', fill = '#984EA3', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc14, color  = '#FF7F00',fill = '#FF7F00',  alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc15, color  = '#FFFF33', fill = '#FFFF33', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc16, color  = '#A65628', fill = '#A65628', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc17, color  = '#F781BF', fill = '#F781BF', alpha = 0.6, inherit.aes = FALSE) +
            geom_sf(data = huc18, color = '#999999', fill = '#999999', alpha = 0.6, inherit.aes = FALSE) +
            coord_sf() +
            geom_point(swe_meta, mapping = aes(x = Longitude, y = Latitude, fill = Elevation), inherit.aes = FALSE, size = 3, shape = 25) +
            labs(x = 'Longitude', y = 'Latitude', title = 'SWE Monitoring Stations', 
                 caption = '**Figure 3**: Map of final set of SNOTEL stations. Total number of stations is 314, colored<br>by elevation. Shaded regions are large river basins (HUC2s) that are numbered here.') +
            theme(plot.title = element_text(size = 20),
                  plot.caption = ggtext::element_markdown(size = 14, hjust = 0), 
                  legend.title = element_text(size = 14), 
                  legend.text = element_text(size = 11), 
                  legend.position = 'right', 
                  legend.key.size = unit(1, 'cm')) +
            geom_label(data = labels, aes(x = longitude, y = latitude, label = label), inherit.aes = FALSE) +
            scale_fill_viridis('Elevation (m)', option = 'magma')
          ggsave(filename = 'filtered_stations_aosc.png', width = 11, height = 8.5, station_map, path = path5)

          
          
          
       



