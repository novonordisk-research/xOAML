##Longitudinal Data Binning Function##

##Inputs##

#longitudinal_data (txt files): 
#needed cols: eid (patient identifier, numeric), trait (biomarker name, character), value (biomarker value, numeric), date (date of measurement, date), (order does not matter)


#path to start_date (txt files):
#needed cols: eid (patient identifier, numeric), start_date (date of first diagnosis, date), (order does matter)

#closest_avg: character indicating whether to use closest value ("closest") to time point or median value from time window ("avg")

#post_bin_years: number of years to bin data post-start date

#pre_bin_years: number of years to bin data pre-start date

#window_months: number of months for time-window prior to yearly time points

#output_dir: directory for output

##################################################################################

longitudinal_binning<-
  function(longitudinal_data_path,
           start_date_path,
           closest_avg,
           post_bin_years,
           pre_bin_years,
           window_months,
           output_dir){
    
    #Load packages
    library(data.table)
    library(dplyr)
    library(lubridate)
    
    ##################################################################################
    
    #Read in start date
    start_date<-
      fread(start_date_path,
            col.names = c("eid", "start_date"))
    
    ##################################################################################
    
    #Create dataframe with full paths to longitudinal data
    full_name<-
      list.files(path = longitudinal_data_path,
                 recursive = FALSE,
                 pattern = "^\\combined_codes.*\\.txt$",
                 full.names = TRUE)
    
    short_name<-
      list.files(path = longitudinal_data_path,
                 recursive = FALSE,
                 pattern = "^\\combined_codes.*\\.txt$",
                 full.names = FALSE)
    
    longitudinal_data<- as.data.frame(cbind(full_name, short_name))
    
    longitudinal_data<-
      longitudinal_data %>%
      filter(short_name != "archive") %>%
      mutate(short_name= as.character(short_name)) %>%
      mutate(short_name= substr(short_name, 1, nchar(short_name)-4)) %>%
      mutate(short_name= substr(short_name, 16, nchar(short_name))) %>%
      as.matrix()
    
    #Read in longitudinal data
    longitudinal_data_list<-
      lapply(longitudinal_data[,"full_name"], fread)
    
    names(longitudinal_data_list)<- longitudinal_data[,"short_name"]
    
    ##################################################################################
    
    #Annotate longitudinal data with start data & only keep people with a start date
    longitudinal_data_annotated<-
      lapply(longitudinal_data_list,
             function(x) inner_join(x,
                                    start_date,
                                    by="eid"))
    
    ##################################################################################
    ##For closest value to time point##
    ##################################################################################
    
    if(closest_avg == "closest") {
      
      for(j in names(longitudinal_data_list)){
        
        print(j)
        
        #Create dataframes for binned results
        longitudinal_data_binned_closest_df<-
          longitudinal_data_annotated[[j]] %>%
          select(-value, -date, -start_date, -trait)
        
        if(missing(post_bin_years)) {
          
          print("post_bin_years missing")
          
        } else {
          
          for(i in 1:post_bin_years){
            
            print(i)
            
            #Bin data into time windows- closest value to time point
            longitudinal_data_binned_closest_entry<-
              longitudinal_data_annotated[[j]] %>%
              dplyr::group_by(eid) %>% #Group by person 
              mutate(start_date= as.Date(start_date),
                     date= as.Date(date)) %>% #Convert all date columns to dates
              mutate(time_start = difftime(date, start_date , units = "days")) %>% #Calculate time between date and start date
              filter(time_start <= (365.25 * i) & time_start >= ((365.25 * i) - (30.4 * window_months))) %>% #Filter dates to be before yr time point & after X months prior to yr time point
              mutate(time_Xyrs = difftime((start_date %m+% years(i)) , date, units = "days")) %>% #Calculate time between date and yr time point
              arrange(abs(time_Xyrs)) %>% #Sort by time between date and yr time point
              mutate(min_time_Xyrs = (time_Xyrs)[1]) %>% #Identify date closest to yr time point
              filter(time_Xyrs == min_time_Xyrs)  #Only keep date closest to yr time point
            
            #Add binned results to dataframe for results
            longitudinal_data_binned_closest_df<-
              left_join(longitudinal_data_binned_closest_df,
                        longitudinal_data_binned_closest_entry %>%
                          select(eid, value), 
                        by= "eid") %>%
              unique()
            
            #Rename columns
            setnames(longitudinal_data_binned_closest_df,
                     old = "value",
                     new = paste(j, "_value_post_", i, "yrs", sep = ""),
                     skip_absent = TRUE)
            
          }
        }
        
        if(missing(pre_bin_years)) {
          
          print("pre_bin_years missing")
          
        } else {
          
          for(i in 1:pre_bin_years){
            
            print(i)
            
            #Bin data into time windows- closest value to time point
            longitudinal_data_binned_closest_entry<-
              longitudinal_data_annotated[[j]] %>%
              dplyr::group_by(eid) %>% #Group by person 
              mutate(start_date= as.Date(start_date),
                     date= as.Date(date)) %>% #Convert all date columns to dates
              mutate(time_start = difftime(date, start_date , units = "days")) %>% #Calculate time between date and start date
              filter(time_start >= (-365.25 * i) & time_start <= ((-365.25 * i) + (30.4 * window_months))) %>% #Filter dates to be before yr time point & after X months prior to yr time point
              mutate(time_Xyrs = difftime((start_date %m+% years(i)) , date, units = "days")) %>% #Calculate time between date and yr time point
              arrange(abs(time_Xyrs)) %>% #Sort by time between date and yr time point
              mutate(min_time_Xyrs = (time_Xyrs)[1]) %>% #Identify date closest to yr time point
              filter(time_Xyrs == min_time_Xyrs)  #Only keep date closest to yr time point
            
            #Add binned results to dataframe for results
            longitudinal_data_binned_closest_df<-
              left_join(longitudinal_data_binned_closest_df,
                        longitudinal_data_binned_closest_entry %>%
                          select(eid, value), 
                        by= "eid") %>%
              unique()
            
            #Rename columns
            setnames(longitudinal_data_binned_closest_df,
                     old = "value",
                     new = paste(j, "_value_pre_", i, "yrs", sep = ""),
                     skip_absent = TRUE)
            
          } 
        }
        
        #Write out results
        filename_closest<-
          paste(output_dir, j, "_", sep = "")
        
        if(missing(post_bin_years)) {
          
          print("post_bin_years missing")
          
        } else {
          
          filename_closest<-
            paste(filename_closest,
                  post_bin_years, "post_bin_years_", sep = "")
        }
        
        if(missing(pre_bin_years)) {
          
          print("pre_bin_years missing")
          
        } else {
          
          filename_closest<-
            paste(filename_closest,
                  pre_bin_years, "pre_bin_years", sep = "")
        }
        
        
        filename_closest<-
          paste(filename_closest,
                window_months, "window_months_", closest_avg, ".txt", sep = "")
        
        write.table(longitudinal_data_binned_closest_df,
                    file= filename_closest,
                    sep = "\t",
                    row.names = F,
                    quote = F)
        
        rm(filename_closest)
        rm(longitudinal_data_binned_closest_df)
        
      }
      
    }
    
    ##################################################################################
    ##For avg value in time window##
    ##################################################################################
    
    if(closest_avg == "avg") {
      
      for(j in names(longitudinal_data_list)){
        
        print(j)
        
        longitudinal_data_binned_avg_df<-
          longitudinal_data_annotated[[j]] %>%
          select(-value, -date, -start_date, -trait)
        
        if(missing(post_bin_years)) {
          
          print("post_bin_years missing")
          
        } else {
          
          for(i in 1:post_bin_years){
            
            print(i)
            
            #Bin data into time windows- avg value in time window
            longitudinal_data_binned_avg_entry<-
              longitudinal_data_annotated[[j]] %>%
              dplyr::group_by(eid) %>% #Group by person 
              mutate(start_date= as.Date(start_date),
                     date= as.Date(date)) %>% #Convert all date columns to dates
              mutate(time_start = difftime(date, start_date , units = "days")) %>% #Calculate time between date and start date
              filter(time_start <= (365.25 * i) & time_start >= ((365.25 * i) - (30.4 * window_months))) %>% #Filter dates to be before yr time point & after X months prior to yr time point
              mutate(avg_value= median(value)) %>% #Calculate median value in time window
              select(-value, 
                     -date) %>% #Remove value column
              unique() #Keep unique rows
            
            #Add binned results to dataframe for results
            longitudinal_data_binned_avg_df<-
              left_join(longitudinal_data_binned_avg_df,
                        longitudinal_data_binned_avg_entry %>%
                          select(eid, avg_value), 
                        by= "eid") %>%
              unique()
            
            #Rename columns
            setnames(longitudinal_data_binned_avg_df,
                     old = "avg_value",
                     new = paste(j, "_value_post_", i, "yrs", sep = ""),
                     skip_absent = TRUE)
            
          }
        }
        
        if(missing(pre_bin_years)) {
          
          print("pre_bin_years missing")
          
        } else {
          
          for(i in 1:pre_bin_years){
            
            print(i)
            
            #Bin data into time windows- avg value in time window
            longitudinal_data_binned_avg_entry<-
              longitudinal_data_annotated[[j]] %>%
              dplyr::group_by(eid) %>% #Group by person 
              mutate(start_date= as.Date(start_date),
                     date= as.Date(date)) %>% #Convert all date columns to dates
              mutate(time_start = difftime(date, start_date , units = "days")) %>% #Calculate time between date and start date
              filter(time_start >= (-365.25 * i) & time_start <= ((-365.25 * i) + (30.4 * window_months))) %>% #Filter dates to be before yr time point & after X months prior to yr time point
              mutate(avg_value= median(value)) %>% #Calculate median value in time window
              select(-value, 
                     -date) %>% #Remove value column
              unique() #Keep unique rows
            
            #Add binned results to dataframe for results
            longitudinal_data_binned_avg_df<-
              left_join(longitudinal_data_binned_avg_df,
                        longitudinal_data_binned_avg_entry %>%
                          select(eid, avg_value), 
                        by= "eid") %>%
              unique()
            
            #Rename columns
            setnames(longitudinal_data_binned_avg_df,
                     old = "avg_value",
                     new = paste(j, "_value_pre_", i, "yrs", sep = ""),
                     skip_absent = TRUE)
            
          } 
        }
        
        #Write out results
        filename_avg<-
          paste(output_dir, j, "_", sep = "")
        
        if(missing(post_bin_years)) {
          
          print("post_bin_years missing")
          
        } else {
          
          filename_avg<-
            paste(filename_avg,
                  post_bin_years, "post_bin_years_", sep = "")
        }
        
        if(missing(pre_bin_years)) {
          
          print("pre_bin_years missing")
          
        } else {
          
          filename_avg<-
            paste(filename_avg,
                  pre_bin_years, "pre_bin_years", sep = "")
        }
        
        
        filename_avg<-
          paste(filename_avg,
                window_months, "window_months_", closest_avg, ".txt", sep = "")
        
        write.table(longitudinal_data_binned_avg_df,
                    file= filename_avg,
                    sep = "\t",
                    row.names = F,
                    quote = F)
        
        rm(filename_avg)
        rm(longitudinal_data_binned_avg_df)
        
        
      }
    }
  }

##################################################################################   



