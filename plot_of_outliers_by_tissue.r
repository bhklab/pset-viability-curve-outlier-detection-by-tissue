# this program finds the dose response curve outliers by tissue and prints cellline-drug pairs 
# that have large variances in viability when plotted

#specify the dataset below
screen <- PSet_gCSI2019

sensitivity_data <- sensitivityRaw(screen)
my_drugs <- screen@treatment$treatmentid
tissue_names <- screen@sample$tissueid
all_tissues <- c(unique(tissue_names))
outliers_all_tissues <- data.frame()

for (tissue in all_tissues) {
  message(paste0("Tissue: ", tissue))
  tissue_sample <- screen@sample[screen@sample$tissueid == tissue, ]
  my_celllines <- tissue_sample$sampleid

  all_variances <- c()
  replicated_drugs <- c()
  replicated_celllines <- c()
  
  for (my_cellline in my_celllines) {
    message(paste0("Processing cell line: ", my_cellline, " in ", tissue))
    for (my_drug in my_drugs) {
      #the pattern is determined by the way replicates are named in the PSet (most start with cellline_drug_)
      #the pattern needs to be determined before using the program (some replicates in PSets start with drug_cellline_ or something else)
      pattern <- paste0("^", my_cellline, "_", my_drug, "_")
      
      #interpreting regex characters as literal for ones found in replicate names
      pattern <- gsub("\\[", "\\\\[", pattern)
      pattern <- gsub("\\]", "\\\\]", pattern)
      pattern <- gsub("\\.", "\\\\.", pattern)
      pattern <- gsub("\\(", "\\\\(", pattern)
      pattern <- gsub("\\)", "\\\\)", pattern)
      
      experiment_ids <- dimnames(sensitivity_data)[[1]]
      
      matches <- grepl(pattern, experiment_ids)
      matching_experiment_ids <- experiment_ids[matches]
      
      if (length(matching_experiment_ids) > 1) {
        replicates_dose_viability_list <- list()
        
        for (experiment_id in matching_experiment_ids) {
          current_experiment_df <- as.data.frame(sensitivity_data[experiment_id, , ])
          replicates_dose_viability_list[[experiment_id]] <- current_experiment_df
        }
        
        if (all(sapply(replicates_dose_viability_list, is.null))) {
          message(paste0("No data in matching experiment ids."))
          next
        }
        
        #puts all the experiments together for each drug and aggregate them
        matches_dose_viability <- na.omit(do.call(rbind, replicates_dose_viability_list))

        #currently doses are rounded to 4 decimals, which can be changed here
        matches_dose_viability$Dose <- round(matches_dose_viability$Dose, 4)

        if (nrow(matches_dose_viability) > 0 && length(unique(matches_dose_viability$Viability[!is.na(matches_dose_viability$Viability)])) >= 2) {
          statistical_variance <- aggregate(Viability ~ Dose, data = matches_dose_viability, FUN = var)
          variance_df <- statistical_variance
          #this is renamed because the array from sensitivity raw has the column named viability
          colnames(variance_df)[colnames(variance_df) == "Viability"] <- "Variance"
          average_variance <- mean(variance_df$Variance, na.rm = TRUE)
          
          #checks for variance data not overall data (so doses need to be identical)
          if (is.finite(average_variance)) {
            all_variances <- c(all_variances, average_variance)
            replicated_drugs <- c(replicated_drugs, my_drug)
            replicated_celllines <- c(replicated_celllines, my_cellline)
          } else {
            message(paste0("Variance for '", pattern, "' is NA/NaN."))
          }
        } else {
          message(paste0(pattern, " was not processed because there is not enough data for it."))
        }
      }
    }
  }

  if (length(all_variances) > 0) {
    average_variance_df <- data.frame(
      all_variances = all_variances,
      replicated_drugs = replicated_drugs,
      replicated_celllines = replicated_celllines
    )
    
    mean_variance <- mean(all_variances)
    sd_variance <- sd(all_variances)
    
    if (!(is.na(sd_variance)) && (sd_variance != 0)) {
      normalized_variances <- (average_variance_df$all_variances - mean_variance) / sd_variance
    } else {
      print(("Standard deviation is 0/NA for this tissue."))
      normalized_variances <- rep(0, nrow(average_variance_df))
    }
    
    average_variance_df$normalized_variances <- normalized_variances

    #change the criteria for an outlier here (right now, an outlier is at least 3 SD)
    outliers_df <- average_variance_df[average_variance_df$normalized_variances > 3, ]
    
    print(outliers_df)
    outliers_all_tissues <- rbind(outliers_all_tissues, outliers_df)
    print(paste0("# of assays with valid replicates in ", tissue, ": ", length(replicated_drugs)))
    print(paste0("Number of outliers: ", nrow(outliers_df)))
    
    hist(
      normalized_variances,
      main = paste0("Histogram of norm. variances in viability for assays in ", tissue, " tissue"),
      xlab = "# of Standard Deviations", 
      col = "mediumpurple1",
      border = "black",
      breaks = 30 
    )
    
    if (nrow(outliers_df) != 0) {
      for (i in 1:nrow(outliers_df)) {
        outlier_drug <- outliers_df$replicated_drugs[i]
        outlier_cellline <- outliers_df$replicated_celllines[i]
        drugDoseResponseCurve(screen, drug = outlier_drug, cellline = outlier_cellline, summarize.replicates = FALSE)
      }
    }
  }
}

#View(outliers_all_tissues) after the program has completed running to see outliers for the entire dataset


