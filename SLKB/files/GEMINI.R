BioCManagerLibraries <- c("gemini", 'ggplot2')
requiredLibraries <- c(BioCManagerLibraries)

for (packageName in requiredLibraries){
  if (!is.element(packageName, installed.packages()[,1])){
    print(paste("Installing package: ", packageName))
    if (packageName %in% BioCManagerLibraries) {
      BiocManager::install(packageName, dependencies = TRUE, INSTALL_opts = '--no-lock') #, INSTALL_opts = '--no-lock') #
    }
    else {
      install.packages(packageName, dependencies = TRUE, type="binary")#, INSTALL_opts = '--no-lock')
      #install.packages("table1", dependencies = TRUE, INSTALL_opts = '--no-lock')#type="binary")#
      
    }
    
  } 
  
  suppressMessages(library(packageName, character.only = TRUE))
  print(paste("Loaded package: ", packageName))
  
}

library(parallel)
detectCores()

################### GEMINI ###################

print('Starting GEMINI Scoring in R...')

# get the current working directory
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2){
    print('Wrong args')
}

# arg 2 is the wd location
setwd(file.path(getwd(), args[2]))

# setwd(file.path("/users/PAS1376/bg12/SyntheticLethality - NewDB/Python_Clean/GEMINI_Files/horlbeck_data/JURKAT"))
# setwd(file.path("/users/PAS1376/bg12/SyntheticLethality - NewDB/Python_Clean/GEMINI_Files/parrish_data/HELA"))

# print the location
print(getwd())

# # set save loc
save_loc <- getwd()
dir.create(save_loc, recursive = TRUE)
save_loc <- file.path(save_loc, "GEMINI_Scores.csv")
    
# get sequences
curr_sequence_ref <- read.csv(file = "sequences.csv")
curr_sequence_ref <- curr_sequence_ref[!duplicated(curr_sequence_ref[, 1]),]

# get counts
curr_counts <- read.csv(file = 'counts.csv')
# remove dups
dups <- duplicated(curr_counts[, 1])
curr_counts <- curr_counts[!dups, ]
rownames(curr_counts) <- curr_counts[, 1]
curr_counts <- curr_counts[, -1]
    
# get annotations
curr_annotations <- read.csv(file = 'guide_gene_annotation.csv')
# remove from annotations as well, if they exist
curr_annotations <- curr_annotations[!dups, ]
colnames(curr_annotations)[1] <- "rowname"

# prepare replicate annotations
curr_replicate_annotations <- data.frame(matrix(ncol = 3, nrow = length(colnames(curr_counts))))
colnames(curr_replicate_annotations) <- c('colname', 'samplename', 'replicate')
curr_replicate_annotations$colname <- colnames(curr_counts)
curr_replicate_annotations$samplename <- sapply(strsplit(colnames(curr_counts), split = '_', fixed = TRUE), '[', 1)
curr_replicate_annotations$replicate <- sapply(strsplit(colnames(curr_counts), split = '_', fixed = TRUE), '[', 2)

ETP_locs <- grep("T0_", colnames(curr_counts))
LTP_locs <- grep("TEnd_", colnames(curr_counts))
    
# create input
gemini_input <- gemini_create_input(counts.matrix = curr_counts,
                                    sample.replicate.annotation = curr_replicate_annotations,
                                    guide.annotation = curr_annotations,
                                    ETP.column = ETP_locs,
                                    LTP.column = LTP_locs,
                                    gene.column.names = c("Gene.1", "Gene.2"),
                                    sample.column.name = "samplename",
                                    samplesAreColumns = TRUE,
                                    verbose = TRUE)
    
    
## follow the pipeline, apply preprocessing
gemini_input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                      CONSTANT = 32)
# initialize gemini model
Model <- gemini_initialize(Input = gemini_input, 
                          nc_gene = "CONTROL",
                          pattern_join = ';',
                          pattern_split = ';', 
                          cores = detectCores(),
                          verbose = TRUE)

#Model$nc_gene <- "CONTROL"

#save.image("horlbeck_jurkat_w_horlbeck_data_06_16_6:52pm.RData")
#n_iterations <- 1
# run inference
Model %<>% gemini_inference(cores = detectCores(),
                            n_iterations = 20,
                            force_results = TRUE,
                            verbose = TRUE,
                            save_iterations = TRUE)

print('Finished running inference')

#save.image("horlbeck_jurkat_w_horlbeck_data_ran.RData")

# save plot
ggplot2::ggsave(filename = file.path(getwd(), "model_mae.png"), plot = gemini_plot_mae(Model))


# nc_pairs <- grep("CONTROL", rownames(Model$s), value = TRUE)
# 
# gemini_scores <- gemini_score_UPDATED(Model,
#                               nc_pairs = c("CONTROL;CONTROL"))

gemini_scores <- gemini_score(Model)

# create dataframe of the scores
gemini_scores <- data.frame(Strong = gemini_scores$strong,
           SensitiveLethality = gemini_scores$sensitive_lethality,
           SensitiveRecovery = gemini_scores$sensitive_recovery)
           
colnames(gemini_scores) <- c('Strong', 'SensitiveLethality', 'SensitiveRecovery')

# save environment and scores
write.csv(gemini_scores, file = save_loc)

save.image(file = file.path(getwd(), "env.RData"))

print('Done!')
