
########################################################################################################################
############################ Segment cells using the Imaris batch function and customized pipeline ########################
########################################################################################################################

# Create model to be used for Imaris segmentation using a representative image of the dataset
# Choose channel for segmentation
# Set Parameters : -smoothing surface details = 1
#                  -Background subtraction for local contrast = 4
#                  -Choose threshold manually
#                 - enable separation (region growing) using a seed points simager of 5 and a morphological split

# Filters the segmented objects based on selected parameters, for example:
# - Area 
# - Intensity of selected channels
# - Sphericity
# - ...

# save model and apply on all samples using batch

########################################################################################################################
############################ Merge exported parameters of interest for segmented cells into a single csv file ########################
########################################################################################################################


# Before running this script: 1: create one folder per sample, and 2. organize imaris output into the respective sample folder
# create one folder in which all the sample folders are located (here "Excel")
# The script will go into each sample folder, then into each segmented cell type, 
# then find the csv files that contain certain values such as "Area", "Postion X", etc. as detailed below
# and merge these csv into one csv file that is saved into the  sample folder.
 
library(dplyr)

Folders <- list.files(path = "Excel/") # this is the folder that contains all the sample folder


for (l2 in Folders){

  Tumour_folder <-l2
  
  folders_inTumour_folder <- list.files(path = paste0("Excel/",Tumour_folder))
  
  for (l3 in folders_inTumour_folder){
    print("New_loop")
    cell_folder <- l3
    print(cell_folder)
    cell_type <- strsplit(as.character(cell_folder), split="_")[[1]][1]
    files <- list.files(path = paste0("Excel/",Tumour_folder,"/",cell_folder))
    file_area <- grep(pattern = "Area.csv",x= files,value = T)
    file_Position_X <- grep(pattern = "Position_X.csv",x= files,value = T)
    file_Position_Y <- grep(pattern = "Position_Y.csv",x= files,value = T)
    file_Sphericity <- grep(pattern = "Circularity.csv",x= files,value = T)
    Intensity_Mean_Ch3 <- grep(pattern = "Intensity_Mean_Ch=3",x= files,value = T)
    Intensity_Mean_Ch4 <- grep(pattern = "Intensity_Mean_Ch=4",x= files,value = T)
    files_list <-  as.list(c(file_area,file_Position_X,file_Position_Y,file_Sphericity,Intensity_Mean_Ch3,Intensity_Mean_Ch4))
    print(files_list)
    
    tables_list <- lapply(files_list,function(l4){
      
      path <- paste0("Excel/",Tumour_folder,"/",cell_folder,"/",l4)
      Test <- read.csv(path,header = F,stringsAsFactors = F)
      colnames(Test) <- Test[3,]
      colnames(Test)[1] <- Test[1,1] # have to give first the name of the "header" to the first column, because otherwise the CHannel is not included in the intensity max
      colnames(Test)[2:ncol(Test)] <- Test[3,2:ncol(Test)]
      Test_sub <-Test[4:nrow(Test),]
      Test_sub_sub  <-select(Test_sub, colnames(Test_sub)[ 1:(which(colnames(Test_sub) == "Unit", arr.ind=TRUE)-1) ] ,"ID") # select param of intrest: columns unil "Unit column, and the one called ID to merge
      return(Test_sub_sub)
    })
    
    merged_table <- Reduce(function(x, y){merge(x, y, by ="ID")},tables_list)
    #head(tables_list[[3]])
    
    colnames(merged_table) <- gsub('=','_',colnames(merged_table))

    
    write.csv(merged_table,paste0("Excel/",Tumour_folder,"/",cell_type,".csv"),row.names = F)
    
  }
}

