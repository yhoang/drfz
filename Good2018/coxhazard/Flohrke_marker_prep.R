# function that checks the differents between two marker lists and then removes that
# difference from dataframe read ins (columns)

# define output path and filename
outfile.training <- "/home/felix/AG_Baumgrass/Scripts/Pri/Pri_good_established/Rds/BCR/BCR_Training_func3_quadrant_absRange_cof0.2.rds"
outfile.validation <- "/home/felix/AG_Baumgrass/Scripts/Pri/Pri_good_established/Rds/BCR/BCR_Validation_func3_quadrant_absRange_cof0.2.rds"

# list names of all markers 
markers.full <- scan("/home/felix/AG_Baumgrass/Data/Good/Marker_combinations/columns_full.txt", what="character", sep="\n")

# list names of markers to extract (by removing difference from full)
markers.func <- scan("/home/felix/AG_Baumgrass/Data/Good/Marker_combinations/columns_func_plus3.txt", what="character", sep="\n")

# list names of markers that are in full but not func (used to remove them from dataframe)
markers.notfunc <- markers.full[!(markers.full %in% markers.func)]

# comparing marker inputs
printf("markers full: %s", length(markers.full))
print(markers.full)
print("######################################")

printf("markers func: %s", length(markers.func))
print(markers.func)
print("######################################")

printf("markers not in func: %s", length(markers.notfunc))
print(markers.notfunc)
print("######################################")

# Training and validation sets from which markers not in func are to be removed
df.training.func <- readRDS("/home/felix/AG_Baumgrass/Scripts/Pri/Pri_good_established/Rds/BCR/BCR_Training_full_quadrant_absRange_cof0.2.rds")
df.validation.func <- readRDS("/home/felix/AG_Baumgrass/Scripts/Pri/Pri_good_established/Rds/BCR/BCR_Validation_full_quadrant_absRange_cof0.2.rds")

# saving rownames
rownames.training <- row.names(df.training.func)
rownames.validation <- row.names(df.validation.func)

##### bug checking #####
#trainingold <- df.training.func
#trainingnew <- readRDS("/home/felix/AG_Baumgrass/Scripts/Pri/Pri_good_established/Rds/BCR/BCR_Training_full_quadrant_absRange_cof0.2.rds")
#print(setdiff(names(trainingold), names(trainingnew)))
#print(grep("CD123", trainingold))
#print(grep("CD123", trainingnew))
#print(names(trainingold[1]))
#print(names(trainingnew[1]))
#print(typeof(trainingold))
#print(typeof(trainingnew))
#print(trainingold[,1])
#print(trainingnew[,1])
#print(names(trainingold[1]))
#print(names(trainingnew[1]))
#diff <- setdiff(names(df.training.func), names(df.validation.func))
#print(length(diff))
#print(diff)
#stop()
##### bug checking #####

# saving marker names of dataframes
names.df.training <- names(df.training.func)
names.df.validation <- names(df.validation.func)

# saving new df from read ins to change
df.training.funconly <- as.data.frame(df.training.func)
df.validation.funconly <- as.data.frame(df.validation.func)

print("#### Working on training set ####")
### training
for (i in 1:length(markers.notfunc)) {
    # index pos of markers to be removed
    #print(names.df.training[1:20])
    printf("Marker that is selected to remove: %s", markers.notfunc[i])
    ind.train <- grep(markers.notfunc[i], names.df.training)
    #ind.val <- grep(markers.notfunc[i], names.df.validation)
    
    # only execute if to removed markers are still in the dataframe
    if (length(ind.train) != 0){
        df.training.funconly <- df.training.funconly[,-ind.train]
        #df.validation.funconly <- df.validation.funconly[,-ind.val]
        names.df.training <- names(df.training.funconly)
        #names.df.validation <- names(df.validation.funconly)

        print(ncol(df.training.funconly))
        print(nrow(df.training.funconly))
        #print(ncol(df.validation.funconly))
        #print(nrow(df.validation.funconly))

    
    }  else {
        print("marker not found")
    }

}
print("#### Working on validation set ####")
### validation
for (i in 1:length(markers.notfunc)) {
    # index pos of markers to be removed
    #print(names.df.training[1:20])
    printf("Marker that is selected to remove: %s", markers.notfunc[i])
    #ind.train <- grep(markers.notfunc[i], names.df.training)
    ind.val <- grep(markers.notfunc[i], names.df.validation)
    #print("validation cols with CD123")
    #print(ind.val[1:100])
    
    
    # only execute if to removed markers are still in the dataframe
    if (length(ind.val) != 0){
        #df.training.funconly <- df.training.funconly[,-ind.train]
        df.validation.funconly <- df.validation.funconly[,-ind.val]
        #names.df.training <- names(df.training.funconly)
        names.df.validation <- names(df.validation.funconly)

        #print(ncol(df.training.funconly))
        #print(nrow(df.training.funconly))
        print(ncol(df.validation.funconly))
        print(nrow(df.validation.funconly))

    
    } else {
        print("marker not found")
    }

}

# showing new dimensions of dataframes
print("Removal succesful")
print("new dimensions of df.training and df.validation: ")
print(ncol(df.training.funconly))
print(nrow(df.training.funconly))
print(ncol(df.validation.funconly))
print(nrow(df.validation.funconly))

# adding rownames
row.names(df.training.funconly) <- rownames.training
row.names(df.validation.funconly) <- rownames.validation

# saving modified training and validation dataframes

saveRDS(df.training.funconly,outfile.training)
saveRDS(df.validation.funconly,outfile.validation)
print("Saved new df as RDs")

