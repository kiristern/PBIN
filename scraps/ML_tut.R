#https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial

otu_table <- read.table("/Users/kiristern/Downloads/otu_table_RF_tutorial.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")  
metadata <- read.table("/Users/kiristern/Downloads/metadata_RF_tutorial.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

dim(otu_table) #1000 rows (each corresponding to a different OTU), 40 samples
dim(metadata) #2 cols of metadata, 40 samples

str(metadata)
summary(metadata)

#### PreProcessing ####
#removing rare features: A basic way to decide which features are unlikely to be informative is if they are non-zero in only a small number of samples. 
#It's a good idea to take a look at this distribution for all features to help choose a reasonable cut-off.
otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
  #most OTUs are rare and found in only a few samples. These OTUs are less likely to help our model so we can throw them out. 
  #Typically researchers discard OTUs that are zero in greater than 75-90% of samples although these cut-offs are somewhat arbitrary. 
  #This cut-off could be optimized if you had an independent dataset or partition of your data that you could use as a validation set.

#removes features that are non-zero less than a specified proportion of the time:
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

#Based on the histogram, remove OTUs that have non-zero values in <= 20% of samples
otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.2)

#check how many OTUs are left
dim(otu_table_rare_removed) #647 OTUs are remaining, which means that 35% of the OTUs were excluded.

#To re-normalize_table so that each sample's column sums to 100 (you can use the ColSums function afterwards for a sanity check):
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
colSums(otu_table_rare_removed_norm)

#might also be interested in excluding OTUs that don't vary much across samples. This usually isn't different from excluding OTUs based on their number of non-zero values. 
#However, this can be useful for other types of feature tables. You may be interested in trying out the nearZeroVar() function from the caret R package to perform this filtering.



#### Transforming data ####

#standardize data by subtracting each sample's mean (center) and dividing by the stdev (scale). Ie. convert to Z-score
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  
#another method: take inverse hyperbolic since and then to mean center by sample
otu_table_asinh_mean_centred <- scale(asinh(otu_table), center=TRUE, scale=FALSE)  
  #think carefully about what transformation is most likely to help your model's performance without detracting too much from the interpretability




#### Runing Model ####
#input tables into the correct format (metadata col is last col of each input df)
  #to prep input tables for classification of state
otu_table_scaled_state <- data.frame(t(otu_table_scaled))  
otu_table_scaled_state$state <- metadata[rownames(otu_table_scaled_state), "state"] 

#prp input tbales for regresison of inflammation score (IS)
otu_table_scaled_IS <- data.frame(t(otu_table_scaled))  
otu_table_scaled_IS$IS <- metadata[rownames(otu_table_scaled_IS), "IS"]

#the 2 parameters for a RF model are the number of trees in the forest (ntree) and the number of features randomly sampled at each node in a tree (mtry). 
#The more trees you run in your forest the better the model will converge. 
#Below I used 501, which is similar to the default, but in practice you may want to use something like 10,001 trees for a robust model (depending on the computational time). 
#Note that I usually choose odd numbers of trees to ensure there are never any ties for binary classification models. 
#Unless you have a reason to change mtry beforehand (or can optimize it with an independent partition of data) it's better to use the default values.

set.seed(151)  

#Run RF to classify inflamed and control samples:
RF_state_classify <- randomForest( x=otu_table_scaled_state[,1:(ncol(otu_table_scaled_state)-1)] , y=otu_table_scaled_state[ , ncol(otu_table_scaled_state)] , ntree=501, importance=TRUE, proximities=TRUE )

#Run RF to regress OTUs against inflammation score (IS):
RF_IS_regress <- randomForest( x=otu_table_scaled_IS[,1:(ncol(otu_table_scaled_IS)-1)] , y=otu_table_scaled_IS[ , ncol(otu_table_scaled_IS)] , ntree=501, importance=TRUE, proximities=TRUE )  


#### Assessing Model Fit ####
RF_state_classify
  #The out-of-bag error rate is 0%, which is perfect! Note that this summary also gives us the mtry parameter, which was 25 in this case.
RF_IS_regress
  #much higher mtry parameter and that instead of out-of-bag error and a confusion matrix, we can use the mean of squared residuals and the % variance explained as performance metrics. 
  #These metrics suggest that this model is also performing extremely well although it isn't quite as clear as for the classification model.




#### Permutation Test ####
#Test model's performance metric to see if it's more extreme than expected by chance based on a permutation test is to test model significance
#To run these significance tests with 1000 permutations
RF_state_classify_sig <- rf.significance( x=RF_state_classify ,  xdata=otu_table_scaled_state[,1:(ncol(otu_table_scaled_state)-1)] , nperm=1000 , ntree=501 )  
RF_IS_regress_sig <- rf.significance( x=RF_IS_regress ,  xdata=otu_table_scaled_IS[,1:(ncol(otu_table_scaled_IS)-1)] , nperm=1000 , ntree=501 )  




#### Acuracy estimated by cross-validation ####
#One method to estimate model performance is to systematically partition the data into training and test sets and repeatedly see how the model performs (cross-validation). 
#There are many different ways to partition the data, but the simplest is leave-one-out cross-validation. The model is trained n times, where n is the number of samples, and one sample is left 
#out each time for testing. This provides estimates of model performance, which tend to be quite similar to the internal measures (out-of-bag error and mean of squared residuals) described above. 
#We will run this cross-validation using the caret R package.

#First step is defining the parameters we want to use for training
fit_control <- trainControl(method="LOOCV") #specifying leave-one-out cross-validation

#run the leave-one-out cross-validation, note the same ntree and mtry parameters are being used
RF_state_classify_loocv <- train( otu_table_scaled_state[,1:(ncol(otu_table_scaled_state)-1)] , 
                                  y=otu_table_scaled_state[, ncol(otu_table_scaled_state)] , 
                                  method="rf", 
                                  ntree=501 , 
                                  tuneGrid=data.frame( mtry=25 ) , 
                                  trControl=fit_control )
RF_IS_regress_loocv <- train( otu_table_scaled_IS[,1:(ncol(otu_table_scaled_IS)-1)] , 
                              y=otu_table_scaled_IS[, ncol(otu_table_scaled_IS)] , 
                              method="rf", 
                              ntree=501 , 
                              tuneGrid=data.frame( mtry=215 ) , 
                              trControl=fit_control )
#check performance metrics:
RF_state_classify_loocv$results
RF_IS_regress_loocv$results
  #Both the Accuracy and Kappa metrics are 100%, which should make us extremely confident in our classification model (you won't see models perform this well usually!). 
  #Also, the R-squared of the regression model is very similar to what we previously found based on the internal testing, which is typical since the internal RF performance metrics should not be biased.


#### ID Important Features ####







