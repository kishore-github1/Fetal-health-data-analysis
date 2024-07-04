#Include the below libraries  
library(caTools)
library(caret)
library(e1071)
library(gpairs)
library(ggplot2)
library(ModelMetrics)
library(ROCR)
library(Metrics)
library(gplots)
library(dplyr)
library(tidyr)


#Indication
#class-1 :Normal
#class-2 :Suspect
#class-3 :Pathological


#Function to calculate mean
mean <- function(data){
  sum <- 0
  for(i in 1:length(data)){
    sum <- sum + data[i]
  }
  return(sum/length(data))
}

#Function to find Standard Deviation
std <- function(data,mean){
  sum <- 0
  for(i in 1:length(data)){
    sum <- sum + (data[i]-mean)^2
  }
  return(sqrt(sum/length(data)))
}

# Normalization Function
norm <- function(x, meanValue, stdValue){
  return (exp(-(x-meanValue)^2/(2*stdValue^2))/(stdValue*sqrt(2*pi)))
}

#Naive Bayes Classifier Function
predict <- function(training_data,test_x,continuous_columns){
  
  training_data.y0 <- subset(training_data,subset=training_data$fetal_health==1)
  training_data.y1 <- subset(training_data,subset=training_data$fetal_health==2)
  training_data.y2 <- subset(training_data,subset=training_data$fetal_health==3)
  
  cont_lkhd <- list()
  for(i in continuous_columns){
    cont_lkhd[[i]] <- list()
    mean <- mean(training_data.y0[[i]])
    std <- sd(training_data.y0[[i]],na.rm=TRUE)
    cont_lkhd[[i]][["0"]] <- list("mean" = mean,"std" = std)
    mean <- mean(training_data.y1[[i]])
    std <- sd(training_data.y1[[i]],na.rm=TRUE)
    cont_lkhd[[i]][["1"]] <- list("mean" = mean,"std" = std)
    mean <- mean(training_data.y2[[i]])
    std <- sd(training_data.y2[[i]],na.rm=TRUE)
    cont_lkhd[[i]][["2"]] <- list("mean" = mean,"std" = std)
  }
  
  predition_y = c()
  prior_target = c(length(which(training_data$fetal_health==1))/nrow(training_data),length(which(training_data$fetal_health==2))/nrow(training_data),length(which(training_data$fetal_health==3))/nrow(training_data))
  for(row in 1:nrow(test_x)){
    prob.0 <- prior_target[1]
    prob.1 <- prior_target[2]
    prob.2 <- prior_target[3]
    fetal_health <- 0
    for(i in continuous_columns){
      test_val <- test_x[row,i]
      prob.0 <- prob.0 * norm(test_val,cont_lkhd[[i]][["0"]][["mean"]],cont_lkhd[[i]][["0"]][["std"]])
      prob.1 <- prob.1 * norm(test_val,cont_lkhd[[i]][["1"]][["mean"]],cont_lkhd[[i]][["1"]][["std"]])
      prob.2 <- prob.2 * norm(test_val,cont_lkhd[[i]][["2"]][["mean"]],cont_lkhd[[i]][["2"]][["std"]])
    }
    
    if (prob.0 > prob.1 & prob.0 > prob.2) {
      fetal_health <- 1
    } else if (prob.1 > prob.2) {
      fetal_health <- 2
    } else {
      fetal_health <- 3
    }
    
    predition_y[row] <- fetal_health
  }
  return(predition_y)
}

#k-fold cross validation function
KFoldCrossValidation <- function(dataset,k) {
  sprintf("Cross Validation with K = %d",k)
  len = nrow(dataset)
  accuracy = c()
  for(i in 1:k){
    from <- 1+as.integer(len/k)*(i-1)
    to <- as.integer(len/k)*i
    testing_data <- dataset[from:to,]
    training_data <- rbind(dataset[0:(from-1),],dataset[(to+1):len,])
    
    test_x <- subset(testing_data,select = -fetal_health)
    test_y <- subset(testing_data,select = fetal_health)
    
    result <- data.frame(test_y)
    result$predicted_target <- predict(training_data,test_x,continuous_columns)
    
    accuracy[i] <- length(which(result$fetal_health==result$predicted_target))/length(result$fetal_health)*100
  }
  return(mean(accuracy))
}


#importing dataset
dataset = read.csv("C:\\Users\\mahan\\Downloads\\IDA Project\\fetal_health.csv") #read the dataset
dataset = na.omit(dataset)
str(dataset)
summary(dataset)

# Handle Missing Values
your_data <- na.omit(dataset)
colSums(is.na(your_data))

# Check for missing values in the entire dataset
if(any(is.na(your_data)) == FALSE){
  print('No missing values')
}

# Assuming 'your_data' is your data frame
num_duplicates <- sum(duplicated(your_data))
cat("Number of Duplicate Rows: ", num_duplicates, "\n")

# Display duplicated rows
duplicated_rows <- your_data[duplicated(your_data), ]
print(duplicated_rows)


# Display all occurrences of duplicated rows
all_duplicates <- your_data[duplicated(your_data) | duplicated(your_data, fromLast = TRUE), ]
print(all_duplicates)

# Assuming 'your_data' is your data frame
your_data <- unique(your_data)

# Remove duplicates and create a new data frame
unique_data <- unique(your_data)
print(unique_data)

## Outlier detection and removal 
column_names <-colnames(unique_data)
threshold <- 1.5
lower_bounds <- c() 
upper_bounds <- c() 

for (column in column_names) { 
  column_iqr <- IQR(your_data[[column]]) 
  column_lower_bound <- quantile(your_data[[column]], 0.25) - threshold * column_iqr 
  column_upper_bound <- quantile(your_data[[column]], 0.75) + threshold * column_iqr 
  
  lower_bounds <- c(lower_bounds, column_lower_bound) 
  upper_bounds <- c(upper_bounds, column_upper_bound) 
}

for (i in 1:length(column_names)) { 
  column <- column_names[i] 
  lower_bound <- lower_bounds[i] 
  upper_bound <- upper_bounds[i] 
  
  your_data <- your_data[your_data[[column]] >= lower_bound & your_data[[column]] <= upper_bound, ] 
} 

count <- table(dataset$fetal_health)
count_percent<-prop.table(count)*100
#class-1 :Normal
#class-2 :Suspect
#class-3 :Pathological
pie(count,labels =paste0("Class-",names(count),": ",round(count_percent,2),"%"),main = "Fetal Health Chart")

#all the columns in the dataset are continuous
continuous_columns = c("baseline.value",                                       
                       "accelerations",                                         
                       "fetal_movement",                                        
                       "uterine_contractions",                                  
                       "light_decelerations",                                   
                       "prolongued_decelerations",                              
                       "abnormal_short_term_variability",                       
                       "mean_value_of_short_term_variability",                  
                       "percentage_of_time_with_abnormal_long_term_variability",
                       "mean_value_of_long_term_variability",                 
                       "histogram_width",                                       
                       "histogram_min",                                         
                       "histogram_max",                                         
                       "histogram_number_of_peaks",                            
                       "histogram_number_of_zeroes",                           
                       "histogram_mode",                                        
                       "histogram_mean",                                        
                       "histogram_median",                                      
                       "histogram_variance",                                    
                       "histogram_tendency"                                  
)#severe_decelerations is ignored since it has low variance

#scatterplot function
ScatterPlot <- function(x,y){
  ggplot(dataset, aes(.data[[x]], .data[[y]])) + geom_point(aes(color = fetal_health)) + 
    scale_x_continuous(x)+
    scale_y_continuous(y)+
    theme_bw() + labs(title=paste("Scatterplot: ",x,"V/S",y))
}

ScatterPlot("baseline.value","accelerations")
ScatterPlot("baseline.value","fetal_movement")
ScatterPlot("baseline.value","uterine_contractions")
ScatterPlot("baseline.value","light_decelerations")
ScatterPlot("baseline.value","prolongued_decelerations")
ScatterPlot("baseline.value","abnormal_short_term_variability")
ScatterPlot("baseline.value","mean_value_of_short_term_variability")
ScatterPlot("baseline.value","percentage_of_time_with_abnormal_long_term_variability")
ScatterPlot("baseline.value","mean_value_of_long_term_variability")
ScatterPlot("baseline.value","histogram_width")
ScatterPlot("baseline.value","histogram_min")
ScatterPlot("baseline.value","histogram_max")
ScatterPlot("baseline.value","histogram_number_of_peaks")
ScatterPlot("baseline.value","histogram_number_of_zeroes")
ScatterPlot("baseline.value","histogram_mode")
ScatterPlot("baseline.value","histogram_mean")
ScatterPlot("baseline.value","histogram_median")
ScatterPlot("baseline.value","histogram_variance")
ScatterPlot("baseline.value","histogram_tendency" )

#K-fold cross validation 
avgAccuracy = c()
least.k = 3
max.k= 15
max_accuracy = 0
maxk=0
for(k in least.k:max.k){
  accuracy <- KFoldCrossValidation(dataset,k)
  if(accuracy>max_accuracy){
    max_accuracy=accuracy
    maxk=k
  }
  avgAccuracy<- append(avgAccuracy,accuracy)
}

print(avgAccuracy)
bestKvalue = which.max(avgAccuracy)+least.k-1
print(bestKvalue)

print(paste0("Highest Accuracy: ", max_accuracy, "for k= ", maxk))
sprintf("So the acuracy with best k = %d value for k fold cross validation is %f",bestKvalue,avgAccuracy[[bestKvalue-least.k+1]])

split <- createDataPartition(dataset$fetal_health,p=0.75,list=FALSE)
td.size <- 0.75*nrow(dataset)
training_data <- dataset[1:td.size,]
testing_data <- dataset[td.size:nrow(dataset),]
test_x <- subset(testing_data, select = -fetal_health)
test_y <- subset(testing_data, select = fetal_health)

#contingency table from an training data set comprising the prior and posterior probabilities.
for(i in continuous_columns){
  prior.ct <- prop.table(table(training_data$fetal_health,training_data[[i]],dnn = list("fetal_health",i)))
  cat("Prior Probabilties of fealth_health and",i,"\n")
  print(prior.ct)
  cat("\n")
  post.ct <- prop.table(table(training_data$fetal_health,training_data[[i]],dnn = list("fetal_health",i)),margin = 2)
  cat("Posterior Probabilties of fetal_health and",i,"\n")
  print(post.ct)
  cat("\n")
}

result <- data.frame(test_y)
result$predicted_target <- predict(training_data,test_x,continuous_columns)
#the predicted labels are stored in a csv file
write.csv(result,"C:\\Users\\mahan\\Downloads\\IDA Project\\result.csv", row.names = FALSE)

# Accuracy 
accuracy <- length(which(result$fetal_health==result$predicted_target))/length(result$fetal_health)*100
sprintf("accuracy: %f",accuracy)

# Confusion Matrix
confusion_matrix <- table(result$fetal_health,result$predicted_target,dnn=names(result))
print(confusion_matrix)
cat("\n")

#Function to find precision,recall,F1-Score
find_PRF1<- function(true_labels, predicted_labels, class_label) {
  # Ensure the inputs have the same length
  if (length(true_labels) != length(predicted_labels)) {
    stop("Input vectors must have the same length.")
  }
  
  # Calculate the number of true positives for the specified class
  true_positives <- sum(true_labels == class_label & predicted_labels == class_label)
  
  # Calculate the number of predicted positives for the specified class
  predicted_positives <- sum(predicted_labels == class_label)
  
  # Calculate the number of actual positives for the specified class
  actual_positives <- sum(true_labels == class_label)
  
  
  # Calculate precision for the specified class
  precision <- true_positives / predicted_positives
  
  # Calculate recall for the specified class
  recall <- true_positives / actual_positives
  
  f1_score_value <- 2 * (precision * recall) / (precision + recall)
  
  cat("Precision for Class", class_label, ":", precision, "\n")
  cat("Recall for Class", class_label, ":", recall, "\n")
  cat("F1 score for Class", class_label, ":", f1_score_value, "\n\n")
  
}

for (class_label in 1:3) {
  find_result <- find_PRF1(result$fetal_health,result$predicted_target, class_label)
}

