
library(glmnet)


#################################################################################################
#
# Evaluate models for estimating rosette dry mass from top-view images
#
#################################################################################################


# Import phenotypic data extracted from ImageJ analysis of rosette photos
#========================================================================
datot <- read.csv("16d_morpho_data.csv", header=T, sep=",") # file containing rosette shape descriptors on all individuals and dry mass measurement on the traning population
datot $idPot <- as.factor(as.character(datot $idPot))
datot $Day <- as.factor(as.character(datot $Day))
dat <- datot[,c(1,11:15,19:21)]
dat <- na.omit(dat)
dat <- droplevels(dat)
counttot <- dat


# Perform and compare models
#--------------------------------------------------
temp0 <- data.frame(nit=NA,ntrain=NA,lm_accuracy=NA,quaf_accuracy=NA,LASSO_accuracy=NA,RIDGE_accuracy=NA)

for(k in c(10,20,30,seq(50,250,50)))
{
  lm_accuracy <- NA
  quad_accuracy <- NA
  LASSO_accuracy <- NA
  RIDGE_accuracy <- NA
  
  for(n in 1:100)
  {
    counttot$set <- NA
    trainsample <- sample(x=levels(counttot$idPot), k)
    
    for(i in 1:dim(counttot)[1])
    {
      if(counttot[i,"idPot"] %in% trainsample) {
        counttot[i,"set"] <- "train"
      } else{
        counttot[i,"set"] <- "test"
      }
    }
    
    testsample <- sample(levels(as.factor(as.character(counttot[counttot$set=="test" & counttot$rosette_DM>0, "idPot"]))), 100)
    
    test <- droplevels(counttot[counttot$idPot %in% testsample, ])
    train <- droplevels(counttot[counttot$idPot %in% trainsample, ])
    train <- na.omit(train)
    test <- na.omit(test)
    test <- droplevels(test)
    train <- droplevels(train)
    
    # lm
    model0 <- lm(rosette_DM ~ Ros_area + Ros_perim + Ros_circ + Ros_AR + Ros_round, data=train)
    test$predDM0 <- predict(model0, newdata=test)
    lm_accuracy[n] <- as.numeric(((cor.test(test$rosette_DM, test$predDM0)$estimate)^2)) 
    
    # quadratic
    model1 <- lm(rosette_DM ~ poly(Ros_area,2)+poly(Ros_perim,2)+poly(Ros_circ,2)+poly(Ros_AR,2)+poly(Ros_round,2), data=train)
    test$predDM1 <- predict(model1, newdata=test)
    quad_accuracy[n] <- as.numeric(((cor.test(test$rosette_DM, test$predDM1)$estimate)^2)) 
    
    # RIDGE model
    model4 <- glmnet(x=as.matrix(train[,c("Ros_area", "Ros_perim", "Ros_circ", "Ros_AR", "Ros_round")]), 
                     y=as.matrix(train[,"rosette_DM"]),
                     family="gaussian",
                     alpha = 0) # for lasso, alpha=1, for ridge alpha=0)
    cv.out2 <- cv.glmnet(x=as.matrix(train[,c("Ros_area", "Ros_perim", "Ros_circ", "Ros_AR", "Ros_round")]), 
                         y=as.matrix(train[,"rosette_DM"]), alpha = 0)
    bestlam2 <- cv.out2$lambda.min
    test$predDM4 <- predict(model4,
                            newx=as.matrix(test[,c("Ros_area", "Ros_perim", "Ros_circ", "Ros_AR", "Ros_round")]), 
                            s = bestlam2)
    RIDGE_accuracy[n] <- as.numeric(((cor.test(test$rosette_DM, test$predDM4)$estimate)^2)) 
    
    if(n %in% seq(10,100,10)){print(n)}
  }
  nit <- c(1:100)
  ntrain <- rep(k, 100)
  tot_accuracy <- cbind(nit, ntrain)
  tot_accuracy <- cbind(tot_accuracy, lm_accuracy)
  tot_accuracy <- cbind(tot_accuracy, quad_accuracy)
  tot_accuracy <- cbind(tot_accuracy, LASSO_accuracy)
  tot_accuracy <- cbind(tot_accuracy, RIDGE_accuracy)
  
  temp0 <- rbind(temp0, tot_accuracy)
  print(k)
  
  write.table(temp0, "Results_models_DMprediction.csv", dec=".", sep=",", row.names = F)
  
}

