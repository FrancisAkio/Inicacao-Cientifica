---
title: "NBLDA, PLDA and PLDA2 example for RNA-Seq"
subtitle: "Two classes"
author: "Francis Akio Ishizaka"
date: "16/09/2020"
output: 
  pdf_document:
    keep_md: yes
---

This implementation classifies non-tumor and tumor.

## Setup

```r
library(MLSeq)
library(tidyverse)
library(S4Vectors)
library(DESeq2)

db <- read_tsv("GSE146889.tsv")
```

### Remove obvious non-expressive genes


```r
db <- db[rowSums(db[,-1]) > 1,]
```


## Define classes non-tumor (N) vs. tumor (T)

```r
classe <- colnames(db) %>% str_detect("normal") %>% ifelse("N", "T") %>% factor()
classe <- classe[-1]
classe <- DataFrame(cnc = classe)
```

## For debuging, select 100 genes

```r
debug <- F
if(debug){
  vars <- apply(db[,-1], 1, var, na.rm = T)
  names(vars) <- db$GeneId
  vars <- vars %>% sort(decreasing = T)
  data <- db[db$GeneId %in% names(vars)[1:100], ]
} else data<- db
```

## Partition data for cross-validation
For reproducibility reasons, it was select 171752 as seed.

 70% for training and 30% for test.


```r
genes <- data[,1]
data <- data[,-1]
set.seed(171752)
ind <- sample(ncol(data), ceiling(ncol(data)* 0.7), F)
dTr <- as.matrix(data[, ind] + 1)
dTs <- as.matrix(data[,-ind] + 1)
cTr <- DataFrame(cnc = classe[ ind,])
cTs <- DataFrame(cnc = classe[-ind,])
data.train <- DESeqDataSetFromMatrix(countData = dTr, colData = cTr,
                                     design = formula(~cnc))
```

```
## converting counts to integer mode
```

```r
data.test <- DESeqDataSetFromMatrix(countData = dTs, colData = cTs,
                                     design = formula(~cnc))
```

```
## converting counts to integer mode
```

## PLDA Classification

```r
mod <- classify(data = data.train, method = "PLDA",
                normalize = "deseq", ref = "N",
                control = discreteControl(method = "repeatedcv", number = 5,
                                          repeats = 10))
trained(mod)
```

```
## 
## Poisson Linear Discriminant Analysis (PLDA)
## 
##    124 samples 
##  48608 predictors 
##      2 classes: 'T', 'N'  (Reference category: 'N') 
## 
## Normalization: DESeq median ratio.
## Power transformation is  NOT performed.
## Resampling: Cross-Validated (5 fold, repeated 10 times)
## Summary of sample sizes: 100, 99, 98, 100, 99, ...
## Summary of selected features:  All features are selected. 
## 
##        rho   Avg.Error   Avg.NonZeroFeat.   Accuracy 
##    0.00000        2.02          48608.00      0.9184 
##  103.39514        2.44           2152.70      0.9015 
##  206.79029        2.46            610.98      0.9008 
##  310.18543        2.48            269.86      0.8999 
##  413.58057        2.74            159.10      0.8894 
##  516.97572        2.88             96.06      0.8838 
##  620.37086        3.02             69.70      0.8781 
##  723.76600        3.18             50.60      0.8716 
##  827.16115        3.32             32.12      0.8661 
##  930.55629        3.30             21.52      0.8669 
## 1033.95143        3.40             14.90      0.8629 
## 1137.34658        3.46             11.30      0.8605 
## 1240.74172        3.50              8.26      0.8588 
## 1344.13686        3.56              6.74      0.8564 
## 1447.53201        3.58              5.86      0.8556 
## 1550.92715        3.58              4.96      0.8556 
## 1654.32229        3.56              4.16      0.8564 
## 1757.71744        3.56              3.58      0.8564 
## 1861.11258        3.56              2.90      0.8564 
## 1964.50772        3.56              2.32      0.8564 
## 2067.90287        3.60              1.48      0.8548 
## 2171.29801        3.64              1.04      0.8531 
## 2274.69315        3.66              1.00      0.8523 
## 2378.08830        3.66              1.00      0.8523 
## 2481.48344        3.66              1.00      0.8523 
## 2584.87858        3.66              1.00      0.8523 
## 2688.27373        3.68              1.00      0.8515 
## 2791.66887        3.74              1.00      0.8491 
## 2895.06402        3.74              1.00      0.8491 
## 2998.45916        3.74              1.00      0.8491 
## 
## The optimum model is obtained when rho = 0.00000 with an overall accuracy of
## Accuracy = 0.9184 over folds. On the average 48608 out of 48608 features was used
## in the classifier.
## 
## NOTES: The optimum model is selected using tuning parameter 'rho' which achieves 
##  the lowest classification error (Avg.Error) or the highest Accuracy. The  
##  classification error is given as the average number of misclassified samples 
##  over cross-validated folds. Similarly, the 'Avg.NonZeroFeat.' is the average  
##  number of non-zero features (the selected variables in the classification task)   
##  over cross-validated folds. As the number of non-zero features decreases, the  
##  model becomes more sparse.
```

### Prediction

```r
pred.mod <- predict(mod, data.test)

pred.mod <- relevel(pred.mod, ref = "N")
verd     <- relevel(cTs$cnc, ref = "N")

tbl <- table(Predito = pred.mod, Verd = verd)
confusionMatrix(tbl, positive = "N")
```

```
## Confusion Matrix and Statistics
## 
##        Verd
## Predito  N  T
##       N 21  3
##       T  7 21
##                                           
##                Accuracy : 0.8077          
##                  95% CI : (0.6747, 0.9037)
##     No Information Rate : 0.5385          
##     P-Value [Acc > NIR] : 4.796e-05       
##                                           
##                   Kappa : 0.6176          
##                                           
##  Mcnemar's Test P-Value : 0.3428          
##                                           
##             Sensitivity : 0.7500          
##             Specificity : 0.8750          
##          Pos Pred Value : 0.8750          
##          Neg Pred Value : 0.7500          
##              Prevalence : 0.5385          
##          Detection Rate : 0.4038          
##    Detection Prevalence : 0.4615          
##       Balanced Accuracy : 0.8125          
##                                           
##        'Positive' Class : N               
## 
```

```r
resultadoPLDA <- confusionMatrix(tbl, positive = "N")
resultadoPLDA$overall['Accuracy']
```

```
##  Accuracy 
## 0.8076923
```

## PLDA2 Classification

```r
mod1 <- classify(data = data.train, method = "PLDA2",
                normalize = "deseq", ref = "N",
                control = discreteControl(method = "repeatedcv", number = 5,
                                          repeats = 10))
trained(mod1)
```

```
## 
## Poisson Linear Discriminant Analysis with Power Transformation (PLDA2)
## 
##    124 samples 
##  48608 predictors 
##      2 classes: 'T', 'N'  (Reference category: 'N') 
## 
## Normalization: DESeq median ratio.
## Power transformation is performed.
## Resampling: Cross-Validated (5 fold, repeated 10 times)
## Summary of sample sizes: 100, 99, 99, 99, 99, ...
## Summary of selected features:  All features are selected. 
## 
##        rho   Avg.Error   Avg.NonZeroFeat.   Accuracy 
##    0.00000        2.50          48608.00      0.8990 
##    1.47312        2.46          12976.82      0.9006 
##    2.94624        2.46           6456.62      0.9006 
##    4.41936        2.48           3645.08      0.8998 
##    5.89248        2.52           2182.34      0.8982 
##    7.36560        2.48           1362.24      0.8998 
##    8.83872        2.48            897.32      0.8998 
##   10.31184        2.56            580.40      0.8965 
##   11.78496        2.64            371.00      0.8933 
##   13.25808        2.52            232.26      0.8981 
##   14.73119        2.48            149.96      0.8998 
##   16.20431        2.42            103.56      0.9023 
##   17.67743        2.48             74.56      0.8998 
##   19.15055        2.44             55.00      0.9014 
##   20.62367        2.60             41.52      0.8950 
##   22.09679        2.72             32.84      0.8902 
##   23.56991        2.82             24.94      0.8861 
##   25.04303        2.82             17.10      0.8861 
##   26.51615        2.86             12.36      0.8846 
##   27.98927        2.82              9.52      0.8862 
##   29.46239        2.82              6.86      0.8862 
##   30.93551        2.86              3.90      0.8846 
##   32.40863        2.88              2.30      0.8838 
##   33.88175        2.88              2.04      0.8838 
##   35.35487        2.88              2.02      0.8838 
##   36.82799        2.92              2.00      0.8822 
##   38.30111        2.94              1.96      0.8814 
##   39.77423        2.94              1.64      0.8814 
##   41.24734        2.94              1.22      0.8814 
##   42.72046        2.92              1.00      0.8822 
## 
## The optimum model is obtained when rho = 16.20431 with an overall accuracy of
## Accuracy = 0.9023 over folds. On the average 103.56 out of 48608 features was used
## in the classifier.
## 
## NOTES: The optimum model is selected using tuning parameter 'rho' which achieves 
##  the lowest classification error (Avg.Error) or the highest Accuracy. The  
##  classification error is given as the average number of misclassified samples 
##  over cross-validated folds. Similarly, the 'Avg.NonZeroFeat.' is the average  
##  number of non-zero features (the selected variables in the classification task)   
##  over cross-validated folds. As the number of non-zero features decreases, the  
##  model becomes more sparse.
```

### Prediction

```r
pred.mod1 <- predict(mod1, data.test)

pred.mod1 <- relevel(pred.mod1, ref = "N")
verd     <- relevel(cTs$cnc, ref = "N")

tbl <- table(Predito = pred.mod1, Verd = verd)
confusionMatrix(tbl, positive = "N")
```

```
## Confusion Matrix and Statistics
## 
##        Verd
## Predito  N  T
##       N 21  1
##       T  7 23
##                                           
##                Accuracy : 0.8462          
##                  95% CI : (0.7192, 0.9312)
##     No Information Rate : 0.5385          
##     P-Value [Acc > NIR] : 2.871e-06       
##                                           
##                   Kappa : 0.6959          
##                                           
##  Mcnemar's Test P-Value : 0.0771          
##                                           
##             Sensitivity : 0.7500          
##             Specificity : 0.9583          
##          Pos Pred Value : 0.9545          
##          Neg Pred Value : 0.7667          
##              Prevalence : 0.5385          
##          Detection Rate : 0.4038          
##    Detection Prevalence : 0.4231          
##       Balanced Accuracy : 0.8542          
##                                           
##        'Positive' Class : N               
## 
```

```r
resultadoPLDA2 <- confusionMatrix(tbl, positive = "N")
resultadoPLDA2$overall['Accuracy']
```

```
##  Accuracy 
## 0.8461538
```

## NBLDA Classification

```r
mod2 <- classify(data = data.train, method = "NBLDA",
                normalize = "deseq", ref = "N", 
                control = discreteControl(method = "repeatedcv", number = 5,
                                          repeats = 10, classProbs = T))
```

### Prediction

```r
pred.mod2 <- predict(mod2, data.test)

pred.mod2 <- relevel(pred.mod2, ref = "N")
verd     <- relevel(cTs$cnc, ref = "N")

tbl <- table(Predito = pred.mod2, Verd = verd)
confusionMatrix(tbl, positive = "N")
```

```
## Confusion Matrix and Statistics
## 
##        Verd
## Predito  N  T
##       N 23  1
##       T  5 23
##                                           
##                Accuracy : 0.8846          
##                  95% CI : (0.7656, 0.9565)
##     No Information Rate : 0.5385          
##     P-Value [Acc > NIR] : 9.883e-08       
##                                           
##                   Kappa : 0.7706          
##                                           
##  Mcnemar's Test P-Value : 0.2207          
##                                           
##             Sensitivity : 0.8214          
##             Specificity : 0.9583          
##          Pos Pred Value : 0.9583          
##          Neg Pred Value : 0.8214          
##              Prevalence : 0.5385          
##          Detection Rate : 0.4423          
##    Detection Prevalence : 0.4615          
##       Balanced Accuracy : 0.8899          
##                                           
##        'Positive' Class : N               
## 
```

```r
resultadoNBLDA <- confusionMatrix(tbl, positive = "N")
resultadoNBLDA$overall['Accuracy']
```

```
##  Accuracy 
## 0.8846154
```
