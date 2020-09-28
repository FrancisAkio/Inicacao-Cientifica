---
title: "NBLDA, PLDA and PLDA2 example for RNA-Seq"
subtitle: "Four classes"
author: "Francis Akio Ishizaka"
date: "16/09/2020"
output: 
  pdf_document:
    keep_md: yes
---

This implementation classifies non-tumor (N), MSS, hypermethylation related MSI (MSI_H) and Lynch Syndrome related MSI (MSI_LS).

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


## Define classes non-tumor (N) vs. MSS vs. hypermethylation related MSI (MSI_H) vs. Lynch Syndrome related MSI (MSI_LS)

```r
classe <- colnames(db) %>%
  str_detect("normal") %>%
  ifelse("N", str_detect(colnames(db),"MSS") %>% 
           ifelse("MSS", str_detect(colnames(db),"MSI_MLH1HM") %>% ifelse("MSI_H", "MSI_LS"))) %>%
   factor()
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

Choosed 70% for training and 30% for test.


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
##      4 classes: 'MSI_LS', 'N', 'MSS', 'MSI_H'  (Reference category: 'N') 
## 
## Normalization: DESeq median ratio.
## Power transformation is  NOT performed.
## Resampling: Cross-Validated (5 fold, repeated 10 times)
## Summary of sample sizes: 98, 99, 100, 99, 100, ...
## Summary of selected features:  All features are selected. 
## 
##        rho   Avg.Error   Avg.NonZeroFeat.   Accuracy 
##    0.00000        9.14          48608.00      0.6315 
##  103.39514       10.68           2507.58      0.5697 
##  206.79029       11.16            712.00      0.5501 
##  310.18543       11.48            317.50      0.5374 
##  413.58057       11.40            180.44      0.5406 
##  516.97572       11.80            111.42      0.5244 
##  620.37086       12.12             77.88      0.5118 
##  723.76600       12.30             56.96      0.5047 
##  827.16115       12.32             37.56      0.5038 
##  930.55629       11.54             25.24      0.5350 
## 1033.95143       11.20             17.76      0.5483 
## 1137.34658       10.98             13.86      0.5571 
## 1240.74172       10.72             10.82      0.5679 
## 1344.13686       10.88              9.38      0.5614 
## 1447.53201       10.64              8.30      0.5711 
## 1550.92715       10.28              6.78      0.5856 
## 1654.32229       10.12              5.58      0.5923 
## 1757.71744        9.96              4.68      0.5987 
## 1861.11258        9.88              3.80      0.6020 
## 1964.50772        9.82              2.90      0.6043 
## 2067.90287        9.80              2.02      0.6051 
## 2171.29801        9.66              1.46      0.6107 
## 2274.69315        9.62              1.14      0.6123 
## 2378.08830        9.62              1.06      0.6123 
## 2481.48344        9.54              1.02      0.6156 
## 2584.87858        9.64              1.00      0.6114 
## 2688.27373        9.80              1.00      0.6050 
## 2791.66887       10.48              1.00      0.5778 
## 2895.06402       12.40              1.00      0.5006 
## 2998.45916       14.22              1.00      0.4267 
## 
## The optimum model is obtained when rho = 0.00000 with an overall accuracy of
## Accuracy = 0.6315 over folds. On the average 48608 out of 48608 features was used
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
##         Verd
## Predito   N MSI_H MSI_LS MSS
##   N      21     1      0   1
##   MSI_H   2     4      6   0
##   MSI_LS  5     1      6   0
##   MSS     0     0      2   3
## 
## Overall Statistics
##                                           
##                Accuracy : 0.6538          
##                  95% CI : (0.5091, 0.7803)
##     No Information Rate : 0.5385          
##     P-Value [Acc > NIR] : 0.06199         
##                                           
##                   Kappa : 0.48            
##                                           
##  Mcnemar's Test P-Value : NA              
## 
## Statistics by Class:
## 
##                      Class: N Class: MSI_H Class: MSI_LS Class: MSS
## Sensitivity            0.7500      0.66667        0.4286    0.75000
## Specificity            0.9167      0.82609        0.8421    0.95833
## Pos Pred Value         0.9130      0.33333        0.5000    0.60000
## Neg Pred Value         0.7586      0.95000        0.8000    0.97872
## Prevalence             0.5385      0.11538        0.2692    0.07692
## Detection Rate         0.4038      0.07692        0.1154    0.05769
## Detection Prevalence   0.4423      0.23077        0.2308    0.09615
## Balanced Accuracy      0.8333      0.74638        0.6353    0.85417
```

```r
resultadoPLDA <- confusionMatrix(tbl, positive = "N")
resultadoPLDA$overall['Accuracy']
```

```
##  Accuracy 
## 0.6538462
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
##      4 classes: 'MSI_LS', 'N', 'MSS', 'MSI_H'  (Reference category: 'N') 
## 
## Normalization: DESeq median ratio.
## Power transformation is performed.
## Resampling: Cross-Validated (5 fold, repeated 10 times)
## Summary of sample sizes: 99, 100, 100, 99, 98, ...
## Summary of selected features:  All features are selected. 
## 
##        rho   Avg.Error   Avg.NonZeroFeat.   Accuracy 
##    0.00000        7.88          48608.00      0.6825 
##    1.47312        8.88          18092.40      0.6422 
##    2.94624        9.58           8022.78      0.6139 
##    4.41936        9.24           4168.18      0.6274 
##    5.89248        8.76           2386.72      0.6469 
##    7.36560        8.80           1463.66      0.6454 
##    8.83872        8.82            959.08      0.6445 
##   10.31184        8.80            616.56      0.6453 
##   11.78496        8.86            394.22      0.6428 
##   13.25808        8.86            248.06      0.6427 
##   14.73119        8.90            161.64      0.6413 
##   16.20431        9.18            111.54      0.6299 
##   17.67743        9.42             79.92      0.6204 
##   19.15055        9.58             57.94      0.6141 
##   20.62367        9.56             43.64      0.6148 
##   22.09679        9.62             34.80      0.6124 
##   23.56991        9.40             26.90      0.6212 
##   25.04303        9.40             18.24      0.6215 
##   26.51615        9.26             13.40      0.6273 
##   27.98927        9.28             10.12      0.6264 
##   29.46239        9.28              7.12      0.6263 
##   30.93551        9.18              4.30      0.6304 
##   32.40863        9.02              2.64      0.6367 
##   33.88175        8.84              2.24      0.6440 
##   35.35487        8.72              2.06      0.6489 
##   36.82799        8.52              2.04      0.6568 
##   38.30111        8.46              1.92      0.6592 
##   39.77423        8.44              1.58      0.6600 
##   41.24734        8.40              1.24      0.6616 
##   42.72046        8.68              1.06      0.6502 
## 
## The optimum model is obtained when rho = 0.00000 with an overall accuracy of
## Accuracy = 0.6825 over folds. On the average 48608 out of 48608 features was used
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
##         Verd
## Predito   N MSI_H MSI_LS MSS
##   N      21     1      0   0
##   MSI_H   1     0      2   0
##   MSI_LS  6     4     10   0
##   MSS     0     1      2   4
## 
## Overall Statistics
##                                           
##                Accuracy : 0.6731          
##                  95% CI : (0.5289, 0.7967)
##     No Information Rate : 0.5385          
##     P-Value [Acc > NIR] : 0.03412         
##                                           
##                   Kappa : 0.4983          
##                                           
##  Mcnemar's Test P-Value : NA              
## 
## Statistics by Class:
## 
##                      Class: N Class: MSI_H Class: MSI_LS Class: MSS
## Sensitivity            0.7500      0.00000        0.7143    1.00000
## Specificity            0.9583      0.93478        0.7368    0.93750
## Pos Pred Value         0.9545      0.00000        0.5000    0.57143
## Neg Pred Value         0.7667      0.87755        0.8750    1.00000
## Prevalence             0.5385      0.11538        0.2692    0.07692
## Detection Rate         0.4038      0.00000        0.1923    0.07692
## Detection Prevalence   0.4231      0.05769        0.3846    0.13462
## Balanced Accuracy      0.8542      0.46739        0.7256    0.96875
```

```r
resultadoPLDA2 <- confusionMatrix(tbl, positive = "N")
resultadoPLDA2$overall['Accuracy']
```

```
##  Accuracy 
## 0.6730769
```
