### 04.10
### chane col's order
### add a simplify table  (table 1,  name == summary)
### rename 'NA' into '-' in pleio/heter cols 

# table 1:  simplified summary table
# table 2: total summary table
# table 3: MR result
# table 4: heterogeneity result
# table 5: pleiotropy result
# table 6: directionality result
 
MR_summary <- function(mr, heter, pleio, steiger){
 
  library(data.table)
#  library(openxlsx)
  library(dplyr)
  library(stringr)
  
  tmpdata <- as.data.frame(mr)
  tmpheter <- as.data.frame(heter)
  tmppleio <- as.data.frame(pleio)

#### 1. remove additional header ...
  tmpdata <- subset(tmpdata, method != 'method')
  tmpheter <- subset(tmpheter, method == 'Inverse variance weighted')
  tmppleio <- subset(tmppleio, exposure != 'exposure' )


#### 2. remove colnames == c('id.exposure', 'id.outcome')
  tmpdata <- tmpdata[, !colnames(tmpdata) %in% c('id.exposure', 'id.outcome')]
  tmpheter <- tmpheter[, !colnames(tmpheter) %in% c('id.exposure', 'id.outcome')]
  tmppleio <- tmppleio[, !colnames(tmppleio) %in% c('id.exposure', 'id.outcome')]

#### 3. order by colnames
  tmpdata <- tmpdata %>% select(exposure, outcome, nsnp, pval, method, ##select函数的作用是对列进行排序，修改列名等操作
                                or, or_lci95, or_uci95, everything()) #everything即剩下的所有列
  tmpheter <- tmpheter %>% select(exposure, outcome, everything())
  tmppleio <- tmppleio %>% select(exposure, outcome, everything())

#### 4. add barcode
  tmpdata1 <- data.frame(tmpdata, barcode = paste(tmpdata$exposure, tmpdata$outcome, sep = '-'))
  tmpheter1 <- data.frame(heterogeneity = tmpheter$Q_pval, barcode = paste(tmpheter$exposure, tmpheter$outcome, sep = '-'))
  tmppleio1 <- data.frame(pleiotropy = tmppleio$pval, barcode = paste(tmppleio$exposure, tmppleio$outcome, sep = '-'))

  tmpheter1[which(as.numeric(tmpheter1$heterogeneity) < 0.05), 'heterogeneity'] <- 'Yes'
  tmpheter1[which(as.numeric(tmpheter1$heterogeneity) >= 0.05), 'heterogeneity'] <- 'No'
  
  tmppleio1[which(as.numeric(tmppleio1$pleiotropy) < 0.05), 'pleiotropy'] <- 'Yes'
  tmppleio1[which(as.numeric(tmppleio1$pleiotropy) >= 0.05), 'pleiotropy'] <- 'No'
  
  
#### optional: deal with steiger result
  if(!is.null(steiger)){
    
    tmpsteiger <- as.data.frame(steiger)
    
    tmpsteiger <- subset(tmpsteiger, exposure != 'exposure')

    tmpsteiger <- tmpsteiger[, ! colnames(tmpsteiger) %in% c('id.exposure', 'id.outcome')]
    tmpsteiger <- arrange(tmpsteiger, exposure, outcome)  ##arrange函数的作用是排序，优先按exposure排序，然后按outcome排序
    tmpsteiger1 <- data.frame(directionality = tmpsteiger$correct_causal_direction, 
                              barcode = paste(tmpsteiger$exposure, tmpsteiger$outcome, sep = '-'))
    tmp_other <- Reduce(function(x, y) merge(x, y, by = 'barcode'), list(tmpheter1, tmppleio1, tmpsteiger1))
  } else { 
      tmp_other <- merge(tmpheter1, tmppleio1, by = 'barcode')
      }   #reduce函数的作用是合并多个数据框，merge的作用是取交集合并
    
  
#### 5. summary table
  tmp_summary <- merge(tmpdata1, tmp_other, by = 'barcode', all.x = T) #all.x=T表示左连接，即保留x的所有，再加上能与y连接的部分
  
  tmp_summary[which(is.na(tmp_summary$pleiotropy) == T), 'pleiotropy'] <- '-'
  tmp_summary[which(is.na(tmp_summary$heterogeneity) == T), 'heterogeneity'] <- '-'
  
  
  if(is.null(steiger)){
    tmp_summary$directionality <- 'Not tested'
  }

  
### heter, pleio, stegier(if exist)  
  tmpMRres <- tmp_summary[, !colnames(tmp_summary) %in% c('heterogeneity', 'pleiotropy', 'directionality')] %>% arrange(exposure, outcome)
  tmpheter_final <- tmpheter %>% arrange(exposure, outcome)
  tmppleio_final <- tmppleio %>% arrange(exposure, outcome)

  if (is.null(steiger)) {
    tmpsteiger_final <- data.frame(`No Steiger result found` = '')
  } else  {
    tmpsteiger_final <- tmpsteiger
  }

  tmptable <- tmp_summary %>% select(exposure, outcome, nsnp, b, se, pval, method, heterogeneity, pleiotropy, 
                                     directionality, or, or_lci95, or_uci95)
  
  colnames(tmptable) <- c('exposure', 'outcome', 'No. of SNPs in MRA', 'Beta', 'Se', 'MR pvalue', 'Method', 'heterogeneity', 'pleiotropy', 
                          'directionality', 'OR', 'LCI95', 'UCI95')

## IV < 2                      
  tmptable[which(tmptable$Method == 'Wald ratio'), 'select'] <- 1
## IV == 2 
  tmptable[which(tmptable$heterogeneity == 'No' & tmptable$pleiotropy == '-' & 
                   tmptable$Method == 'Inverse variance weighted (fixed effects)'), 'select'] <- 1
  tmptable[which(tmptable$heterogeneity == 'Yes' & tmptable$pleiotropy == '-' & 
                   tmptable$Method == 'Inverse variance weighted'), 'select'] <- 1
  
## IV > 2
  tmptable[which(tmptable$heterogeneity == 'No' & tmptable$pleiotropy == 'No' & 
                   tmptable$Method == 'Inverse variance weighted (fixed effects)'), 'select'] <- 1
  
  tmptable[which(tmptable$heterogeneity == 'Yes' & tmptable$pleiotropy == 'No' & 
                   tmptable$Method == 'Inverse variance weighted'), 'select'] <- 1
  
if(unique(tmptable$pleiotropy) == 'Yes')
{
  tmptable[which(tmptable$Method == 'MR Egger'), 'select'] <- 1
  tmptable[which(tmptable$Method == 'Weighted median'), 'select'] <- 1
}

# ## select_all == 1 + 2; select_sig == 1  
#   tmptable[which(tmptable$heterogeneity == 'Yes' & tmptable$pleiotropy == 'No' &
#                    tmptable$Method == 'MR Egger'), 'select'] <- 2
#   tmptable[which(tmptable$heterogeneity == 'Yes' & tmptable$pleiotropy == 'No' &
#                    tmptable$Method == 'Weighted median'), 'select'] <- 2
#   
#   tmptable[which(tmptable$heterogeneity == 'Yes' & tmptable$pleiotropy == 'Yes' & 
#                    tmptable$Method == 'Weighted median'), 'select'] <- 2
#   
#   tmptable[which(tmptable$heterogeneity == 'No' & tmptable$pleiotropy == 'Yes' & 
#                    tmptable$Method == 'MR Egger'), 'select'] <- 2
#   
  
  tmptable_final <- subset(tmptable, select == 1)
  tmptable_final <- tmptable_final[, which(colnames(tmptable_final) != 'select')]
  tmptable_final <- arrange(tmptable_final, exposure, outcome)
  
  # tmptable_total <- subset(tmptable, select %in% c(1, 2))
  # tmptable_total <- tmptable_total[, which(colnames(tmptable_total) != 'select')]
  # tmptable_total <- arrange(tmptable_total, exposure, outcome)
  
  
  # result <- list(tmptable_final, tmptable_total, tmpMRres, tmpheter_final, tmppleio_final, tmpsteiger_final)
  # names(result) <- c('summary', 'summaty_total', 'MRresult', 'heterogeneity', 'pleiotropy', 'directionality')

  result <- list(tmptable_final, tmpMRres, tmpheter_final, tmppleio_final, tmpsteiger_final)
  names(result) <- c('summary', 'MRresult', 'heterogeneity', 'pleiotropy', 'directionality')
  
  return(result)
}


  