# ADRIFT Fin whale click processing
# Updated 2023-09-21: First draft

# Change to wherever all the pre-val databases and binaries are for you
dbFolder <- './Data/Pre_Validation/'
binFolder <- './Data/ADRIFT_FIN/Binaries'
library(PAMpal)
library(here)
# path to finClick.R file I sent you
source(here('finClick.R'))

pps <- PAMpalSettings(db=dbFolder, binaries=binFolder,
                      sr_hz='auto', #change this if it was decimated
                      filterfrom_khz=0,
                      filterto_khz=NULL,
                      winLen_sec=1,
                      threshold=-5,
                      functions=list('ClickDetector'=finClick))

data <- processPgDetections(pps, mode='db', id='All_ADRIFT_Fin_Preval')
# please send Taiki this file!
saveRDS(data, file='All_ADRIFT_Fin_Preval.rds')

# Next step we'll get the validated UIDs from the set of validated DBs
valFolder <- '../Data/ADRIFT_FIN/Validated/'
valFiles <- list.files(valFolder, full.names=TRUE)
library(dplyr)
valData <- lapply(valFiles, function(x) {
    # this just pulls the event/click info from the DB
    result <- PAMpal:::getDbData(x)
    if(is.null(result) || nrow(result) == 0) {
        return(result)
    }
    result$db <- basename(x)
    result
}) %>% bind_rows()
data <- readRDS('../Data/ADRIFT_FIN/All_ADRIFT_Fin_Preval.rds')
clickData <- getClickData(data)
# clicks are fin whales if their IDs are in the validated set
clickData <- bind_rows(lapply(split(clickData, clickData$db), function(x) {
    thisDbVal <- valData[valData$db == basename(x$db[1]), ]
    if(is.null(thisDbVal) || nrow(thisDbVal) == 0) {
        x$isFin <- FALSE
        return(x)
    }
    x$isFin <- x$UID %in% thisDbVal$UID
    x
}))

# clickData$isFin <- clickData$UID %in% valData$UID
# should be like 45k not, 1400 fin
table(clickData$isFin)
library(ggplot2)
gslope <- ggplot(clickData) +
    geom_density(aes(x=freqSlope, col=isFin)) +
    xlim(-.05, .05)
table(clickData$isFin, clickData$detectorName)
library(patchwork)
gpeak <- ggplot(clickData) +
    geom_density(aes(x=peak, col=isFin)) +
    xlim(0, .05)
gslope+gpeak
ggplot(clickData) +
    geom_density(aes(x=fmax_10dB, col=isFin))

library(caret)
library(randomForest)

hm <- slice_sample(clickData, by=isFin, n=1400)
hm <- PAMpal:::dropCols(hm, c('UID', 'UTC', 'Channel', 'angle',
                              'angleError', 'BinaryFile', 'eventLabel',
                              'eventId', 'db', 'species', 'predicted', 'peakTime'))
hm$isFin <- factor(hm$isFin)
library(randomForest)

mtryTest <- tuneRF(x=hm[-26], y=hm[[26]], stepFactor = 1.5, improve=1e-5, ntreeTry=500)
mtryBest <- mtryTest[which.min(mtryTest[, 'OOBError']), 'mtry']
rf <- randomForest(isFin ~ .,
                   data=hm, mtry=10, ntree=1.5e3)

summary(rf)
rfp <- rfPermute(isFin ~ ., data=hm, proximity=TRUE)
plotImportance(rfp)
summary(rf)

preds <- predict(rfp, newdata=PAMpal:::dropCols(clickData, c('UID', 'UTC', 'Channel', 'angle',
                                                            'angleError', 'BinaryFile', 'eventLabel',
                                                            'eventId', 'db', 'species', 'isFin')))
table(preds, clickData$isFin)
varImpPlot(rf)
