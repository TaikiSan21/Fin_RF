db <- 'Data/Pre_Validation/'
bin <- 'Data/ADRIFT/Binaries/'
library(PAMpal)
dbFiles <- list.files(db, full.names=TRUE)
binDirs <- list.dirs(bin, full.names=TRUE, recursive=FALSE)
studyList <- vector('list', length=length(dbFiles))
for(i in seq_along(dbFiles)) {
    cat(i)
    pps <- PAMpalSettings(db=dbFiles[i], binaries = binDirs[i], sr_hz=200,
                          filterfrom_khz=0, filterto_khz=NULL,
                          winLen_sec=1,
                          threshold=-5,
                          functions=list('ClickDetector'=finClick))
    binList <- pps@binaries$list
    whichCopy <- grepl('\\(1\\)', binList)
    binList <- binList[!whichCopy]
    pps@binaries$list <- binList
    allData <- processPgDetections(pps, mode='recording', id=basename(binDirs[i]))
    studyList[[i]] <- allData
}
saveRDS(studyList, file='SeparateStudyListAll.rds')
allData <- bindStudies(readRDS('SeparateStudyListAll.rds'))
allPps <- PAMpalSettings(db=db, binaries = bin, sr_hz=200,
                         filterfrom_khz=0, filterto_khz=NULL,
                         winLen_sec=1,
                         threshold=-5,
                         functions=list('ClickDetector'=finClick))
# appear to be lots of copies of binaries with (1) in them
sum(grepl('\\(1\\)', allPps@binaries$list))
binList <- allPps@binaries$list
whichCopy <- grepl('\\(1\\)', binList)
binList <- binList[!whichCopy]
# these are all "1" indicates they are all copies
sapply(allPps@binaries$list[whichCopy], function(x) {
    x <- gsub('\\(1\\)', '', x)
    sum(grepl(x, binList))
})
allPps@binaries$list <- binList
allData <- processPgDetections(allPps, mode='recording', id='ADRIFT_FinAll')
# very weird duplicate click UID 2240000011 in event '5296.210611015734.wav'
# is in two binary files, second one is clearly mistake since it starts at 11
# then switches to proper bins.
# Same problem 2835000008 / "5296.210613095813.wav"
allData[['5296.210611015734.wav']]$Click_Detector_1 <- allData[['5296.210611015734.wav']]$Click_Detector_1[-9,]
allData[["5296.210613095813.wav"]]$Click_Detector_1 <- allData[["5296.210613095813.wav"]]$Click_Detector_1[-8, ]
allData[["5296.210616005830.wav"]]$Click_Detector_1 <- allData[["5296.210616005830.wav"]]$Click_Detector_1[-12, ]
allData <- calculateICI(allData, time='peakTime')
saveRDS(allData, file='Data/AllADRIFTDupRemove_ici_Study.rds')
allData <- readRDS('Data/AllADRIFTDupRemove_ici_Study.rds')
allIci <- getICI(allData, type='data')
allClick <- getClickData(allData)
allClick <- joinIci(allClick, allIci)


allPdf <- makePredDf(newModel, allClick)
# can mark our pre-val wav files from the database event names
library(RSQLite)
preDb <- list.files('Data/Pre_Validation/', full.names=TRUE)
preEv <- lapply(preDb, function(x) {
    con <- dbConnect(x, drv=SQLite())
    on.exit(dbDisconnect(con))
    ev <- dbReadTable(con, 'Click_Detector_OfflineEvents')
    wavs <- gsub('Added by PAMpal, event ID: ', '', ev$comment)
    wavs <- gsub(' ', '', wavs)
    wavs
})
names(preEv) <- basename(preDb)
allPdf <- bind_rows(lapply(split(allPdf, allPdf$db), function(x) {
    thisWavs <- preEv[[basename(x$db[1])]]
    x$validated <- x$eventId %in% thisWavs
    x
}))
allPdf$drift <- gsub('(ADRIFT_[0-9]{1,3})_.*', '\\1', basename(allPdf$db))
allPdf$hour <- floor_date(allPdf$UTC, unit='1hour')
allPdf$sixMin <- floor_date(allPdf$UTC, unit='6min')
allPdf <- allPdf %>%
    group_by(eventId) %>%
    mutate(sixMin=as.numeric(difftime(UTC, min(UTC), units='mins')) / 6,
           sixMin=floor(sixMin),
           sixMin=paste0(eventId, '_', sixMin)) %>%
    ungroup()

allSix <- summarisePredEvent(filter(allPdf, !validated), thresh=.5, evCol=c('drift', 'thirty','hour', 'sixMin'))

ggplot() + geom_density(data=allSix, aes(x=meanPos)) +
    geom_density(data=sixTest, aes(x=meanPos), col='red')
ggplot() +
    geom_density(data=filter(allSix, nPreds>0), aes(x=nPreds)) +
    geom_density(data=filter(sixTest, nPreds>0), aes(x=nPreds), col='red')

maybeHour <- nThreshSummary(allSix, evCol='hour')
maybeThirty <- nThreshSummary(allSix, evCol='thirty')
View(maybeHour$threshSummary)
View(maybeThirty$threshSummary)
library(ggplot2)
allPdf %>%
    mutate(pred = pTrue > .5) %>%
    ggplot() +
    geom_density(aes(x=pTrue, col=validated)) +
    facet_wrap(pred~drift, scales='free')
library(lubridate)
allSumm <- allPdf %>%
    mutate(estimate = pTrue > 0.5) %>%
    group_by(hour) %>%
    mutate(n=n(),
           nPreds = sum(estimate==TRUE),
           predEvent = nPreds > 0,
           sumPreds = sum(pTrue),
           meanPred = mean(pTrue),
           sumPos = sum(pTrue * (pTrue > thresh)),
           meanPos = sumPos/nPreds,
           ici = mean(All_ici)) %>%
    ungroup()

ggplot(allSumm) +
    geom_point(aes(x=nPreds, y=meanPos, col=validated))

lapply(split(grp, basename(grp$db)), function(x) {
    toNext <- as.numeric(x$start[2:nrow(x)]) - as.numeric(x$start[1:(nrow(x)-1)])
    dur <- as.numeric(x$interval)
    list(toNext=median(toNext), duration=median(dur))
})
# 54, 31, 30, 27, 28, 17 are 22 minute
# there are times where the same detection gets put in mutliple events if
# its right at the end because wav file start times in the db shift by
# a second, so you can get unlucky. This is from event-defined overlap.
sapply(split(allClick, allClick$db), function(x) {
    distCols <- c('UID', 'detectorName', 'BinaryFile', 'Channel', 'eventId')
    nrow(distinct(x[distCols])) == nrow(x)
    # length(unique(x$UID)) == nrow(x)
})

# make 6 min events
allPdf <- allPdf %>%
    group_by(eventId) %>%
    mutate(sixMin=as.numeric(difftime(UTC, min(UTC), units='mins')) / 6,
           sixMin=round(sixMin, 0),
           sixMin=paste0(eventId, '_', sixMin)) %>%
    ungroup()
allSum <- summarisePredEvent(allPdf, thresh=.5, evCol='sixMin')
allSum$train <- FALSE
trainSum <- summarisePredEvent(newPdf, thresh=.5, evCol='sixMin')
trainSum$train <- TRUE
combSum <- bind_rows(allSum, trainSum)
ggplot(combSum) +
    geom_histogram(aes(x=meanPos, fill=train))

