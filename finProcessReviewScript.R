# Fin whale random forest review full workflow
# First version created 2023-12-26
source('finFunctions.R')
## PAMpal PROCESSING ####
library(PAMpal)
# dbFiles needs to be a vector of all the databases you want to process
db <- 'Data/Pre_Validation/'
dbFiles <- list.files(db, full.names=TRUE)
# binDirs needs to be a vector of the corresponding binary file folders for the dbs
bin <- 'Data/ADRIFT/Binaries/'
binDirs <- list.dirs(bin, full.names=TRUE, recursive=FALSE)
# We need to make sure that these are in the same order
basename(dbFiles)
basename(binDirs)

# Here we process the data one drift at a time, then we'll combine them after
# Since drifts have overlapping times, we can't just process everything all at
# once with mode='recording', we have to split it up like this
# Saving intermediate products with "saveRDS" just in case something
# crashes during processing so you don't have to start all over again
studyList <- vector('list', length=length(dbFiles))
for(i in seq_along(dbFiles)) {
    cat(i)
    pps <- PAMpalSettings(db=dbFiles[i], binaries = binDirs[i], sr_hz=200,
                          filterfrom_khz=0, filterto_khz=NULL,
                          winLen_sec=1,
                          threshold=-5,
                          functions=list('ClickDetector'=finClick))
    binList <- pps@binaries$list
    # the data I processed had some binaries that were accidentally
    # copied twice (had a "... (1).pgdf" in the name. We don't want that
    whichCopy <- grepl('\\(1\\)', binList)
    binList <- binList[!whichCopy]
    pps@binaries$list <- binList
    allData <- processPgDetections(pps, mode='recording', id=basename(binDirs[i]))
    studyList[[i]] <- allData
    saveRDS(studyList, file=paste0('SeparateStudyList_to_', i, '.rds'))
}
# saving the raw data once we're done processing so that we have a clean starting
# place in the future. If you make it to here you can delete the "..._to_X.rds"
# files created above
saveRDS(studyList, file='SeparateStudyListAll.rds')
allData <- bindStudies(readRDS('SeparateStudyListAll.rds'))
allData <- calculateICI(allData, time='peakTime')
# connecting wav files, I used the 12kHz decimated files. The "recDirs"
# must be in the same order as the databases
basename(files(allData)$db)
# point this to folder containing your recordings. If they aren't all
# in the same folder, you'll have to type them out in order
recDirs <- list.dirs('D:/FinClick/ADRIFT/Recordings/', recursive=FALSE)
allData <- addRecordings(allData, folder=recDirs)
# we'll save our data at this point because the above steps take some time
saveRDS(allData, file='CombinedStudy_RecIci.rds')
# allData <- readRDS('CombinedStudy_RecIci.rds')

## MODEL PREDICTION ####
# Now we're building our model prediction dataset
allIci <- getICI(allData, type='data')
allClick <- getClickData(allData)
allClick <- joinIci(allClick, allIci)

model <- readRDS('FinModelDec2023.rds')
allPdf <- makePredDf(model, allClick)
# This should work for ADRIFT, will need to be adjusted for others
driftPattern <- '(ADRIFT_[0-9]{1,3})_.*'
allPdf$drift <- gsub(driftPattern, '\\1', basename(allPdf$db))
allPdf$hour <- floor_date(allPdf$UTC, unit='1hour')
allPdf$sixMin <- floor_date(allPdf$UTC, unit='6min')
# the choices of "3" and ".85" seem to work, but could be adjusted in the future
allReview <- markReviewStatus(allPdf, nMin=3, revThresh = .85, evCol=c('drift', 'hour', 'sixMin'))

# we need to reorganize our AcousticStudy to match the hour-long events we are interested in
hourGroups <- bind_rows(lapply(split(allPdf, allPdf$drift), function(x) {
    tRange <- floor_date(range(x$UTC), unit='1hour')
    tSeq <- seq(from=tRange[1], to=tRange[2], by=3600)
    group <- data.frame(start=tSeq, end=tSeq + 3600, db=x$db[1])
    group$id <- paste0(x$drift[1], '_', 1:nrow(group))
    group$drift <- x$drift[1]
    group$hour <- group$start
    group
}))

allHour <- reorganizeEvents(allData, grouping=hourGroups)
# save this as our final dataset to create images from
saveRDS(allHour, 'Data/ADRIFT_AllHourStudy.rds')
# allHour <- readRDS('Data/ADRIFT_AllHourStudy.rds')

# this will go through and create the folders of review images
# it takes a really long time, but should print out your progress along the way
createReviewFolder(allHour, allPdf)

# This shouldn't be needed, but if you want to look at a specific drift/number
# the "id" matches the number rows in the summary graphics for each drift
plotReviewSpec(allHour, allPdf, allReview, drift='ADRIFT_027', id=95, plotAll=TRUE)
