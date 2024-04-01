dbRec <- c(paste0('ADRIFT_0', c(27,28, 49, 50, 51), '_PG_2_02_02_FinClick.sqlite3'),
           paste0('ADRIFT_00', c(2, 3), '_PG_2_02_02_FinClick.sqlite3'),
           paste0('ADRIFT_0', c(30), '_PG_2_02_09_FinClick.sqlite3'))
dataRec <- filter(allData, basename(database) %in% dbRec)
files(dataRec)$db
dataRec <- addRecordings(dataRec, folder=c('D:/FinClick/ADRIFT/Recordings/ADRIFT_002_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_003_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_027_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_028_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_030_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_049_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_050_CENSOR_12kHz/',
                                           'D:/FinClick/ADRIFT/Recordings/ADRIFT_051_CENSOR_12kHz/'))

# saveRDS(dataRec, 'Data/ADRIFT_AllRecStudy.rds')
dataRec <- readRDS('Data/ADRIFT_AllRecStudy.rds')

recPdf <- filter(allPdf, drift %in% paste0('ADRIFT_0', c('02', '03', '27', '28', '30', '49', '50', '51')))
saveRDS(recPdf, file='Data/ADRIFT_AllRecPdf.rds')
recReview <- markReviewStatus(recPdf, nMin=3, revThresh = .85, evCol=c('drift', 'hour', 'sixMin'))

hourGroups <- bind_rows(lapply(split(recPdf, recPdf$drift), function(x) {
    tRange <- floor_date(range(x$UTC), unit='1hour')
    tSeq <- seq(from=tRange[1], to=tRange[2], by=3600)
    group <- data.frame(start=tSeq, end=tSeq + 3600, db=x$db[1])
    group$id <- paste0(x$drift[1], '_', 1:nrow(group))
    group$drift <- x$drift[1]
    group$hour <- group$start
    group
}))

hourRec <- reorganizeEvents(dataRec, grouping=hourGroups)
saveRDS(hourRec, 'Data/ADRIFT_AllRecHourStudy.rds')
hourRec <- readRDS('Data/ADRIFT_AllRecHourStudy.rds')

# acst, pdf,
# df with id, start, end (should be full hour)
plotReviewSpec(hourRec, recPdf, recReview, drift='ADRIFT_003', id=70, length=2, thresh=.5)
# drift, hour, sixMin, reviewStatus, hourStatus, nPreds, id
g <- plotReviewStatus(recReview, drift=='ADRIFT_003')


# pdf %>%  review (has hourStatus-hour, reviewStatus-sixMin)
# pdf %>% hourGroups (has Id) %>% reorgStudy

createReviewFolder(hourRec, recPdf)
# 51, 65, 87, 95 were dropped but fin
plotReviewSpec(hourRec, recPdf, recReview, drift='ADRIFT_027', id=95, plotAll=TRUE)
