# finRF functions
library(PAMpal)
library(dplyr)
library(ggplot2)
library(randomForest)
library(lubridate)
library(patchwork)
library(rfPermute)
library(scales)

fitRFLoop <- function(data, outList, ntree=5e3, mtry=7, style=c('strat', 'sub'), drop=NULL) {
    style <- match.arg(style)
    data$isFin <- factor(data$isFin, levels=c(TRUE, FALSE))


    result <- lapply(outList, function(x) {
        # browser()
        rf <- fitOneRF(data, x, ntree, mtry, style, drop)
        # list(model=rf, pred=makePredDf(rf, test))
        makePredDf(rf, test)
    })
    result
}

fitOneRF <- function(data, test, ntree, mtry=7, style=c('strat', 'sub'), drop=NULL, nodesize=1,
                     seed=11223388, mtryTune=FALSE, rfPermute=FALSE, sampProp=1) {
    style <- match.arg(style)
    data$isFin <- factor(data$isFin, levels=c(TRUE, FALSE))
    nonModelCols <- c('UID', 'UTC', 'Channel', 'angle',
                      'angleError', 'BinaryFile', 'eventLabel',
                      'eventId', 'db', 'species', 'predicted', 'peakTime', 'drift',
                      'Click_Detector_1_ici', 'Click_Detector_0_ici',
                      # adding these bc they are always 0
                      'duration', 'peak2', 'peak3', 'trough', 'trough2',
                      'peakToPeak2', 'peakToPeak3', 'peak2ToPeak3')
    nonModelCols <- c(nonModelCols, drop)
    train <- data[!data$drift %in% test, ]
    train <- PAMpal:::dropCols(train, nonModelCols)
    test <- data[data$drift %in% test, ]
    test <- PAMpal:::dropCols(test, nonModelCols)
    if(mtryTune) {
        whichFin <- which(colnames(train) == 'isFin')
        result <- tuneRF(x=train[-whichFin], y=train[[whichFin]],
                         stepFactor = 1.5, improve=1e-6, ntree=ntree,
                         mtryStart=mtry,
                         strata=train$isFin, sampsize=rep(sum(train$isFin == TRUE), 2))
        return(result)
    }
    switch(style,
           'strat' = {
               set.seed(seed)
               if(rfPermute) {
                   rf <- rfPermute(isFin ~ ., data=train, ntree=ntree, mtry=mtry,
                                   strata=train$isFin, sampsize=rep(sum(train$isFin == TRUE)*sampProp, 2),
                                   nodesize=nodesize)
               } else {
                   rf <- randomForest(isFin ~ ., data=train, ntree=ntree, mtry=mtry,
                                      strata=train$isFin, sampsize=rep(sum(train$isFin == TRUE)*sampProp, 2),
                                      nodesize=nodesize)
               }
           },
           'sub' = {
               set.seed(seed)
               train <- slice_sample(train, n=sum(train$isFin == TRUE)*sampProp, by=isFin)
               if(rfPermute) {
                   rf <- rfPermute(isFin ~ ., data=train, ntree=ntree, mtry=mtry,
                                   nodesize=nodesize)
               } else {
                   rf <- randomForest(isFin ~ ., data=train, ntree=ntree, mtry=mtry,
                                      nodesize=nodesize)
               }

           }
    )
    rf
}

markValidClick <- function(data, dbDir, ici=NULL) {
    dbFiles <- list.files(dbDir, full.names=TRUE, pattern='sqlite3$')
    valData <- lapply(dbFiles, function(x) {
        # this just pulls the event/click info from the DB
        result <- PAMpal:::getDbData(x)
        if(is.null(result) || nrow(result) == 0) {
            return(result)
        }
        result$db <- basename(x)
        result
    }) %>% bind_rows()

    result <- lapply(split(data, data$db), function(x) {
        thisVal <- valData[valData$db == basename(x$db[1]), ]
        # Case if no validated clicks for this DB
        if(is.null(thisVal) || nrow(thisVal) == 0) {
            x$isFin <- FALSE
            return(x)
        }
        # otherwise mark those present in valid DB
        x$isFin <- x$UID %in% thisVal$UID
        x
    })
    result <- bind_rows(result)
    if(is.null(ici)) {
        return(result)
    }
    result <- joinIci(result, ici)
    result
}

joinIci <- function(x, ici) {
    ici <- filter(ici, detectorName == 'All') %>%
        group_by(eventId, Channel) %>%
        mutate(ici_next = c(ici[-1], 0)) %>%
        ungroup() %>%
        rename(ici, 'ici_prev'='ici')

    beforeJoin <- nrow(x)
    x <- left_join(x, ici[c('eventId', 'UID', 'BinaryFile','ici_prev', 'ici_next')])
    afterJoin <- nrow(x)
    if(beforeJoin != afterJoin) {
        warning('Multi join bad')
    }
    x
}

makePredDf <- function(model, test, oob=FALSE) {
    hasTruth <- 'isFin' %in% colnames(test)
    if(hasTruth) {
        test$isFin <- factor(test$isFin, levels=c(TRUE, FALSE))
    }
    # can only use OOB estimates if test==train data, do some checks
    if(!hasTruth ||
       nrow(test) != nrow(model$votes) ||
       sum(model$y == test$isFin) != nrow(test)) {
        if(oob) {
            warning('Data appear different than training data, cannot use OOB')
        }
        oob <- FALSE
    }
    if(oob) {
        probs <- model$votes
    } else {
        probs <- predict(model, newdata=test, type='prob')
    }
    probs <- data.frame(probs)
    colnames(probs) <- c('pTrue', 'pFalse')
    if(hasTruth) {
        test$truth <- test$isFin
    }
    cbind(test, probs)
}

markTF <- function(pdf, thresh=.5, evCol='hour') {
    hasTruth <- 'truth' %in% colnames(pdf)
    pdf$estimate <- factor(pdf$pTrue >= thresh, levels = c(TRUE, FALSE))
    if(hasTruth) {
        pdf$TP <- (pdf$truth == TRUE) & (pdf$estimate == TRUE)
        pdf$FP <- (pdf$truth == FALSE) & (pdf$estimate == TRUE)
        pdf$FN <- (pdf$truth == TRUE) & (pdf$estimate == FALSE)
        pdf$TN <- (pdf$truth == FALSE) & (pdf$estimate == FALSE)
        # pdf$detResult <- NA
        # for(i in 1:nrow(pdf)) {
        #     pdf$detResult[i] <- typeLab(as.logical(pdf$truth[i]), as.logical(pdf$estimate[i]))
        # }
        pdf$detResult <- typeLab(as.logical(pdf$truth), as.logical(pdf$estimate))
        pdf <- pdf %>%
            group_by(.data[[evCol]]) %>%
            mutate(finEvent = any(truth == TRUE),
                   predEvent = any(estimate == TRUE),
                   evResult = typeLab(finEvent, predEvent),
                   nPreds = sum(estimate == TRUE)) %>%
            ungroup()
    } else {
        pdf <- pdf %>%
            group_by(.data[[evCol]]) %>%
            mutate(predEvent = any(estimate == TRUE),
                   nPreds = sum(estimate == TRUE)) %>%
            ungroup()
    }
    pdf
}

summarisePredEvent <- function(pdf, thresh=.5, evCol='eventId') {
    hasTruth <- 'truth' %in% colnames(pdf)
    pdf <- markTF(pdf, thresh) %>%
        group_by(across(all_of(evCol)))
    if(hasTruth) {
        pdf <- pdf %>%
            summarise(n=n(),
                      nFins = sum(truth==TRUE),
                      nPreds = sum(estimate==TRUE),
                      nTP = sum(TP),
                      nFP = sum(FP),
                      nFN = sum(FN),
                      nTN = sum(TN),
                      finEvent = nFins > 0,
                      predEvent = nPreds > 0,
                      evResult = typeLab(finEvent, predEvent),
                      sumPreds = sum(pTrue),
                      meanPred = mean(pTrue),
                      sumPos = sum(pTrue * (pTrue > thresh)),
                      meanPos = sumPos/(sum(TP)+sum(FP)),
                      ici = mean(All_ici))
    } else {
        pdf <- pdf %>%
            summarise(n=n(),
                      nPreds = sum(estimate==TRUE),
                      predEvent = nPreds > 0,
                      sumPreds = sum(pTrue),
                      meanPred = mean(pTrue),
                      sumPos = sum(pTrue * (pTrue > thresh)),
                      meanPos = sumPos/nPreds,
                      ici = mean(All_ici))
    }
    ungroup(pdf)
}

compareThresh <- function(pdf, int=.01, evCol='event') {
    ints <- seq(from=0, to=1, by=int)

    result <- lapply(ints, function(i) {
        tf <- summarisePredEvent(pdf, thresh=i, evCol=evCol)
        results <- table(tf$evResult)
        list(nTP = sum(tf$nTP),
             nFP = sum(tf$nFP),
             nFN = sum(tf$nFN),
             nTN = sum(tf$nTN),
             evTP = results['TP'],
             evFP = results['FP'],
             evFN = results['FN'],
             evTN = results['TN'])
    })
    bind_rows(result) %>%
        mutate(threshold = ints,
               evPrecision = evTP / (evTP + evFP),
               evRecall = evTP / (evTP + evFN),
               precision = nTP / (nTP + nFP),
               recall = nTP / (nTP + nFN))
}

positiveICI <- function(x, evCol='hour', thresh=.5) {
    x <- markTF(x, evCol=evCol, thresh=thresh)
    x <- bind_rows(
        lapply(
            split(x, x[[evCol]]), function(e) {
                bind_rows(lapply(split(e, e$estimate), function(p) {
                    if(nrow(p) == 0) {
                        return(p)
                    }
                    p$predICI <- 0
                    if(p$estimate[1] == FALSE) {
                        return(p)
                    }
                    if(nrow(p) <= 1) {
                        return(p)
                    }
                    p <- arrange(p, UTC)
                    p$predICI[2:nrow(p)] <- as.numeric(
                        difftime(
                            p$UTC[2:nrow(p)],
                            p$UTC[1:(nrow(p)-1)],
                            units='secs')
                    )
                    p
                }
                ))
            }
        )
    )
    x
}

plotComps <- function(models, names, test, thresh=.5) {
    if(length(models) != length(names)) {
        stop('Must name each experiment')
    }
    preds <- lapply(models, function(x) {
        if(is.data.frame(x)) {
            return(x)
        }
        markTF(makePredDf(x, test=test), thresh=thresh)
    })
    names(preds) <- names
    prData <- bind_rows(lapply(preds, function(x) {
        pr_curve(x, truth, pTrue)
    }), .id='Experiment')
    prPlot <- ggplot(prData, aes(x=recall, y=precision, color=Experiment)) +
        geom_path() +
        coord_equal()
    rocData <- bind_rows(lapply(preds, function(x) {
        roc_curve(x, truth, pTrue)
    }), .id='Experiment')
    rocPlot <- ggplot(rocData, aes(x=1-specificity, y=sensitivity, col=Experiment)) +
        geom_path() +
        coord_equal()
    cmPlots <- lapply(seq_along(preds), function(x) {
        cm <- conf_mat(preds[[x]], truth, estimate)
        cmRec <- cm
        cmRec$table <- t(t(cmRec$table) / colSums(cmRec$table))
        cmPrec <- cm
        cmPrec$table <- cmPrec$table / rowSums(cmPrec$table)
        (autoplot(cm, type='heatmap') + ggtitle(names[x])) /
            (autoplot(cmRec, type='heatmap') + ylab('Recall')) /
            (autoplot(cmPrec, type='heatmap') + ylab('Precision'))


    })
    cmPlots <- purrr::reduce(cmPlots, `|`)
    list(prPlot=prPlot, rocPlot=rocPlot, cmPlot=cmPlots)
}




typeLab <- function(truth, estimate) {
    # if(length(truth) > 1) {
    #     return(sapply(seq_along(truth), function(x) {
    #         typeLab(truth[x], estimate[x])
    #     }))
    # }
    # if(isTRUE(truth)) {
    #     if(isTRUE(estimate)) {
    #         return('TP')
    #     }
    #     return('FN')
    # }
    # if(isTRUE(estimate)) {
    #     return('FP')
    # }
    # 'TN'
    scores <- rep(1, length(truth))
    firstPos <- ifelse(truth, 2, 0)
    secondPos <- ifelse(estimate, 1, 0)
    scoreMap <- c('TN', 'FP', 'FN', 'TP')
    scoreMap[scores+firstPos+secondPos]
}

plotThreshComp <- function(x, min=.48) {
    x <- filter(x, threshold > min)
    (ggplot(x, aes(x=threshold)) +
            geom_line(aes(y=nFN, col='FN')) +
            geom_line(aes(y=nFP, col='FP')) +
            geom_line(aes(y=nTP, col='TP')) +
            ggtitle('Detection Level')) +
        (ggplot(x, aes(x=threshold)) +
             geom_line(aes(y=evFN, col='FN')) +
             geom_line(aes(y=evFP, col='FP')) +
             geom_line(aes(y=evTP, col='TP')) +
             ggtitle('Event Level'))
}

# model summary plots
plotModelSummary <- function(pdf, evCol='hour', thresh=.5, title=NULL, file=NULL) {
    pdf <- positiveICI(pdf, evCol=evCol, thresh=thresh)
    detScore <- ggplot(pdf, aes(x=pTrue, fill=isFin)) +
        geom_histogram(binwidth=.025, position='identity', alpha=.7) +
        xlim(thresh-.05, 1.05) +
        ggtitle('Detection Score vs Class') +
        scale_fill_manual(values=scales::hue_pal()(2), breaks=c(FALSE, TRUE))

    evSum <- summarisePredEvent(pdf, thresh=thresh, evCol=evCol)
    results <- table(evSum$evResult)
    resultChar <- paste0(names(results), ':', results, collapse=', ')
    title <- paste(title, resultChar)
    maxPreds <-  max(evSum$nPreds[!evSum$finEvent])+1
    evSum$finEvent <- factor(evSum$finEvent, levels=c(TRUE, FALSE))
    avgScore <- ggplot(evSum, aes(x=meanPos, fill=finEvent)) +
        geom_histogram(binwidth=.025, position='identity', alpha=.7) +
        ggtitle(paste0('Average Positive Score vs Class (thresh=', thresh, ')')) +
        xlim(thresh-.05, 1.05) +
        scale_fill_manual(values=scales::hue_pal()(2), breaks=c(FALSE, TRUE))

    numPos <- ggplot(evSum, aes(x=nPreds, y=meanPos, col=finEvent)) +
        geom_point(size=2) +
        xlim(-.1, maxPreds) +
        ggtitle(paste0('Num Positive vs Average Score by Class (thresh=', thresh, ')')) +
        ylim(thresh-.05, 1.05) +
        scale_color_manual(values=scales::hue_pal()(2), breaks=c(FALSE, TRUE))

    # threshComp <- plotThreshComp(compareThresh(pdf, evCol=evCol))
    ici <- ggplot(pdf, aes(x=predICI, col=detResult)) +
        geom_density() +
        xlim(1, 120) +
        ggtitle(paste0('ICI Between Positive Predictions (thresh=', thresh, ')')) +
        scale_color_manual(values=scales::hue_pal()(2), breaks=c('FP', 'TP'))
    reviewPlots <- checkPctReview(evSum, ns=1:5, tCheck=seq(from=.6, to=.85, by=.05))
    plot <- (detScore | avgScore) /
        (numPos | ici) /
        (reviewPlots) +
        plot_annotation(title)
    if(!is.null(file)) {
        ggsave(file, plot=plot, width=14, height=8, units = 'in')
    }
    plot
}

# plot nPred threshold vs
checkNCutoff <- function(evSum, ns=1:10) {
    data <- bind_rows(lapply(ns, function(n) {
        # dropIx <- evSum$nPreds <= n & evSum$nPreds > 0
        dropIx <- evSum$nPreds == n
        count <- table(evSum$finEvent[dropIx])
        list(
            count = as.numeric(count),
            finEvent = as.logical(names(count)),
            numPredictions = rep(n, length(count))
        )
    }))
    finDropLab <- group_by(data, numPredictions) %>%
        summarise(height=sum(count),
                  nFin=sum(count * finEvent),
                  n=mean(numPredictions))
    ggplot() +
        geom_bar(data=data, aes(x=numPredictions, y=count, fill=finEvent), stat='identity') +
        geom_text(data=finDropLab, aes(x=n, y=height, label=nFin)) +
        ggtitle('Num Events with N Preds (Num Fin Events)')

}

checkPctReview <- function(evSum, ns=1:5, tCheck=c(.5, .6, .65, .7, .75, .8)) {
    nFins <- sum(evSum$finEvent == TRUE)
    evSum$finEvent <- as.logical(evSum$finEvent)
    nEvents <- nrow(evSum)
    cutOffPlot <- checkNCutoff(evSum, ns=1:10)
    evSum <- evSum[evSum$nPreds > 0, ]
    pctData <- bind_rows(
        lapply(
            ns, function(n) {
                dropIx <- evSum$nPreds < n
                thisSum <- evSum[!dropIx, ]
                bind_rows(lapply(tCheck, function(t) {
                    revIx <- thisSum$meanPos <= t
                    result <- list(n=n,
                                   revThresh = t,
                                   revFin = sum(thisSum$finEvent[revIx]),
                                   revNotFin = sum(!thisSum$finEvent[revIx]),
                                   acceptFin = sum(thisSum$finEvent[!revIx]),
                                   acceptNotFin = sum(!thisSum$finEvent[!revIx])
                    )
                    result$revPct <- (result$revFin + result$revNotFin) / nEvents
                    result
                }))
            }))
    pctPlot <- ggplot(pctData) +
        geom_tile(aes(x=n, y=revThresh, width=1, height=.05, fill=revPct)) +
        geom_text(aes(x=n, y=revThresh, label=round(revPct, 2))) +
        scale_fill_gradientn(colors=rev(scales::viridis_pal()(25))) +
        guides(fill='none') +
        ggtitle(paste0('Pct of Data To Review (', nEvents, ' events)')) +
        xlab('Min Predictions') +
        ylab('Review Under Thresh')
    falsePosPlot <- ggplot(pctData) +
        geom_tile(aes(x=n, y=revThresh, fill=acceptNotFin, width=1, height=.05)) +
        geom_text(aes(x=n, y=revThresh, label=acceptNotFin)) +
        scale_fill_gradientn(colors=rev(scales::viridis_pal()(25))) +
        guides(fill='none') +
        ggtitle(paste0('Num FP Accepted (', nFins, ' fin events)')) +
        xlab('Min Predictions') +
        ylab('Review Under Thresh')

    cutOffPlot + pctPlot + falsePosPlot +
        plot_layout(widths = c(1, 1.4, 1.4))
}

reorganizeEvents <- function(x, grouping, progress=TRUE) {# REGROUPER DEVEL
    grouping <- PAMpal:::checkGrouping(grouping)
    acousticEvents <- vector('list', length = nrow(grouping))
    evName <- as.character(grouping$id)

    colsToDrop <- c('Id', 'comment', 'sampleRate', 'detectorName', 'parentID',
                    'sr', 'callType', 'newUID')
    names(acousticEvents) <- evName
    noDetEvent <- character(0)

    binData <- getDetectorData(x, measures=FALSE)
    for(d in seq_along(binData)) {
        binData[[d]]$callType <- names(binData)[d]
    }
    binData <- lapply(binData, function(d) split(d, d$detectorName))
    names(binData) <- NULL
    binData <- unlist(binData, recursive = FALSE)
    binData <- squishList(binData)
    hasDb <- 'db' %in% colnames(grouping)
    if(hasDb) {
        binData <- lapply(binData, function(x) {
            split(x, basename(x$db))
        })
        binData <- purrr::transpose(binData)
        hasGrpDb <- basename(grouping$db) %in% names(binData)
        if(!all(hasGrpDb)) {
            missDb <- unique(basename(grouping$db[!hasGrpDb]))
            warning('Databases ', paste0(missDb, collapse=', '),
                    ' are in the "grouping" file but not the AcousticStudy ',
                    'to reorganize')
            grouping <- grouping[hasGrpDb, ]
        }
    }

    if(progress) {
        cat('Organizing events...\n')
        pb <- txtProgressBar(min=0, max=length(acousticEvents), style=3)
    }
    for(i in seq_along(acousticEvents)) {
        if(hasDb) {
            thisData <- binData[[basename(grouping$db[i])]]
            thisData <- lapply(thisData, function(x) {
                data <- x[PAMpal:::withinLHS(x$UTC, grouping$interval[i]), ]
                if(nrow(data) == 0) return(NULL)
                data
            })
        } else {
            thisData <- lapply(binData, function(x) {
                data <- x[PAMpal:::withinLHS(x$UTC, grouping$interval[i]), ]
                if(nrow(data) == 0) return(NULL)
                data
            })
        }
        # Check possible DBs by start/end time of events in sa list earlier
        hasDat <- sapply(thisData, function(x) !is.null(x))
        # looks unnecessary but subset breaks if thisData is already list()
        if(!any(hasDat)) {
            thisData <- list()
        } else {
            thisData <- thisData[hasDat]
        }

        oldEvents <- lapply(thisData, function(x) unique(x$eventId)) %>%
            unlist(recursive=FALSE) %>% unique()
        if(length(oldEvents) == 0) {
            next
        }
        fullBins <- lapply(events(x[oldEvents]), function(e) {
            files(e)$binaries
        }) %>% unlist() %>% unique()

        binariesUsed <- lapply(thisData, function(x) unique(x$BinaryFile)) %>%
            unlist(recursive = FALSE) %>% unique()
        binariesUsed <- sapply(binariesUsed, function(b) {
            full <- grep(b, fullBins, value=TRUE)
            if(length(full) == 0) {
                warning('No matching binary for ', b)
            }
            full
        })
        # Check and warning here for empty event
        if(length(thisData) == 0) {
            noDetEvent <- c(noDetEvent, names(acousticEvents)[i])
        }
        thisData <- lapply(thisData, function(x) {
            x$BinaryFile <- basename(x$BinaryFile)
            thisType <- unique(x$callType)
            x <- PAMpal:::dropCols(x, colsToDrop)
            attr(x, 'calltype') <- thisType
            x
        })
        thisSr <- grouping$sr[[i]] # XXX pull from settings probably?
        oldSr <- lapply(events(x[oldEvents]), function(e) {
            settings(e)$sr
        }) %>% unlist() %>% unique()

        # thisDb <- allDbs[grepl(basename(grouping$db[i]), allDbs)] # XXX this was previously matched, pull from evs
        oldDb <- lapply(events(x[oldEvents]), function(e) {
            files(e)$db
        }) %>% unlist() %>% unique()
        acousticEvents[[i]] <-
            PAMpal:::AcousticEvent(id=evName[i], detectors = thisData, settings = list(sr = oldSr[1]),# source = thisSource),
                                   files = list(binaries=binariesUsed, db=oldDb, calibration=NULL))
        ancillary(acousticEvents[[i]])$grouping <- grouping[i, ]
        if(progress) {
            setTxtProgressBar(pb, value=i)
        }
    }
    nullEv <- sapply(acousticEvents, is.null)
    acousticEvents <- acousticEvents[!nullEv]
    events(x) <- acousticEvents
    ancillary(x)$grouping <- grouping
    x
}

markReviewStatus <- function(x, nMin=1, revThresh=.75, evCol=c('drift', 'eventId', 'sixMin'), predThresh=.5) {
    # check if this is summarised first
    if(!all(c('nPreds', 'meanPos') %in% colnames(x))) {
        x <- summarisePredEvent(x, evCol=evCol, thresh=predThresh)
    }
    x$reviewStatus <- NA
    dropIx <- x$nPreds < nMin
    revIx <- !dropIx & (x$meanPos <= revThresh)
    accIx <- !dropIx & (x$meanPos > revThresh)
    x$reviewStatus[dropIx] <- 'Dropped'
    x$reviewStatus[revIx] <- 'Review'
    x$reviewStatus[accIx] <- 'Accepted'
    if(all(c('drift', 'hour') %in% colnames(x))) {
        x <- x %>%
            group_by(drift, hour) %>%
            mutate(nReview = sum(reviewStatus == 'Review'),
                   nAccept = sum(reviewStatus == 'Accepted'),
                   hourStatus = ifelse(nReview > 0, 'Review', 'Dropped'),
                   hourStatus = ifelse(nAccept > 0, 'Accepted', hourStatus)) %>%
            ungroup()
        # THIS BAD WHEN NOT FULL RANGE IN REVIEW
        # x <- group_by(x, drift) %>%
        #     arrange(hour) %>%
        #     mutate(id = paste0(drift, '_', as.numeric(as.factor(hour)))) %>%
        #     ungroup()
        hourGroups <- bind_rows(lapply(split(x, x$drift), function(d) {
            # tRange <- floor_date(range(d$UTC), unit='1hour')
            tRange <- range(d$hour)
            tSeq <- seq(from=tRange[1], to=tRange[2], by=3600)
            group <- list(hour=tSeq)
            group$id <- paste0(d$drift[1], '_', 1:length(tSeq))
            group$drift <- d$drift[1]
            group
        }))
        x <- left_join(x, hourGroups, by=c('drift', 'hour'))
    }
    x
}

# study is for plotGram, pdf is for UIDs to filter to, review has
# start/end and event id for small-ing the study
# review has six minute and hour level review status markings
plotReviewSpec <- function(acst, pdf, review, drift, id=1, length=2, thresh=.5,
                           sr=6e3, wl=4096, hop=.125, freqRange=c(0, .1),
                           detections=c('click'), color='red', size=3, filename=NULL,
                           outDir=NULL, plotAll = FALSE, progress=TRUE) {
    if(!all(c('reviewStatus', 'hourStatus') %in% colnames(pdf))) {
        pdf <- left_join(pdf,
                         distinct(review[c('drift', 'hour', 'sixMin', 'id', 'reviewStatus', 'hourStatus')]),
                         by=c('drift', 'hour', 'sixMin'))
    }
    review <- distinct(select(review, any_of(c('id', 'start', 'end', 'drift','hour'))))
    if(!all(c('start', 'end') %in% colnames(review))) {
        review$start <- review$hour
        review$end <- review$start + 3600
    }
    review <- review[review$id == paste0(drift, '_', id), ]
    if(nrow(review) == 0) {
        warning('No event ID ', drift, '_', id)
        return(NULL)
    }
    if(nrow(review) > 1) {
        warning('More than one row matches event ID ', drift, '_', id)
        return(NULL)
    }
    acst <- acst[review$id]

    pdf <- filter(pdf, id == review$id)
    if(isFALSE(plotAll)) {
        pdf <- filter(pdf, pTrue > thresh, reviewStatus != 'Dropped')
    }
    acst <- filter(acst, UID %in% pdf$UID)
    revTimes <- seq(from=review$start, to=review$end, by=60 * length)
    revPlan <- data.frame(start=revTimes[1:(length(revTimes)-1)], end=revTimes[2:length(revTimes)])
    revPlan$id <- 1:nrow(revPlan)
    pdf$revMatch <- rep(NA, nrow(pdf))
    for(i in 1:nrow(revPlan)) {
        thisMatch <- pdf$UTC >= revPlan$start[i] &
            pdf$UTC < revPlan$end[i]
        pdf$revMatch[thisMatch] <- revPlan$id[i]
    }
    if(anyNA(pdf$revMatch)) {
        noMatch <- is.na(pdf$revMatch)
        warning(sum(noMatch), ' detections were outside of event time range')
        pdf <- pdf[!noMatch, ]
    }
    plotIx <- unique(pdf$revMatch)
    nPlots <- length(plotIx)
    revPlan <- filter(revPlan, id %in% plotIx)
    if(progress) {
        pb <- txtProgressBar(min=0, max=nPlots, style=3)
        setTxtProgressBar(pb, value=0)
    }
    if(is.null(filename)) {
        filename <- paste0(review$id, '.png')
    }
    if(is.null(outDir)) {
        outDir <- '.'
    }
    if(!dir.exists(outDir)) {
        dir.create(outDir)
    }
    png(filename=file.path(outDir, filename), width=1600, height=600 * nPlots, units='px')
    on.exit(dev.off())
    # par(cex.main=2, adj=0)
    par(mfrow=c(nPlots, 1), cex.axis=2, cex.lab=2)
    # can only set this if only plotting click points
    # otherwise it overlaps with contour lwd
    if(identical(detections, 'click')) {
        lwd <- 2
    } else {
        lwd <- NULL
    }
    for (i in 1:nPlots) {
        tryCatch({plotGram(acst, evNum=1, start=revPlan$start[i], end=revPlan$end[i],
                           sr=sr, wl=wl, hop=hop, freqRange=freqRange, detections=detections,
                           detCol=color, size=size,
                           title=NA, lwd=lwd)
            title(paste0('Review ID: ', review$id), cex.main=2, adj=0)
        },
        warning = function(w) {
            cat('\n', 'Issue in', review$id, 'at time', revPlan$start[i], '\n')
            plotGram(acst, evNum=1, start=revPlan$start[i], end=revPlan$end[i],
                     sr=sr, wl=wl, hop=hop, freqRange=freqRange, detections=detections,
                     detCol=color, size=size,
                     title=NA, lwd=lwd)
            title(paste0('Review ID: ', review$id), cex.main=2, adj=0)
        })
        if(progress) {
            setTxtProgressBar(pb, value=i)
        }
    }
    revPlan
}
# drift, hour, sixMin, reviewStatus, hourStatus, nPreds, id
plotReviewStatus <- function(x, drift, filename=NULL) {
    x <- filter(x, drift == drift)
    x <- x %>%
        mutate(xmin = minute(sixMin)+.15,
               xmax = xmin+6-.15,
               ymin = hour+360,
               ymax = ymin + 3600-360)
    x$id <- gsub(paste0(drift, '_'), '', x$id)
    nHours <- length(unique(x$hour))

    g <- ggplot() +
        geom_rect(data=x, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, col=reviewStatus, fill=nPreds)) +
        geom_rect(data=filter(x, hourStatus == 'Accepted'),
                  aes(xmin=60, xmax=62, ymin=ymin, ymax=ymax), fill=hue_pal()(3)[2]) +
        geom_rect(data=filter(x, hourStatus == 'Accepted'),
                  aes(xmin=-2, xmax=0, ymin=ymin, ymax=ymax), fill=hue_pal()(3)[2]) +
        geom_rect(data=filter(x, hourStatus == 'Review'),
                  aes(xmin=60, xmax=62, ymin=ymin, ymax=ymax), fill=hue_pal()(3)[1]) +
        geom_rect(data=filter(x, hourStatus == 'Review'),
                  aes(xmin=-2, xmax=0, ymin=ymin, ymax=ymax), fill=hue_pal()(3)[1]) +
        geom_text(data=x, aes(x=-4, y=hour + 3600*.5, label=id), size=3) +
        # geom_tile(aes(x=xmin, y=ymin, col=reviewStatus, fill=nPreds), alpha=.5) +
        scale_fill_gradient(low='white', high='black', limits=c(0, 10)) +
        scale_color_manual(values=c('lightgrey', scales::hue_pal()(3)[c(2,1)]), breaks=c('Dropped', 'Accepted', 'Review')) +
        theme_bw() +
        ggtitle(drift) +
        scale_y_datetime(expand=c(0, 0))
    if(is.null(filename)) {
        filename <- paste0(drift, '_ReviewSumamry.png')
    }
    if(is.na(filename)) {
        return(g)
    }
    ggsave(plot=g, filename=filename, width=2000,
           height = 30 * nHours + 200, units='px', dpi=300, limitsize=FALSE)
    g
}

createReviewFolder <- function(acst, pdf, drift=NULL, folder=NULL, nMin=3, revThresh=.85, specLen=2, progress=TRUE,
                               status=c('Review', 'Accepted')) {
    if(is.null(drift)) {
        drift <- unique(pdf$drift)
    }
    if(is.null(folder)) {
        folder <- paste0(drift, '_ReviewProducts')
    }
    if(any(!dir.exists(folder))) {
        sapply(folder[!dir.exists(folder)], dir.create)
    }
    pdf <- pdf[pdf$drift %in% drift, ]
    review <- markReviewStatus(pdf, nMin=nMin, revThresh=revThresh, evCol=c('drift', 'hour', 'sixMin'))
    pdf <- left_join(pdf,
                     distinct(review[c('drift', 'hour', 'sixMin', 'id', 'reviewStatus', 'hourStatus')]),
                     by=c('drift', 'hour', 'sixMin'))
    for(d in seq_along(drift)) {
        if(progress) {
            cat('Starting drift', drift[d], '...\n')
        }
        thisPdf <- pdf[pdf$drift == drift[d], ]
        thisReview <- review[review$drift == drift[d], ]
        revSumPlot <- plotReviewStatus(thisReview,
                                       drift=drift[d],
                                       filename = file.path(folder[d], paste0(drift[d], '_ReviewSummary.png')))
        thisReview <- distinct(select(thisReview, any_of(c('id', 'start', 'end', 'drift','hour', 'hourStatus'))))
        # don't fail if nothing to review
        nToReview <- sum(thisReview$hourStatus %in% status)
        if(nToReview == 0) {
            revCsv <- distinct(thisReview[c('drift', 'hour', 'hourStatus', 'id')])
            names(revCsv)[3] <- 'reviewStatus'
            revCsv$manualReview <- NA
            write.csv(revCsv, file = file.path(folder[d], paste0(drift[d], '_ReviewLog.csv')), row.names=FALSE)
            next
        }
        if(progress) {
            pb <- txtProgressBar(min=0, max=nToReview, style=3)
            pbix <- 1
        }
        for(i in 1:nrow(thisReview)) {
            if(!thisReview$hourStatus[i] %in% status) {
                next
            }
            plotReviewSpec(acst, thisPdf, thisReview, drift=drift[d],
                           id=as.numeric(gsub(paste0(drift[d], '_'), '', thisReview$id[i])), length=specLen,
                           filename=paste0(thisReview$id[i], '.png'),
                           outDir=file.path(folder[d], thisReview$hourStatus[i]), progress=FALSE)
            if(progress) {
                setTxtProgressBar(pb, value=pbix)
                pbix <- pbix + 1
            }
        }
        if(progress) {
            close(pb)
            cat('\n')
        }
        revCsv <- distinct(thisReview[c('drift', 'hour', 'hourStatus', 'id')])
        names(revCsv)[3] <- 'reviewStatus'
        revCsv$manualReview <- NA
        write.csv(revCsv, file = file.path(folder[d], paste0(drift[d], '_ReviewLog.csv')), row.names=FALSE)
    }
}

finClick <- function(data, sr_hz='auto', winLen_sec=1, threshold=-5, calibration=NULL) {
    # List names of all required packages here to guarantee they will be
    # installed when others try to use your function, ex:
    # packageList <- c('seewave', 'signal')
    packageList <- c()
    for(p in packageList) {
        if(!require(p, character.only = TRUE)) {
            install.packages(p)
            require(p, character.only = TRUE)
        }
    }
    # prepping output
    result <- list()
    # Do same stuff for each channel
    for(chan in 1:ncol(data$wave)) {
        # create storage for this channels outputs
        thisResult <- list()
        if(sr_hz == 'auto') {
            sr <- data$sr
        } else {
            sr <- sr_hz
        }
        # browser()
        ##### DO STUFF AFTER HERE #######
        thisWave <- data$wave[, chan]
        fftSize <- round(sr * winLen_sec, 0)
        fftSize <- fftSize + (fftSize %% 2)
        wherePeak <- which.max(abs(thisWave))
        starts <- seq(from=-.5, to=.5, by=.05) * fftSize + wherePeak - fftSize/4
        windows <- vector('list', length=length(starts))
        for(s in seq_along(starts)) {
            start <- max(starts[s], 1)
            end <- min(starts[s] + fftSize / 2, length(thisWave))
            if(start > length(thisWave) ||
               end < 1) {
                windows[[s]] <- list(dB=NA_real_, freq=0, time=starts[s]/sr)
                next
            }
            fftVec <- start:end
            winWave <- thisWave[start:end]
            if(length(winWave) < fftSize/2) {
                winWave <- c(winWave, rep(0, fftSize/2 - length(winWave)))
            }
            if(all(winWave == 0)) {
                windows[[s]] <- list(dB=NA_real_, freq=0, time=starts[s]/sr)
                next
            }
            spec <- seewave::spec(winWave, wl=fftSize/2, f=sr, norm=FALSE, correction='amplitude', plot=FALSE)
            dB <- 20*log10(spec[, 2])
            whichMax <- which.max(dB)[1]
            windows[[s]] <- list(dB=dB[whichMax], freq=spec[whichMax, 1], time=starts[s]/sr)
        }
        windows <- bind_rows(windows)
        windows$dB[is.na(windows$dB)] <- min(windows$dB, na.rm=TRUE) - 10
        windows$dB <- windows$dB - max(windows$dB)
        useCalc <- windows[windows$dB > threshold, ]
        peak <- useCalc$freq[which.max(useCalc$dB)[1]]
        center <- mean(range(useCalc$freq))
        if(nrow(useCalc) > 1) {
            lm <- stats::lm(freq ~ time, data=useCalc)
            slope <- unname(lm$coefficients[2])
        } else {
            slope <- 0
        }
        # thisSpec <- spec(thisWave, f=data$sr, wl=512, norm=FALSE, correction='amplitude', plot=FALSE)
        # `spec` has problems with high samplerates producing NA frequencies (integer overflow)
        # so we recreate the frequency vector in a more stable way
        # freq <- seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec))
        # dB <- 20*log10(thisSpec[,2])

        # Store results you want output in `thisResult`:
        # wavMax <- max(thisWave)
        # thisResult$wavMax <- wavMax
        #### STORE THIS CHANNEL ####
        result[[chan]] <- list(freqSlope=slope, centerFreq=center, cenPeakDiff=unname(center-peak))
    }
    # Combine results for each channel into a single dataframe
    # Theres no need to label each channel, this handles in earlier steps within PAMpal
    result <- bind_rows(result)
    result
    # starts
}
