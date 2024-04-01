# testing to see if adding partial 28Hz data lets us predict correctly
fakeClick <- clickData
(fakeClick %>%
    filter(drift == 'ADRIFT_027') %>%
    ggplot(aes(x=UTC, y=peak, col=isFin)) +
    geom_point()) %>%
ggplotly()

# 7-28 22:00 cutoff to split AD27

is27 <- fakeClick$drift == 'ADRIFT_027'
firstHalf <- fakeClick$UTC <= as.POSIXct('2022-07-28 22:00:00', tz='UTC')
sum(is27 & firstHalf)
sum(is27 & !firstHalf)
fakeClick$drift[is27 & firstHalf] <- 'ADRIFT_027.1'
fakeClick$drift[is27 & !firstHalf] <- 'ADRIFT_027.2'

fakeSet <- c('ADRIFT_003', 'ADRIFT_028', 'ADRIFT_027.2', 'ADRIFT_047', 'ADRIFT_053')
fakeModel <- fitOneRF(fakeClick, test=fakeSet,
                      mtry=5, ntree=4e3, style='strat', drop=c('rand'),
                         seed=882111)
fakePred <- makePredDf(fakeModel, test=fakeClick)

fakePred$test <- fakePred$drift %in% fakeSet
fakePred$hour <- paste0(fakePred$drift, '_', floor_date(fakePred$UTC, unit='1hour'))
fakePred <- positiveICI(fakePred, evCol='hour', thresh=.5)
fakeTest <- fakePred$drift %in% fakeSet
table(fakePred$truth[fakeTest], fakePred$estimate[fakeTest])
table(fakePred$detResult[fakeTest], fakePred$drift[fakeTest])
# result: it does not
ggplot(fakePred, aes(x=peak, col=detResult)) +
    geom_density() +
    facet_wrap(~test, ncol=1) +
    xlim(0, .05)


fakePred %>%
    filter(drift %in% c('ADRIFT_027.1', 'ADRIFT_027.2')) %>%
    ggplot(aes(x=pTrue, y=peak, col=detResult)) +
    geom_point()
