library(nnls)
library(Hmisc)
library(parallel)

options(width = 220)
options(scipen = 999)
options(digits = 4)

printDims = function(X, text){print(paste(text, ':', paste(dim(X), collapse = 'x')))}

rye.scale = function(X = NULL) {
    return(apply(X, 2, function(i){i = i - min(i); i / max(i)}))
}

rye.populationMeans = function(X = NULL, fam = NULL, alpha = NULL, weight = NULL, fn = median, referenceGroups = NULL) {

    ## Find the mean of each reference popuatlion
    if (!is.null(referenceGroups)) {
        means = aggregate(X, by = list(referenceGroups[fam[ , 'population']]), fn)
    } else {
        means = aggregate(X, by = list(fam[ , 'population']), fn)
    }

    ## Reformat
    rownames(means) = means[ , 1]
    means = means[ , 2:ncol(means)]

    ## Apply shrinkage by given method and alpha
    means = apply(means, 2, function(i)  i + (((1/2 - i)**2) * (((i > 1/2) * -1) + (i <= 1/2)) * alpha))
    
    ## Weight each feature
    means = t(t(means) * weight)

    return(means)
}

rye.predict = function(X = NULL, means = NULL, weight = NULL, referenceGroups = NULL) {
    
    estimates = t(apply(t(t(X) * weight), 1, function(i){c = coef(nnls(A= as.matrix(t(means)), b = i)); c / sum(c)}))
    colnames(estimates) = rownames(means)
    
    if (!is.null(referenceGroups)) {
        estimates = do.call(cbind, lapply(unique(referenceGroups), function(i) apply(estimates[ , names(referenceGroups)[referenceGroups == i], drop = FALSE], 1, sum)))
        colnames(estimates) = unique(referenceGroups)
    }
    return(estimates)
}

ruye.squaredError = function(expected = NULL, predicted = NULL) {
    return((expected - predicted) ** 2)
}

rye.absoluteError = function(expected = NULL, predicted = NULL) {
    return(abs(expected - predicted))
}

rye.gibbs = function(X = NULL, fam = NULL, referenceGroups = NULL, 
                 alpha = NULL, optimizeAlpha = TRUE,
                 weight = NULL, optimizeWeight = TRUE,
                 iterations = 100, sd = 0.0001) {

    pops = names(alpha)

    ## Assume the correct ref assignment is 100% their population
    expected = matrix(0, nrow = nrow(X), ncol = length(pops), dimnames = list(rownames(fam), pops))
    expected[fam[ , c('id', 'population')]] <- 1

    ## Make each pop its own group if groups aren't given
    if (is.null(referenceGroups)) {
        referenceGroups = pops
        names(referenceGroups) = pops
    }
    
    expected = matrix(0, nrow = nrow(X), ncol = length(unique(referenceGroups)), dimnames = list(rownames(fam), unique(referenceGroups)))
    expected[cbind(fam[ , 'id'], referenceGroups[fam[ , 'population']])] = 1
    
    fam = cbind(fam, referenceGroups[fam[ , 'population']])
    colnames(fam)[ncol(fam)] = 'group'

    ## Get the starting error
    means = rye.populationMeans(X = X, fam = fam, alpha = alpha, weight = weight, referenceGroups = referenceGroups)[pops, ]
    predicted = rye.predict(X = X, means = means, weight = weight, referenceGroups = referenceGroups)
    oldError = rye.absoluteError(expected = expected, predicted = predicted)
    oldError = cbind(apply(oldError, 1, mean))
    oldError = aggregate(oldError, by = list(fam[ , 'group']), mean)
    oldError = oldError[ , -1]
    oldError = mean(oldError)
    
    ## Return values
    minError = oldError
    minParams = list(minError, alpha, weight, means, predicted)

    ## Momentum between iterations
    alphaMomentum = rep(0, length(alpha))
    weightMomentum = rep(0, length(weight))
    momentum = 1/10
    
    for (iteration in seq(iterations)) {
        
        ## Pick new alpha and weight for this iteration
        newAlpha = alpha
        if (optimizeAlpha) {
            toUpdate = sample(seq(length(newAlpha)))[1]
            newAlpha[toUpdate] = newAlpha[toUpdate] + rnorm(n = 1, sd = (abs(newAlpha[toUpdate]) + 0.001) * sd) + alphaMomentum[toUpdate]
            newAlpha[newAlpha < 0] = 0
        }

        newWeight = weight
        if (optimizeWeight) {
            toUpdate = sample(seq(length(newWeight)))[1]
            newWeight[toUpdate] = newWeight[toUpdate] + rnorm(n = 1, sd = (newWeight[toUpdate] + 0.001) * sd) + weightMomentum[toUpdate]
            newWeight[newWeight < 0] = 0
        }

        ## Find the new errors
        means = rye.populationMeans(X = X, fam = fam, alpha = newAlpha, weight = newWeight, referenceGroups = referenceGroups)[pops, ]
        predicted = rye.predict(X = X, means = means, weight = newWeight, referenceGroups = referenceGroups)
        newError = rye.absoluteError(expected = expected, predicted = predicted)
        newError = cbind(apply(newError, 1, mean))
        newError = aggregate(newError, by = list(fam[ , 'group']), mean)
        newError = newError[ , -1]
        newError = mean(newError)
        
        ## Find the jump odds
        odds = pnorm(newError, mean = oldError, sd = oldError / 1000)
        odds = c(1 - odds, odds)

        ## If this is the best error we've seen, then keep it
        if (newError < minError) {
            minError = newError
            minParams = list(minError, alpha, weight, means, predicted)
        }

        ## See if we jump
        if (runif(n = 1, min = 0, max = 1) < odds[1]) {
            oldError = newError
            alphaMomentum = (alphaMomentum / 2) + ((newAlpha - alpha) * momentum)
            weightMomentum = (weightMomentum / 2) + ((newWeight - weight) * momentum)
            alpha = newAlpha
            weight = newWeight
        }
        
    }
    
    return(minParams)
    
}


rye.optimize = function(X = NULL, fam = NULL,
                    referencePops = NULL, referenceGroups = NULL,
                    alpha = NULL, optimizeAlpha = TRUE,
                    weight = NULL, optimizeWeight = TRUE,
                    iterations = 100, rounds = 25, threads = 1, startSD = 0.005, endSD = 0.001,
                    populationError = FALSE) {

    ## Pull out the reference PCs
    referenceFAM = fam[fam[ , 'population'] %in% referencePops , ]
    referenceX = X[rownames(referenceFAM), ]

    ## Start with the shrinking at 0.05 for all pops by default
    if (is.null(alpha)) {
        alpha = rep(0.001, length(referencePops))
    } 
    names(alpha) = referencePops

    ## Weights
    if (is.null(weight)) {
        weight = 1 / seq(ncol(X))
    }

    allErrors = c()

    for (round in seq(rounds)) {

        sd = startSD - (startSD - endSD) * log(round)/log(rounds)
        if (threads > 1) {
            params = mclapply(seq(threads), function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
                                                                  iterations = iterations,
                                                                  alpha = alpha, weight = weight, sd = sd,
                                                                  optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight), mc.cores = 32)
        } else {
            params = lapply(seq(threads), function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
                                                                  iterations = iterations,
                                                                  alpha = alpha, weight = weight, sd = sd,
                                                                  optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight))
        } 


        errors = unlist(lapply(params, function(i) i[[1]]))
        
        bestError = which.min(errors)
        meanError = mean(errors)
        print(paste0('Round ', round, '/', rounds, ' Mean error: ', sprintf("%.6f", meanError), ', Best error: ', sprintf('%.6f', errors[bestError])))
        
        bestParams = params[[bestError]]
        alpha = bestParams[[2]]
        weight = bestParams[[3]]

        allErrors = c(allErrors, errors[bestError])

        ## See if our error hasn't decreased substantially in 5 rounds
        if (round > 5) {
            errorChange = allErrors[(round - 5):round]
            errorChange = max(errorChange) - min(errorChange)
            if (errorChange <= 0.000025) {
                break
            }
        }
        
    }

    return(bestParams)
}

