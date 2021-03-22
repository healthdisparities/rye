library(nnls)
library(Hmisc)
library(parallel)

options(width = 220)
options(scipen = 999)
options(digits = 4)

printDims = function(X, text){print(paste(text, ':', paste(dim(X), collapse = 'x')))}
hmean = function(i) {1/mean(1/i)}

scale = function(X = NULL) {
    return(apply(X, 2, function(i){i = i - min(i); i / max(i)}))
}

populationMeans = function(X = NULL, fam = NULL, alpha = NULL, weight = NULL, fn = median) {

    ## Find the mean of each reference popuatlion
    means = aggregate(X, by = list(fam[ , 'population']), FUN = fn)
    rownames(means) = means[ , 1]
    means = means[ , 2:ncol(means)]

    ## Apply shrinkage by given method and alpha
    means = apply(means, 2, function(i)  i + (((1/2 - i)**2) * (((i > 1/2) * -1) + (i <= 1/2)) * alpha))
    
    ## Weight each feature
    means = t(t(means) * weight)

    return(means)
}

predict = function(X = NULL, means = NULL, weight = NULL, referenceGroups = NULL) {
    estimates = t(apply(t(t(X) * weight), 1, function(i){c = coef(nnls(A= as.matrix(t(means)), b = i)); c / sum(c)}))
    colnames(estimates) = rownames(means)
    if (!is.null(referenceGroups)) {
        estimates = do.call(cbind, lapply(unique(referenceGroups), function(i) apply(estimates[ , names(referenceGroups)[referenceGroups == i], drop = FALSE], 1, sum)))
        colnames(estimates) = unique(referenceGroups)
    }
    return(estimates)
}

squaredError = function(expected = NULL, predicted = NULL) {
    #return((expected - predicted) ** 2)
        return(abs(expected - predicted))
}

meanSquaredError = function(expected = NULL, predicted = NULL) {
    error = squaredError(expected = expected, predicted = predicted)
    return(sum(error) / nrow(error)) 
}

populationSquaredError = function(expected = NULL, predicted = NULL, fam = NULL) {
    error = squaredError(expected = expected, predicted = predicted)
    error = aggregate(error, by = list(fam[rownames(error), 'population']), mean)
    rownames(error) = error[ , 1]
    error = error[ , -1]
    return(sum(error) / nrow(error))
}

rye.gibbs = function(X = NULL, fam = NULL, referenceGroups = NULL, 
                 alpha = NULL, optimizeAlpha = TRUE,
                 weight = NULL, optimizeWeight = TRUE,
                 iterations = 100, sd = 0.0001) {

    iteration = 1
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
    means = populationMeans(X = X, fam = fam, alpha = alpha, weight = weight)[pops, ]
    predicted = predict(X = X, means = means, weight = weight, referenceGroups = referenceGroups)
    oldError = squaredError(expected = expected, predicted = predicted)
    oldError = cbind(apply(oldError, 1, sum))
    oldError = aggregate(oldError, by = list(fam[ , 'group']), mean)
    oldError = oldError[ , -1]
    oldError = mean(oldError)
    
    ## Return values
    minError = oldError
    minParams = list(minError, alpha, weight, means, predicted)

    ## Momentum between iterations
    alphaMomentum = rep(0, length(alpha))
    weightMomentum = rep(0, length(weight))
    momentum = 1
    
    while (iteration <= iterations) {
        
        ## Pick new alpha and weight for this iteration
        newAlpha = alpha
        if (optimizeAlpha) {
            toUpdate = sample(seq(length(newAlpha)))[1]
            newAlpha[toUpdate] = newAlpha[toUpdate] + rnorm(n = 1, sd = (newAlpha[toUpdate] + 0.001) * sd) + alphaMomentum[toUpdate]
            newAlpha[newAlpha < 0] = 0
        }

        newWeight = weight
        if (optimizeWeight) {
            toUpdate = sample(seq(length(newWeight)))[1]
            newWeight[toUpdate] = newWeight[toUpdate] + rnorm(n = 1, sd = (newWeight[toUpdate] + 0.001) * sd) + weightMomentum[toUpdate]
            newWeight[newWeight < 0] = 0
        }

        ## Find the new errors
        means = populationMeans(X = X, fam = fam, alpha = newAlpha, weight = newWeight)[pops, ]
        predicted = predict(X = X, means = means, weight = newWeight, referenceGroups = referenceGroups)
        newError = squaredError(expected = expected, predicted = predicted)
        newError = cbind(apply(newError, 1, sum))
        newError = aggregate(newError, by = list(fam[ , 'group']), mean)
        newError = newError[ , -1]
        newError = mean(newError)
        
        ## Find the jump odds
        odds = c(oldError, newError) - (min(c(oldError, newError)) - (mean(c(oldError, newError)) * 0.0001))
        odds = odds / sum(odds)

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
        
        iteration = iteration + 1
    }
    
    return(minParams)
    
}

rye.optimize = function(X = NULL, fam = NULL,
                    referencePops = NULL, referenceGroups = NULL,
                    alpha = NULL, maxAlpha = 0.1, optimizeAlpha = TRUE,
                    weight = NULL, maxWeight = 2, optimizeWeight = TRUE,
                    iterations = 10, rounds = 10, threads = 20, startSD = 0.25, endSD = 0.025,
                    populationError = FALSE) {

    ## Pull out the reference PCs
    referenceFAM = fam[fam[ , 'population'] %in% referencePops , ]
    referenceX = X[rownames(referenceFAM), ]


    ## Start with the shrinking at 0.05 for all pops by default
    if (is.null(alpha)) {
        alpha = rep(0.05, length(referencePops))
    } 
    names(alpha) = referencePops

    ## Weights
    if (is.null(weight)) {
        weight = 1 / seq(ncol(X))
    }

    round = 1
    allParams = list()
    allErrors = c()
    while (round <= rounds) {

        sd = startSD - (startSD - endSD) * log(round)/log(rounds)
        params = mclapply(seq(threads), function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
#        params = lapply(1, function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
                                                          iterations = iterations,
                                                          alpha = alpha, weight = weight, sd = sd,
                                                          optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight), mc.cores = 32)
#                                                                 optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight))


        errors = unlist(lapply(params, function(i) i[[1]]))
        
        bestError = which.min(errors)
        meanError = mean(errors)
        print(paste('Round', round, '/', rounds, 'Mean MSE:', sprintf("%.8f", meanError), ', Best MSE:', sprintf('%.8f', errors[bestError])))
        
        bestParams = params[[bestError]]
        alpha = bestParams[[2]]
        weight = bestParams[[3]]

        allParams = c(allParams, list(bestParams))
        allErrors = c(allErrors, errors[bestError])
        
        round = round + 1

    }

    lastChange = length(allParams)
    return(allParams[[lastChange]])

        
}

