#!/usr/bin/env Rscript
script_author = "Andrew Conley, Lavanya Rishishwar"
script_copyright = "Copyright 2021, Andrew Conley, Lavanya Rishishwar"
script_credits = c("Andrew Conely", "Lavanya Rishishwar", "Maria Ahmad", "Shivam Sharma", "Emily Norris")
script_license = "GPL"
script_version = "0.1"
script_maintainer = "Andrew Conley, Lavanya Rishishwar"
script_email = "aconley@ihrc.com; lrishishwar@ihrc.com"
script_status = "Development"
script_title = "rye.R"


################################################################################
### Load libraries
requiredPackages = c('nnls','Hmisc','parallel', 'optparse', 'crayon')
for(p in requiredPackages){
  if(!suppressMessages(require(p,character.only = TRUE, quietly = T))){
    stop(paste0("Library ", p, " is required, I can't seem to find it."))
  } 
}

options(width = 220)
options(scipen = 999)
options(digits = 4)

################################################################################
### Function definition
printDims = function(X, msg){cat(paste(msg, ':', paste(dim(X), collapse = 'x'), "\n"))}

printError = function(msg){cat(red(paste(msg, collapse = " "), "\n"))}
printWarn = function(msg){cat(yellow(paste(msg, collapse = " "), "\n"))}
logmsg = function(msg){cat(green(paste(format(Sys.time(), "[ %b %d %Y - %X ]"), msg, collapse = " "), "\n"))}
progressmsg = function(msg){cat(magenta(paste(msg, collapse = " "), "\n"))}
pretty_time = function(time){
  out_string = ""
  if(time > 60*60*24){
    days = round(time/(60*60*24))
    time = time %% (60*60*24)
    out_string = paste0(out_string, days, " days, ")
  }
  if(time > 60*60){
    hours = round(time/(60*60))
    time = time %% (60*60)
    out_string = paste0(out_string, hours, " hours, ")
  }
  if(time > 60){
    mins = round(time/60)
    time = time %% 60
    out_string = paste0(out_string, mins, " mins, ")
  }
  out_string = paste0(out_string, round(time, 2), " seconds")
  return(out_string)
}

rye.scale = function(X = NULL) {
  return(apply(X, 2, function(i){i = i - min(i); i / max(i)}))
}

rye.populationMeans = function(X = NULL, fam = NULL, alpha = NULL, weight = NULL, fn = median, referenceGroups = NULL) {
  
  ## Find the mean of each reference population
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
                        weight = NULL, optimizeWeight = TRUE, attempts = 4,
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
      params = mclapply(seq(attempts), function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
                                                            iterations = iterations,
                                                            alpha = alpha, weight = weight, sd = sd,
                                                            optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight), mc.cores = threads)
    } else {
      params = lapply(seq(attempts), function(i) rye.gibbs(X = referenceX, fam = referenceFAM, referenceGroups = referenceGroups,
                                                          iterations = iterations,
                                                          alpha = alpha, weight = weight, sd = sd,
                                                          optimizeAlpha = optimizeAlpha, optimizeWeight = optimizeWeight))
    } 
    
    
    errors = unlist(lapply(params, function(i) i[[1]]))
    
    bestError = which.min(errors)
    meanError = mean(errors)
    progressmsg(paste0('Round ', round, '/', rounds, ' Mean error: ', sprintf("%.6f", meanError),
               ', Best error: ', sprintf('%.6f', errors[bestError])))
    
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

rye = function(eigenvec_file = NULL, eigenval_file = NULL,
               pop2group_file = NULL, output_file = NULL,
               threads = 4, pcs = 20, optim_rounds = 200,
               optim_iter = 100, attempts=4){
  ## Perform core operation
  #TODO: Change file reading method to data.table
  logmsg("Reading in Eigenvector file")
  fullPCA = read.table(eigenvec_file, header = FALSE, row.names = NULL)
  rownames(fullPCA) = fullPCA[ , 2]
  logmsg("Reading in Eigenvalue file")
  fullEigenVal = read.table(eigenval_file, header = FALSE, row.names = NULL)[,1]
  logmsg("Reading in pop2group file")
  pop2group = read.table(pop2group_file, header = T, stringsAsFactors = F)
  referenceGroups = pop2group$Group
  names(referenceGroups) = pop2group$Pop
  
  ## Regenerate the FAM from the PCA input
  logmsg("Creating individual mapping")
  fam = as.matrix(fullPCA[ , c(1, 2)])
  colnames(fam) = c('population', 'id')
  rownames(fam) = fam[ , 'id']
  allPops = unique(fam[ , 'population'])
  
  ## Cast PCA to a matrix & scale the PCs
  logmsg("Scaling PCs")
  fullPCA = fullPCA[ , 3:ncol(fullPCA)]
  fullPCA = as.matrix(fullPCA)
  fullPCA = rye.scale(fullPCA)
  
  ## Weight the PCs by their eigenvalues
  logmsg("Weighting PCs")
  weight = fullEigenVal / max(fullEigenVal)
  
  ## Using each region as a population, e.g., combine British and French to WesternEuropean
  logmsg("Aggregating individuals to population groups")
  regionFAM = fam
  regionFAM[fam[,1] %in% names(referenceGroups),1] = referenceGroups[fam[fam[,1] %in% names(referenceGroups), 1]]
  referenceGroups = unique(referenceGroups)
  names(referenceGroups) = referenceGroups
  referencePops = referenceGroups
  
  ## Optimize estimates using NNLS
  logmsg("Optimizing estimates using NNLS")
  scaledWeight = weight[seq(pcs)]
  unifAlpha = rep(0.001, length(referencePops))
  names(unifAlpha) = referencePops
  optParams = rye.optimize(X = fullPCA[,seq(pcs)], fam = regionFAM,
                           referencePops = referencePops, referenceGroups = referenceGroups,
                           startSD = 0.01, endSD = 0.005,
                           threads = threads, iterations = optim_iter,
                           rounds = optim_rounds, attempts=attempts,
                           weight = scaledWeight, alpha = unifAlpha, 
                           optimizeWeight = TRUE, optimizeAlpha = TRUE)
  optWeight = optParams[[3]]
  optMeans = optParams[[4]]
  
  ## Calculate ancestry estimates
  logmsg("Calculate per-individual ancestry estimates")
  optEstimates = rye.predict(X = fullPCA[,seq(pcs)], means = optMeans, weight = optWeight)
  optEstimates = t(apply(optEstimates, 1, function(i) i /sum(i)))
  
  ## Find the mean of each population
  # logmsg("Calculate per-population mean ancestry estimates")
  # optEstimateMeans = do.call(cbind, lapply(allPops, function(i) cbind(apply(t(optEstimates[fam[,1] == i, ,drop = FALSE]), 1, mean))))
  # colnames(optEstimateMeans) = allPops
  optEstimatesAgg = NULL
  for(group in referenceGroups){
    optEstimatesAgg = cbind(optEstimatesAgg, apply(optEstimates[ , group, drop = FALSE], 1, sum))
  }
  colnames(optEstimatesAgg) = as.character(referenceGroups)
  
  ## Create output files
  logmsg("Create output files")
  write.table(x = optEstimatesAgg, 
              file = paste0(output_file, '-', pcs, '.', length(referenceGroups),'.Q'), 
              col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  write.table(x = optEstimates, 
              file = paste0(output_file, '-', pcs, '.', ncol(optEstimates), '.Q'),
              col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  write.table(x = fam[rownames(optEstimatesAgg), ], 
              file = paste0(output_file, '-', pcs, '.fam'), col.names = TRUE, 
              row.names = TRUE, quote = FALSE, sep = '\t')
  
}


validate_arguments <- function(opt){
  ## Verify the arguments
  #TODO: Move argument validation to its own function
  argumentsGood = TRUE
  if (is.null(opt$eigenvec)) {
    argumentsGood = FALSE
    printError('Eigenvector file not given (--eigenvec)')
  } else if (!file.exists(opt$eigenvec)) {
    argumentsGood = FALSE
    printError(paste('Eigenvector file (--eigenvec=', opt$eigenvec, ') not found'))
  }
  if (is.null(opt$eigenval)) {
    argumentsGood = FALSE
    printError('Eigenvalue file not given (--eigenval)')
  } else if (!file.exists(opt$eigenval)) {
    argumentsGood = FALSE
    printError(paste('Eigenvalue file (--eigenval=', opt$eigenval, ') not found'))
  }
  if (is.null(opt$pop2group)) {
    argumentsGood = FALSE
    printError('Population-to-group mapping file not given (--pop2group)')
  } else if (!file.exists(opt$pop2group)) {
    argumentsGood = FALSE
    printError(paste('Population-to-group mapping file (--pop2group=', opt$pop2group, ') not found'))
  }
  if (is.null(opt$output)) {
    argumentsGood = FALSE
    printError('Output prefix not given (--output)')
  }
  
  #TODO: Implement check for file dimensions
  #TODO: Ensure number of threads don't exceed machine capacity
  
  if (!argumentsGood){
    printError('Incomplete/incorrect arguments were observed, cannot continue.')
    printError(c("Run", script_title, "-h for usage information"))
    # stop("Exiting...") # I don't like the error message
    q(save = "no", status = 1)
  }
}

################################################################################

optionList = list(
  make_option('--eigenval', type = 'character', default = NULL,
              help = 'Eigenvalue file [REQUIRED]', metavar = '<EVAL_FILE>'),
  make_option('--eigenvec', type = 'character', default = NULL,
              help = 'Eigenvector file [REQUIRED]', metavar = '<EVEC_FILE>'),
  make_option('--pop2group', type = 'character', default = NULL,
              help = 'Population-to-group mapping file [REQUIRED]', metavar = '<P2G_FILE>'),
  make_option('--output', type = 'character', default = "output",
              help = 'Output prefix (Default = output)', metavar = '<OUTPUT_PREFIX>'),
  make_option('--threads', type = 'numeric', default = 4,
              help = 'Number of threads to use (Default = 4)', metavar = '<THREADS>'),
  make_option('--pcs', type = 'numeric', default = 20,
              help = 'Number of PCs to use (Default = 20)', metavar = '<#PCs>'),
  make_option('--rounds', type = 'numeric', default = 200,
              help = 'Number of rounds to use for optimization (higher number = more accurate but slower; Default=200)',
              metavar = '<optim-rounds>'),
  make_option('--iter', type = 'numeric', default = 100,
              help = 'Number of iterations to use for optimization (higher number = more accurate but slower; Default=100)',
              metavar = '<optim-iters>'),
  make_option('--attempts', type = 'numeric', default = 4,
              help = 'Number of attempts to find the optimum values (Default = 4)', metavar = '<ATTEMPTS>')
)

optParser = OptionParser(option_list = optionList)
opt = parse_args(optParser)
# Debug only
# opt = parse_args(optParser, args = c("--eigenvec=extractedChrAllPrunedNoSan.25.eigenvec.gz",
#                                      "--eigenval=extractedChrAllPrunedNoSan.25.eigenval",
#                                      "--pop2group=pop2group.txt"))
# print(opt)
start_time <- proc.time()
logmsg("Parsing user supplied arguments...")
validate_arguments(opt)
logmsg("Arguments passed validation")
logmsg(paste0("Running core rye with ", opt$threads, " threads"))
rye(eigenvec_file = opt$eigenvec,
    eigenval_file = opt$eigenval,
    pop2group_file = opt$pop2group,
    output_file = opt$output,
    threads = opt$threads,
    attempts = opt$attempts,
    pcs = opt$pcs,
    optim_rounds = opt$rounds,
    optim_iter = opt$iter)
logmsg("Process completed")
end_time <- proc.time() - start_time
logmsg(end_time)
logmsg(paste0("The process took ", pretty_time(end_time["user.self"])))
