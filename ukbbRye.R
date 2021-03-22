library(nnls)
library(Hmisc)
library(parallel)
source('rye.R')

options(width = 220)
options(scipen = 999)
options(digits = 4)

fullPCA = read.table('./unphased/1000Genomeshgdpukbbpruned/extractedChrAllPruned.eigenvec', header = FALSE, row.names = NULL)
rownames(fullPCA) = fullPCA[ , 2]

## Regenerate the FAM from the PCA inut
fam = as.matrix(fullPCA[ , c(1, 2)])
colnames(fam) = c('population', 'id')
rownames(fam) = fam[ , 'id']

## Cast PCA to a matrix & scale
fullPCA = fullPCA[ , 3:ncol(fullPCA)]
fullPCA = as.matrix(fullPCA)
fullPCA = scale(fullPCA)

## Our reference populations
european = c('FIN', 'Orcadian', 'GBR', 'French', 'French_Basque', 'IBS', 'TSI', 'Tuscan', 'North_Italian')
asian = c('CHB', 'Dai', 'She', 'Tujia', 'Miaozu', 'JPT', 'KHV', 'Cambodians')
amerindian = c('Karitiana', 'Pima', 'Surui')
southAsian = c('ITU', 'STU')
northIndian = c('GIH')
african = c('GWD', 'MSL', 'YRI', 'ESN', 'LWK', 'Biaka_Pygmies', 'Mbuti_Pygmies')
northAfrican = c('Bedouin', 'Druze', 'Mozabite', 'Palestinian')
persian = c('Kalash')
referencePops = c(european, asian, amerindian, southAsian, northIndian, african, northAfrican, persian)

referenceGroups = referencePops
names(referenceGroups) = referencePops
europeanGroup = european
asianGroup = asian
amerindianGroup = amerindian
southAsianGroup = southAsian
northIndianGroup = northIndian
africanGroup = african
northAfricanGroup = northAfrican
persianGroup = persian

referenceGroups = c('NorthernEuropean', rep('WesternEuropean', 3), rep('Iberian', 2), rep('Italian', 3),
                    rep('Asian', length(asian)),
                    rep('Amerindian', length(amerindian)),
                    rep('SouthAsian', 2),
                    rep('NorthIndian', 1),
                    rep('Senegambian', 2), rep('Nigerian', 2), 'EastAfrican', rep('Pygmies', 2),
                    rep('NorthAfrican', 4),
                    rep('Persian', 1))
names(referenceGroups) = referencePops
europeanGroup = c('NorthernEuropean', 'WesternEuropean', 'Iberian', 'Italian')
asianGroup = c('Asian')
amerindianGroup = c('Amerindian')
southAsianGroup = c('SouthAsian')

northIndianGroup = c('NorthIndian')
africanGroup = c('Senegambian', 'Nigerian', 'EastAfrican', 'Pygmies')
northAfricanGroup = c('NorthAfrican')
persianGroup = c('Persian')


rfMixEuropean = c(europeanGroup, northAfricanGroup)
rfMixAsian = c(asianGroup, amerindianGroup)
rfMixAfrican = africanGroup

## Compare each to RFMix
rfMixEstimates = read.table('rfMixFractions.3.Q', row.names = 1, header = TRUE)
colnames(rfMixEstimates)[3] = 'Asian'

common = rownames(rfMixEstimates)[rownames(rfMixEstimates) %in% rownames(fullPCA)]
admixedInd = grep(paste(allAdmixed, collapse = '|'), common, val = TRUE)

## Across PCs

for (pcs in seq(from = 25, to = 25, by = 1)) {

    unifAlpha = rep(0.005, length(referencePops))
    names(unifAlpha) = referencePops

    optAlphaParams = rye.optimize(X = fullPCA[,seq(pcs)], fam = fam,
                                       referencePops = referencePops, referenceGroups = referenceGroups,
                                       startSD = 0.25, endSD = 0.25,
                                       threads = 100, iterations = 200, rounds = 100, weight = scaledWeight[seq(pcs)],
                                       alpha = unifAlpha, optimizeWeight = TRUE, optimizeAlpha = TRUE, populationError = FALSE)
    optAlphaWeight = optAlphaParams[[3]]
    optAlphaMeans = optAlphaParams[[4]]
    optAlphaEstimates = predict(X = fullPCA[,seq(pcs)], means = optAlphaMeans, weight = optAlphaWeight, referenceGroups = referenceGroups)
    write.table(x = optAlphaEstimates, file = paste0('unphased/1000Genomeshgdpukbbpruned/populationEstimatesAlpha', pcs, '.7.Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    optAlphaEstimates[optAlphaEstimates < 0.005] = 0
    optAlphaEstimates = t(apply(optAlphaEstimates, 1, function(i) i /sum(i)))
    optAlphaEstimateMeans = do.call(cbind, lapply(allPops, function(i) cbind(apply(t(optAlphaEstimates[fam[,1] == i, ,drop = FALSE]), 1, mean))))
    colnames(optAlphaEstimateMeans) = allPops
    estimates = optAlphaEstimates
    optAlphaEstimatesAgg = cbind(apply(estimates[ , europeanGroup, drop = FALSE], 1, sum), 
                        apply(estimates[ , asianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , amerindianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , southAsianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , africanGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , northAfricanGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , persianGroup, drop = FALSE], 1, sum))
    colnames(optAlphaEstimatesAgg) = c('European', 'Asian', 'Amerindian', 'SouthAsian', 'African', 'NorthAfrican', 'Persian')
    optAlphaAsianAmerindian = optAlphaEstimatesAgg
    optAlphaAsianAmerindian[ , 'Asian'] = optAlphaAsianAmerindian[ , 'Asian'] + optAlphaAsianAmerindian[ , 'Amerindian']
    optAlphaAsianAmerindian[ , 'European'] = optAlphaAsianAmerindian[ , 'European'] + optAlphaAsianAmerindian[ , 'NorthAfrican']
    optAlphaSSE = meanSquaredError(expected = rfMixEstimates[admixedInd, ], predicted = optAlphaAsianAmerindian[admixedInd, colnames(rfMixEstimates)])
    print(paste('Optimized alpha, scaled weight,', pcs, 'PCs, Reference MSE:', paste0(sprintf('%.8f', optAlphaParams[[1]]), ','), 'Admixed Continental MSE:', sprintf('%.8f', optAlphaSSE))) 

    optAlphaRValues = matrix(0, nrow = length(allAdmixed), ncol = length(c('European', 'African', 'Amerindian')), dimnames = list(allAdmixed, c('European', 'African', 'Amerindian')))
    for (pop in allAdmixed) {
        optAlphaRValues[pop, 'European'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixEuropean, drop = FALSE], 1, sum),
                                                              rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'European'])$r[1,2]
        optAlphaRValues[pop, 'African'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAfrican, drop = FALSE], 1, sum),
                                                             rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'African'])$r[1,2]
        optAlphaRValues[pop, 'Amerindian'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAsian, drop = FALSE], 1, sum),
                                                                rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'Asian'])$r[1,2]
        
    }
    print(optAlphaRValues)
    write.table(x = optAlphaEstimatesAgg, file = paste0('unphased/1000Genomeshgdpukbbpruned/continentalEstimatesAlpha', pcs, '.7Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    write.table(x = optAlphaEstimates, file = paste0('unphased/1000Genomeshgdpukbbpruned/subcontinentalEstimatesAlpha', pcs, '.', ncol(optAlphaEstimates), '.Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    write.table(x = fam[rownames(optAlphaEstimatesAgg), ], file = paste0('unphased/1000Genomeshgdpukbbpruned/continentalEstimatesAlpha', pcs, '.fam'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    

    unifAlphaScaledWeightMeans = populationMeans(X = fullPCA[fam[fam[ , 'population'] %in% referencePops, 'id'], seq(pcs)],
                           fam = fam[fam[fam[ , 'population'] %in% referencePops, 'id'], ],
                           alpha = unifAlpha, weight = scaledWeight[seq(pcs)])
    unifAlphaScaledWeightEstimates = predict(X = fullPCA[,seq(pcs)], means = unifAlphaScaledWeightMeans, weight = scaledWeight[seq(pcs)], referenceGroups = referenceGroups)
    unifAlphaScaledWeightEstimates[unifAlphaScaledWeightEstimates < 0.005] = 0
    unifAlphaScaledWeightEstimates = t(apply(unifAlphaScaledWeightEstimates, 1, function(i) i /sum(i)))
    unifAlphaScaledWeightEstimateMeans = do.call(cbind, lapply(allPops, function(i) cbind(apply(t(unifAlphaScaledWeightEstimates[fam[,1] == i, ,drop = FALSE]), 1, mean))))
    colnames(unifAlphaScaledWeightEstimateMeans) = allPops
    estimates = unifAlphaScaledWeightEstimates
    unifAlphaScaledWeightEstimatesAgg = cbind(apply(estimates[ , europeanGroup, drop = FALSE], 1, sum), 
                        apply(estimates[ , asianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , amerindianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , southAsianGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , africanGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , northAfricanGroup, drop = FALSE], 1, sum),
                        apply(estimates[ , persianGroup, drop = FALSE], 1, sum))
    colnames(unifAlphaScaledWeightEstimatesAgg) = c('European', 'Asian', 'Amerindian', 'SouthAsian', 'African', 'NorthAfrican', 'Persian')
    unifAlphaScaledWeightAsianAmerindian = unifAlphaScaledWeightEstimatesAgg
    unifAlphaScaledWeightAsianAmerindian[ , 'Asian'] = unifAlphaScaledWeightAsianAmerindian[ , 'Asian'] + unifAlphaScaledWeightAsianAmerindian[ , 'Amerindian']
    unifAlphaScaledWeightAsianAmerindian[ , 'European'] = unifAlphaScaledWeightAsianAmerindian[ , 'European'] + unifAlphaScaledWeightAsianAmerindian[ , 'NorthAfrican']
    unifAlphaScaledWeightSSE = meanSquaredError(expected = rfMixEstimates[admixedInd, ], predicted = unifAlphaScaledWeightAsianAmerindian[admixedInd, colnames(rfMixEstimates)])
    print(paste('Uniform alpha, scaled weigh,', pcs, 'PCs, Continetnal MSE:', sprintf('%.8f', unifAlphaScaledWeightSSE)))
    unifAlphaScaledWeightRValues = matrix(0, nrow = length(allAdmixed), ncol = length(c('European', 'African', 'Amerindian')), dimnames = list(allAdmixed, c('European', 'African', 'Amerindian')))
    for (pop in allAdmixed) {
        unifAlphaScaledWeightRValues[pop, 'European'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixEuropean, drop = FALSE], 1, sum),
                                                              rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'European'])$r[1,2]
        unifAlphaScaledWeightRValues[pop, 'African'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAfrican, drop = FALSE], 1, sum),
                                                             rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'African'])$r[1,2]
        unifAlphaScaledWeightRValues[pop, 'Amerindian'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAsian, drop = FALSE], 1, sum),
                                                                rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'Asian'])$r[1,2]
        
    }
    write.table(x = unifAlphaScaledWeightEstimates, file = paste0('unphased/1000Genomeshgdpukbbpruned/subcontinentalEstimatesUnif', pcs, '.', ncol(unifAlphaScaledWeightEstimates), '.Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    print(unifAlphaScaledWeightRValues)

}

barf
colors = c('#DFB028', '#400B8E',  '#C70F00', '#BB33C7')
names(colors) = c('European', 'African', 'Asian', 'SouthAsian')

pdf(file = 'plots/pcaVsRFMix.pdf', width = 3.5*3, height = 3.5*length(aggregateEstimates), useDingbats = FALSE)
par(mfrow = c(length(aggregateEstimates), 3))
par(mar = c(4,4,1,1))

for (pop in names(aggregateEstimates)) {

    estimatesAgg = aggregateEstimates[[pop]]
    rValues = unlist(lapply(colnames(rfMixEstimates), function(i)rcorr(estimatesAgg[common, i], rfMixEstimates[common, i], 'pearson')$r[1,2]))
    names(rValues) = colnames(rfMixEstimates)
   
    for (i in colnames(rfMixEstimates)) {
        plot(type = 'n', bty = 'n', las = 1,
             xlab = paste('RFMix', i), ylab = paste('PCA', i),
             x = c(0, 1), y = c(0, 1))

        error = squaredError(expected = rfMixEstimates[common, i, drop = FALSE], predicted = estimatesAgg[common, i, drop = FALSE])
        print(error)
 
        
	lines(x = c(0, 1), y = c(0, 1), lwd = 1/2, lty = 2)
	points(x = rfMixEstimates[common, i], y = estimatesAgg[common, i], pch = 19, cex = 1/3, col = colors[i])
        
        text(x = .1, y = .8, label = paste('r = ', sprintf("%.3f", rValues[i])), adj = c(0, 0))
    	text(x = .1, y = .7, label = paste('SSE = ', sprintf("%.3f", error)), adj = c(0, 0))
    }
}
dev.off()


