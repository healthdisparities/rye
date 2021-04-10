library(nnls)
library(Hmisc)
library(parallel)
source('rye/rye.R')

options(width = 220)
options(scipen = 999)
options(digits = 4)

## PCA data from PLINK.  Re-do with eigenstrat removing outliers?
fullPCA = read.table('./unphased/1000Genomeshgdpukbbpruned/extractedChrAllPrunedNoSanHWE40.eigenvec', header = TRUE, row.names = NULL)
fullEigenVal = read.table('./unphased/1000Genomeshgdpukbbpruned/extractedChrAllPrunedNoSanHWE40.eigenval', header = FALSE, row.names = NULL)[,1]
rownames(fullPCA) = fullPCA[ , 2]

## Regenerate the FAM from the PCA ipnut
fam = as.matrix(fullPCA[ , c(1, 2)])
colnames(fam) = c('population', 'id')
rownames(fam) = fam[ , 'id']

## Cast PCA to a matrix & scale the PCs
fullPCA = fullPCA[ , 3:ncol(fullPCA)]
fullPCA = as.matrix(fullPCA)
fullPCA = rye.scale(fullPCA)

## Our reference populations
european = c('FIN', 'GBR', 'French', 'French_Basque', 'IBS', 'TSI', 'Tuscan', 'North_Italian')
asian = c('CHB', 'Dai', 'She', 'Tujia', 'Miaozu', 'JPT', 'KHV', 'Cambodians')
amerindian = c('Karitiana', 'Pima', 'Surui')
southAsian = c('STU')
northIndian = c('GIH')
african = c('GWD', 'YRI', 'ESN', 'LWK', 'Biaka_Pygmies', 'Mbuti_Pygmies')
northAfrican = c('Bedouin', 'Druze', 'Mozabite', 'Palestinian')
persian = c('Kalash')
referencePops = c(european, asian, amerindian, southAsian, northIndian, african, northAfrican, persian)

allPops = unique(fam[ , 'population'])

referenceGroups = c('NorthernEuropean', rep('WesternEuropean', 2), rep('Iberian', 2), rep('Italian', 3),
                    rep('Asian', length(asian)),
                    rep('Amerindian', length(amerindian)),
                    rep('SouthAsian', 1),
                    rep('NorthIndian', 1),
                    rep('Senegambian', 1), rep('Nigerian', 2), 'EastAfrican', rep('Pygmies', 2),
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

## Populations for comparing to RFMix
rfMixEuropean = c(europeanGroup, northAfricanGroup, southAsianGroup, northIndianGroup)
rfMixAsian = c(asianGroup, amerindianGroup)
rfMixAfrican = africanGroup

## RFMix results for comparing admixed individuals
rfMixEstimates = read.table('rfMixFractions.3.Q', row.names = 1, header = TRUE)
colnames(rfMixEstimates)[3] = 'Asian'

allAdmixed = c('ASW', 'ACB', 'MXL', 'CLM', 'PUR', 'PEL')
common = rownames(rfMixEstimates)[rownames(rfMixEstimates) %in% rownames(fullPCA)]
admixedInd = grep(paste(allAdmixed, collapse = '|'), common, val = TRUE)

## Weight the PCs by their eigenvalues
weight = fullEigenVal / max(fullEigenVal)

pcs = 20

## Using each region as a population, e.g., combine British and French to WesternEuropean
regionFAM = fam
regionFAM[fam[,1] %in% names(referenceGroups),1] = referenceGroups[fam[fam[,1] %in% names(referenceGroups), 1]]
referenceGroups = unique(referenceGroups)
names(referenceGroups) = referenceGroups
referencePops = referenceGroups
    
scaledWeight = weight[seq(pcs)]
unifAlpha = rep(0.001, length(referencePops))
names(unifAlpha) = referencePops

optParams = rye.optimize(X = fullPCA[,seq(pcs)], fam = regionFAM,
                         referencePops = referencePops, #referenceGroups = referenceGroups,
                         startSD = 0.01, endSD = 0.005,
                         threads = 25, iterations = 100, rounds = 200, weight = scaledWeight,
                         alpha = unifAlpha, optimizeWeight = TRUE, optimizeAlpha = TRUE)
optWeight = optParams[[3]]
optMeans = optParams[[4]]
    
optEstimates = rye.predict(X = fullPCA[,seq(pcs)], means = optMeans, weight = optWeight)

## Toss out anything under 1%
#optEstimates[optEstimates < 0.01] = 0
optEstimates = t(apply(optEstimates, 1, function(i) i /sum(i)))

## Find the mean of each popultion
optEstimateMeans = do.call(cbind, lapply(allPops, function(i) cbind(apply(t(optEstimates[fam[,1] == i, ,drop = FALSE]), 1, mean))))
colnames(optEstimateMeans) = allPops
optEstimatesAgg = cbind(apply(optEstimates[ , europeanGroup, drop = FALSE], 1, sum), 
                             apply(optEstimates[ , asianGroup, drop = FALSE], 1, sum),
                             apply(optEstimates[ , amerindianGroup, drop = FALSE], 1, sum),
                             apply(optEstimates[ , southAsianGroup, drop = FALSE], 1, sum),
                             apply(optEstimates[ , africanGroup, drop = FALSE], 1, sum),
                             apply(optEstimates[ , northAfricanGroup, drop = FALSE], 1, sum),
                             apply(optEstimates[ , persianGroup, drop = FALSE], 1, sum))
colnames(optEstimatesAgg) = c('European', 'Asian', 'Amerindian', 'SouthAsian', 'African', 'NorthAfrican', 'Persian')

## Combine Asian and NA
optAsianAmerindian = optEstimates
optAsianAmerindian[ , 'Asian'] = apply(optAsianAmerindian[ , rfMixAsian, drop = FALSE], 1, sum)
optAsianAmerindian= cbind(optAsianAmerindian, apply(optAsianAmerindian[ , rfMixEuropean, drop = FALSE], 1, sum))
colnames(optAsianAmerindian)[ncol(optAsianAmerindian)] = 'European'
optAsianAmerindian = cbind(optAsianAmerindian, apply(optAsianAmerindian[ , rfMixAfrican, drop = FALSE], 1, sum))
colnames(optAsianAmerindian)[ncol(optAsianAmerindian)] = 'African'
optError = as.matrix(rye.absoluteError(expected = rfMixEstimates[admixedInd, ], predicted = optAsianAmerindian[admixedInd, colnames(rfMixEstimates)]))
optError = mean(optError)
print(paste('Optimized Reference error:', paste0(sprintf('%.6f', optParams[[1]]), ','), 'Admixed Continental error:', sprintf('%.6f', optError))) 

## Correlate with RFMix values
optRValues = matrix(0, nrow = length(allAdmixed), ncol = length(c('European', 'African', 'Amerindian')), dimnames = list(allAdmixed, c('European', 'African', 'Amerindian')))
for (pop in allAdmixed) {
    optRValues[pop, 'European'] = rcorr(apply(optEstimates[grep(pop, rownames(optEstimates), val = TRUE), rfMixEuropean, drop = FALSE], 1, sum),
                                                              rfMixEstimates[grep(pop, rownames(optEstimates), val = TRUE), 'European'])$r[1,2]
    optRValues[pop, 'African'] = rcorr(apply(optEstimates[grep(pop, rownames(optEstimates), val = TRUE), rfMixAfrican, drop = FALSE], 1, sum),
                                                             rfMixEstimates[grep(pop, rownames(optEstimates), val = TRUE), 'African'])$r[1,2]
    optRValues[pop, 'Amerindian'] = rcorr(apply(optEstimates[grep(pop, rownames(optEstimates), val = TRUE), rfMixAsian, drop = FALSE], 1, sum),
                                                                rfMixEstimates[grep(pop, rownames(optEstimates), val = TRUE), 'Asian'])$r[1,2]
}
print(optRValues)

## Save stupid results of stupid people
write.table(x = optEstimatesAgg, file = paste0('unphased/1000Genomeshgdpukbbpruned/continentalEstimatesOpt', pcs, '.7.Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
write.table(x = optEstimates, file = paste0('unphased/1000Genomeshgdpukbbpruned/subcontinentalEstimatesOpt', pcs, '.', ncol(optEstimates), '.Q'),
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
write.table(x = fam[rownames(optEstimatesAgg), ], file = paste0('unphased/1000Genomeshgdpukbbpruned/continentalEstimatesOpt', pcs, '.fam'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
    

## Without optimization for comparison
unifMeans = rye.populationMeans(X = fullPCA[regionFAM[ , 'population'] %in% referencePops, seq(pcs)],
                                fam = regionFAM[regionFAM[ , 'population'] %in% referencePops, ],
                                alpha = unifAlpha, weight = scaledWeight[seq(pcs)])

unifEstimates = rye.predict(X = fullPCA[,seq(pcs)], means = unifMeans, weight = scaledWeight[seq(pcs)])

## Roll up the means for each population
unifEstimateMeans = do.call(cbind, lapply(allPops, function(i) cbind(apply(t(unifEstimates[fam[,1] == i, ,drop = FALSE]), 1, mean))))
colnames(unifEstimateMeans) = allPops

## Aggregate the continental ancestries
estimates = unifEstimates
unifEstimatesAgg = cbind(apply(estimates[ , europeanGroup, drop = FALSE], 1, sum), 
                         apply(estimates[ , asianGroup, drop = FALSE], 1, sum),
                         apply(estimates[ , amerindianGroup, drop = FALSE], 1, sum),
                         apply(estimates[ , southAsianGroup, drop = FALSE], 1, sum),
                         apply(estimates[ , africanGroup, drop = FALSE], 1, sum),
                         apply(estimates[ , northAfricanGroup, drop = FALSE], 1, sum),
                         apply(estimates[ , persianGroup, drop = FALSE], 1, sum))
colnames(unifEstimatesAgg) = c('European', 'Asian', 'Amerindian', 'SouthAsian', 'African', 'NorthAfrican', 'Persian')

## Get the merged resultsfor RFMix
unifAsianAmerindian = unifEstimates
unifAsianAmerindian[ , 'Asian'] = unifAsianAmerindian[ , 'Asian'] + unifAsianAmerindian[ , 'Amerindian']
unifAsianAmerindian = cbind(unifAsianAmerindian, apply(unifAsianAmerindian[ , rfMixEuropean, drop = FALSE], 1, sum))
colnames(unifAsianAmerindian)[ncol(unifAsianAmerindian)] = 'European'
unifAsianAmerindian = cbind(unifAsianAmerindian, apply(unifAsianAmerindian[ , rfMixAfrican, drop = FALSE], 1, sum))
colnames(unifAsianAmerindian)[ncol(unifAsianAmerindian)] = 'African'
unifError = as.matrix(rye.absoluteError(expected = rfMixEstimates[admixedInd, ], predicted = unifAsianAmerindian[admixedInd, colnames(rfMixEstimates)]))
unifError = mean(unifError)
print(paste('Uniform alpha, scaled weight Admixed Continetnal error:', sprintf('%.6fe', unifError)))

## AFR/EUR/EAS correlation between Rye and RFMix
unifRValues = matrix(0, nrow = length(allAdmixed), ncol = length(c('European', 'African', 'Amerindian')), dimnames = list(allAdmixed, c('European', 'African', 'Amerindian')))
for (pop in allAdmixed) {
    unifRValues[pop, 'European'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixEuropean, drop = FALSE], 1, sum),
                                         rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'European'])$r[1,2]
    unifRValues[pop, 'African'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAfrican, drop = FALSE], 1, sum),
                                        rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'African'])$r[1,2]
    unifRValues[pop, 'Amerindian'] = rcorr(apply(estimates[grep(pop, rownames(estimates), val = TRUE), rfMixAsian, drop = FALSE], 1, sum),
                                           rfMixEstimates[grep(pop, rownames(estimates), val = TRUE), 'Asian'])$r[1,2]
        
}
print(unifRValues)

## Same stupid results for dumb idiots
write.table(x = unifEstimates, file = paste0('unphased/1000Genomeshgdpukbbpruned/subcontinentalEstimatesUnif', pcs, '.', ncol(unifEstimates), '.Q'), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')

## Plot the Rye-RFMix correlation for the admixed 1KGP individuals
colors = c('#DFB028', '#400B8E',  '#C70F00', '#BB33C7')
names(colors) = c('European', 'African', 'Asian', 'SouthAsian')

aggregateEstimates = list(optAsianAmerindian, unifAsianAmerindian)
names(aggregateEstimates) = c('Optimized', 'Uniform')

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

        error = mean(rye.absoluteError(expected = rfMixEstimates[admixedInd, i, drop = FALSE], predicted = estimatesAgg[admixedInd, i, drop = FALSE])[,1])
        print(error)
 
        
	lines(x = c(0, 1), y = c(0, 1), lwd = 1/2, lty = 2)
	points(x = rfMixEstimates[common, i], y = estimatesAgg[common, i], pch = 19, cex = 1/3, col = colors[i])
        
        text(x = .1, y = .8, label = paste('r = ', sprintf("%.3f", rValues[i])), adj = c(0, 0))
    	text(x = .1, y = .7, label = paste('Mean Error = ', sprintf("%.3f", error)), adj = c(0, 0))
    }
}
dev.off()


