#Testing Script
require(abc)

datlen = 10000
peturbAmt = 0.05

loc2 = 3
x = c(rnorm(datlen* (1-peturbAmt)), rnorm(datlen*peturbAmt, loc2, 1))

sortx = sort(x)

rg = seq(min(x), max(x), 0.03)

actualP = function(xg){
  return(0.95* dnorm(xg) + 0.05*dnorm(xg, loc2))
} 

nbins = nclass.FD(x)
origHist = hist(x, nbins, freq=FALSE, main = "Contamination Model")
lines(rg, 0.95*dnorm(rg), lw=2, col='red', lt=3)
lines(rg, 0.05*dnorm(rg, loc2), lw=2, col='blue', lt=3)


priors = rnorm(datlen, 0, 3)

sumstats = function(pmean){
  rdata = rnorm(datlen, pmean)
  return(c(median(rdata), IQR(rdata)))
}

sumCalc = sapply(priors,  sumstats)

mstat = c(median(x), IQR(x))
test = abc(target = mstat, param=priors, sumstat = t(sumCalc), tol = .05, method = 'neuralnet')
hist(test$adj.values, "FD", main = "Estimates of Mean", xlab='mu')
sample(test$adj.values, 10)

pnew = sample(test$adj.values, 1)
pnew
xnewSim = rnorm(datlen, pnew)

filteredX = xnewSim[xnewSim >=origHist$breaks[1]]
filteredX = filteredX[filteredX <= tail(origHist$breaks,1)]
b = length(filteredX)
b
newHist = hist(filteredX, origHist$breaks, xlim = range(origHist$breaks), main="Simulated Main Model")
plot(origHist)

empiricalPeturb=294/datlen



h = 2 * IQR(x) * length(x)^(-1/3)

varianceOrig = (origHist$counts + 1.0 * b^2/datlen^2 * newHist$counts)/h

empiricalPeturb = 294
empiricalPeturb= empiricalPeturb/datlen

rem = (origHist$counts - (1-empiricalPeturb)* newHist$counts)
barplot(rem, h, main = "Differenced Histogram", ylab="Count")

remainder = (origHist$counts - (1-empiricalPeturb)* newHist$counts)/sqrt(varianceOrig)
ind = 1:length(remainder)
remainder[is.nan(remainder) | is.infinite(remainder)] =0
plot(origHist$breaks[-1], remainder, 'l', main = "Standardised Remainder")
remsq = remainder^2
plot(origHist$breaks[-1], remsq, 'l', main = "Standardised Square Remainder")

#plot(varianceOrig)

ls = loess(remsq~ind, span=0.2)
lsrem = loess(abs(remainder)~ind, span=0.2)
lines(origHist$breaks[-1], ls$fitted, col = 'red')

sum(abs(rem[ls$fitted>1]))


# Resample from remainder
inds =  seq(66, 90)
rem0 = pmax(0, rem)

newEmp = sum(rem0[inds])

sampleAmt = rem0[inds]*newEmp/sum(rem0[inds])

sampleAmtorig = sampleAmt

IndexinSort = cumsum(origHist$counts)[inds-1]
IndexinSort
upperIndexinSort = cumsum(origHist$counts)[inds]
upperIndexinSort
origHist$counts[inds]

sampleAmtFinal = pmin(sampleAmt, origHist$counts[inds])
sampleAmt = round(sampleAmt)
sampleAmt

len = length

fullSample = c()
for(i in 1:length(IndexinSort)){
  strata = sortx[(IndexinSort[i]+1):upperIndexinSort[i]]
  currSamp = sample(strata, sampleAmt[i])
  fullSample = c(fullSample, currSamp)
  print(c(len(strata), len(currSamp), len(fullSample)))
}

hist(fullSample, "FD", main = "Resampled values")



#abc on fullsample
contamstat = c(median(fullSample))
contPrior = rnorm(1000, median(fullSample))

contsumstats = function(pmean){
  rdata = rnorm(350, pmean)
  return(c(median(rdata)))
}

contsumCalc = sapply(contPrior,  contsumstats)
cont = abc(target = contamstat, param=contPrior, sumstat = contsumCalc, tol = .1, method = 'neuralnet')
hist(cont$adj.values, "FD", main = "Contamination Parameter")


contnew = sample(cont$adj.values, 1)
contnew

contSim = rnorm(sum(sampleAmt), contnew)


contSimhist = hist(contSim, origHist$breaks)
barplot(origHist$counts - contSimhist$counts, h)

newDiff = origHist$counts - contSimhist$counts

estDat = rep(origHist$breaks[-length(origHist$breaks)], pmax(newDiff,0)) + h * runif(sum(pmax(newDiff,0)))

hist(estDat, "FD")
npriors = rnorm(datlen/10, 0, 3)
nsumCalc = sapply(npriors,  median)
test2 = abc(target = median(estDat), param=npriors, sumstat = nsumCalc, tol = .05, method = 'neuralnet')
hist(test2$adj.values, "FD")
