graphics.off()
rm(list=ls(all=TRUE))
require(rjags)
require(runjags)
load.module('dic')
load.module('lba')

# -----------------------------------------------------------------------------
# MODEL SPECIFICATION
# -----------------------------------------------------------------------------

modelstring = "
data {
    for (c in 1:nCondition) {
        for (s in 1:nSubject) {
            inCondition[c, s] <- 1 * (condition[s] == c)
        }
        nInCondition[c] <- sum(inCondition[c,])
    }
}
model {
    for (i in 1:nData) {
        A[i] <- sim[subject[i]] * sum(M[subject[i], 1:setSize[i]]) + (1 - sim[subject[i]]) * M[subject[i], lag[i]] + background[subject[i]]
        
        drift[i] <- A[i] / (A[i] + u[subject[i]] + v[subject[i]] * setSize[i])
        
        rt[i] ~ dlba(drift[i], 1 - drift[i], boundsep[subject[i]] * (1 - bias[subject[i]]) + startvar[subject[i]], boundsep[subject[i]] * bias[subject[i]] + startvar[subject[i]], startvar[subject[i]], startvar[subject[i]], eta[subject[i]], eta[subject[i]], ter[subject[i]], pContam[subject[i]], 0.5, 0.2, 5)
        
        # Predictive 
        predRTRaw[i] ~ dlba(drift[i], 1 - drift[i], boundsep[subject[i]] * (1 - bias[subject[i]]) + startvar[subject[i]], boundsep[subject[i]] * bias[subject[i]] + startvar[subject[i]], startvar[subject[i]], startvar[subject[i]], eta[subject[i]], eta[subject[i]], ter[subject[i]], pContam[subject[i]], 0.5, 0.2, 5)
        predRT[i] <- abs(predRTRaw[i]) * (((predRTRaw[i] > 0) * (corrResp[i] == 1)) + ((predRTRaw[i] < 0) * (corrResp[i] == 0))) - abs(predRTRaw[i]) * (((predRTRaw[i] > 0) * (corrResp[i] == 0)) + ((predRTRaw[i] < 0) * (corrResp[i] == 1)))
    }
    
    # Summed strengths
    
    for (s in 1:nSubject) {
        for (l in 1:maxSetSize) {
            M[s, l] <- pow(l, -beta[s])
        }
        M[s, maxSetSize + 1] <- 0
    }
    
    # Boundary separation
    
    logBoundsepMean ~ dnorm(0, 0.001)
    
    logBoundsepCondPrec ~ dgamma(1.001, 0.001)
    logBoundsepSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logBoundsepCond[c] ~ dnorm(0, logBoundsepCondPrec)
    }
    
    for (s in 1:nSubject) {
        logBoundsepSubj[s] ~ dnorm(0, logBoundsepSubjPrec)
        
        boundsep[s] <- exp(logBoundsepMean + logBoundsepCond[condition[s]] + logBoundsepSubj[s])
    }
    
    for (c in 1:nCondition) {
        boundsepMean[c] <- inprod(boundsep, inCondition[c,]) / nInCondition[c]
        boundsepSD[c] <- sqrt(inprod((boundsep - boundsepMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Startpoint variability
    
    logStartvarMean ~ dnorm(0, 0.001)
    
    logStartvarCondPrec ~ dgamma(1.001, 0.001)
    logStartvarSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logStartvarCond[c] ~ dnorm(0, logStartvarCondPrec)
    }
    
    for (s in 1:nSubject) {
        logStartvarSubj[s] ~ dnorm(0, logStartvarSubjPrec)
        
        startvar[s] <- exp(logStartvarMean + logStartvarCond[condition[s]] + logStartvarSubj[s])
    }
    
    for (c in 1:nCondition) {
        startvarMean[c] <- inprod(startvar, inCondition[c,]) / nInCondition[c]
        startvarSD[c] <- sqrt(inprod((startvar - startvarMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Bias
    
    logitBiasMean ~ dnorm(0, 0.001)
    
    logitBiasCondPrec ~ dgamma(1.001, 0.001)
    logitBiasSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logitBiasCond[c] ~ dnorm(0, logitBiasCondPrec)
    }
    
    for (s in 1:nSubject) {
        logitBiasSubj[s] ~ dnorm(0, logitBiasSubjPrec)
        
        bias[s] <- ilogit(logitBiasMean + logitBiasCond[condition[s]] + logitBiasSubj[s])
    }
    
    for (c in 1:nCondition) {
        biasMean[c] <- inprod(bias, inCondition[c,]) / nInCondition[c]
        biasSD[c] <- sqrt(inprod((bias - biasMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Eta
    
    logEtaMean ~ dnorm(0, 0.001)
    
    logEtaCondPrec ~ dgamma(1.001, 0.001)
    logEtaSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logEtaCond[c] ~ dnorm(0, logEtaCondPrec)
    }
    
    for (s in 1:nSubject) {
        logEtaSubj[s] ~ dnorm(0, logEtaSubjPrec)
        
        eta[s] <- exp(logEtaMean + logEtaCond[condition[s]] + logEtaSubj[s])
    }
    
    for (c in 1:nCondition) {
        etaMean[c] <- inprod(eta, inCondition[c,]) / nInCondition[c]
        etaSD[c] <- sqrt(inprod((eta - etaMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Residual time
    
    logitTerMean ~ dnorm(0, 0.001)
    
    logitTerCondPrec ~ dgamma(1.001, 0.001)
    logitTerSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logitTerCond[c] ~ dnorm(0, logitTerCondPrec)
    }
    
    for (s in 1:nSubject) {
        logitTerSubj[s] ~ dnorm(0, logitTerSubjPrec)
        
        ter[s] <- ilogit(logitTerMean + logitTerCond[condition[s]] + logitTerSubj[s]) * 5
    }
    
    for (c in 1:nCondition) {
        terMean[c] <- inprod(ter, inCondition[c,]) / nInCondition[c]
        terSD[c] <- sqrt(inprod((ter - terMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Contaminant probability
    
    logitPContamMean ~ dnorm(0, 0.001)
    
    logitPContamCondPrec ~ dgamma(1.001, 0.001)
    logitPContamSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logitPContamCond[c] ~ dnorm(0, logitPContamCondPrec)
    }
    
    for (s in 1:nSubject) {
        logitPContamSubj[s] ~ dnorm(0, logitPContamSubjPrec)
        
        pContam[s] <- ilogit(logitPContamMean + logitPContamCond[condition[s]] + logitPContamSubj[s])
    }
    
    for (c in 1:nCondition) {
        pContamMean[c] <- inprod(pContam, inCondition[c,]) / nInCondition[c]
        pContamSD[c] <- sqrt(inprod((pContam - pContamMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Beta
    
    logBetaMean ~ dnorm(0, 0.001)
    
    logBetaCondPrec ~ dgamma(1.001, 0.001)
    logBetaSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logBetaCond[c] ~ dnorm(0, logBetaCondPrec)
    }
    
    for (s in 1:nSubject) {
        logBetaSubj[s] ~ dnorm(0, logBetaSubjPrec)
        
        beta[s] <- exp(logBetaMean + logBetaCond[condition[s]] + logBetaSubj[s])
    }
    
    for (c in 1:nCondition) {
        betaMean[c] <- inprod(beta, inCondition[c,]) / nInCondition[c]
        betaSD[c] <- sqrt(inprod((beta - betaMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # U
    
    logUMean ~ dnorm(0, 0.001)
    
    logUCondPrec ~ dgamma(1.001, 0.001)
    logUSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logUCond[c] ~ dnorm(0, logUCondPrec)
    }
    
    for (s in 1:nSubject) {
        logUSubj[s] ~ dnorm(0, logUSubjPrec)
        
        u[s] <- exp(logUMean + logUCond[condition[s]] + logUSubj[s])
    }
    
    for (c in 1:nCondition) {
        uMean[c] <- inprod(u, inCondition[c,]) / nInCondition[c]
        uSD[c] <- sqrt(inprod((u - uMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # V
    
    logVMean ~ dnorm(0, 0.001)
    
    logVCondPrec ~ dgamma(1.001, 0.001)
    logVSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logVCond[c] ~ dnorm(0, logVCondPrec)
    }
    
    for (s in 1:nSubject) {
        logVSubj[s] ~ dnorm(0, logVSubjPrec)
        
        v[s] <- exp(logVMean + logVCond[condition[s]] + logVSubj[s])
    }
    
    for (c in 1:nCondition) {
        vMean[c] <- inprod(v, inCondition[c,]) / nInCondition[c]
        vSD[c] <- sqrt(inprod((v - vMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Background memory strength
    
    logBackgroundMean ~ dnorm(0, 0.001)
    
    logBackgroundCondPrec ~ dgamma(1.001, 0.001)
    logBackgroundSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logBackgroundCond[c] ~ dnorm(0, logBackgroundCondPrec)
    }
    
    for (s in 1:nSubject) {
        logBackgroundSubj[s] ~ dnorm(0, logBackgroundSubjPrec)
        
        background[s] <- exp(logBackgroundMean + logBackgroundCond[condition[s]] + logBackgroundSubj[s])
    }
    
    for (c in 1:nCondition) {
        backgroundMean[c] <- inprod(background, inCondition[c,]) / nInCondition[c]
        backgroundSD[c] <- sqrt(inprod((background - backgroundMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
    
    # Similarity
    
    logitSimMean ~ dnorm(0, 0.001)
    
    logitSimCondPrec ~ dgamma(1.001, 0.001)
    logitSimSubjPrec ~ dgamma(1.001, 0.001)
    
    for (c in 1:nCondition) {
        logitSimCond[c] ~ dnorm(0, logitSimCondPrec)
    }
    
    for (s in 1:nSubject) {
        logitSimSubj[s] ~ dnorm(0, logitSimSubjPrec)
        
        sim[s] <- ilogit(logitSimMean + logitSimCond[condition[s]] + logitSimSubj[s])
    }
    
    for (c in 1:nCondition) {
        simMean[c] <- inprod(sim, inCondition[c,]) / nInCondition[c]
        simSD[c] <- sqrt(inprod((sim - simMean[c])^2, inCondition[c,]) / nInCondition[c])
    }
}
"

# -----------------------------------------------------------------------------
# DATA
# -----------------------------------------------------------------------------

load('../data/nosofsky_sternberg_data_raw.rdata')

data$rt <- data$rt / 1000

data <- data[data$block > 1 & data$rt > 0.2 & data$rt < 5 & data$setsize > 1,]

maxSetSize <- max(data$setsize)
data$lag[data$lag == 0] <- maxSetSize + 1

subject <- as.numeric(factor(data$subject))
setSize <- data$setsize
lag <- data$lag
condition <- as.numeric(tapply(data$condition, INDEX=subject, FUN=function(x) x[1]))

rt <- data$rt * (2 * data$resp - 1)

corrResp <- data$studied

nCondition <- max(condition)
nSubject <- max(subject)
nData <- nrow(data)

datalist <- list(
    rt=rt, corrResp=corrResp,
    maxSetSize=maxSetSize,
    subject=subject, setSize=setSize, lag=lag,
    condition=condition,
    nCondition=nCondition, nSubject=nSubject, nData=nData
)


# -----------------------------------------------------------------------------
# RUN MODEL
# -----------------------------------------------------------------------------

parNames <- c('boundsep', 'startvar', 'bias', 'ter', 'pContam', 'eta', 'beta', 'u', 'v', 'background', 'sim')

parameters <- c(
    paste0(parNames, 'Mean'), paste0(parNames, 'SD'),
    'deviance', 'predRT'
)

writeLines(modelstring, con='model.txt')

adaptSteps <- 100
burnInSteps <- 100
nChains <- 5
thinSteps <- 1
numSteps <- 100

if (nChains == 1) {
    jModel <- jags.model('model.txt', data=datalist, n.chains=nChains, n.adapt=adaptSteps)
    update(jModel, n.iter=burnInSteps)
    samples <- coda.samples(jModel, variable.names=parameters, n.iter=numSteps*thinSteps, thin=thinSteps)
} else {
    samples <- as.mcmc.list(run.jags(model = 'model.txt', monitor=parameters, data=datalist, n.chains=nChains, burnin=burnInSteps, sample=numSteps, adapt=adaptSteps, modules=c('lba', 'dic'), thin=thinSteps, method='parallel', summarise=FALSE, silent.jags=FALSE))
}

# -----------------------------------------------------------------------------
# ANALYZE
# -----------------------------------------------------------------------------

plot(samples[, c(grep('Mean', varnames(samples), value=TRUE), grep('SD', varnames(samples), value=TRUE), 'deviance')], ask=TRUE, density=FALSE, smooth=FALSE)

# Summary plots

samples <- as.matrix(samples)

source('~/Documents/MaPLab/plotPost-mod.R')
source('~/Documents/MaPLab/plotmeansggplot.r')

for (par in parNames) {
    x11(width=7, height=6)
    par(mfrow=c(nCondition, nCondition))
    for (i in 1:nCondition) {
        for (j in 1:nCondition) {
            if (j < i) {
                plot.new()
            } else if (j == i) {
                plotPost(samples[,paste0(par, 'Mean[',i, ']')], xlab=paste(par, 'mean'), main=levels(data$condition)[i], showMode='mode')
            } else {
                plotPost(samples[,paste0(par, 'Mean[',i, ']')] - samples[,paste0(par, 'Mean[',j, ']')], xlab='', main=paste(levels(data$condition)[i], levels(data$condition)[j], sep=' - '), compVal=0, showMode='mode')
            }
        }
    }
}

# Posterior predictive

predCorr <- matrix(1, nrow=nrow(samples), ncol=nData)
predCorr[samples[, paste('predRT[', 1:nData, ']', sep='')] < 0] <- NA
predErr <- matrix(1, nrow=nrow(samples), ncol=nData)
predErr[samples[, paste('predRT[', 1:nData, ']', sep='')] > 0] <- NA

data$predPCorr <- colMeans(samples[, paste0('predRT[',1:nData,']')] > 0)
data$predCorrRT <- colMeans(samples[, paste0('predRT[',1:nData,']')] * predCorr, na.rm=TRUE)
data$predErrRT <- colMeans(samples[, paste0('predRT[',1:nData,']')] * predErr, na.rm=TRUE)

data$setsizefactor <- factor(data$setsize)

x11(width=11, height=4.5)
print(plot.means.ggplot('corr', 'lag', 'setsizefactor', 'condition', 'subject', data=data) + ylim(0, 1), newpage=F)
x11(width=11, height=4.5)
print(plot.means.ggplot('predPCorr', 'lag', 'setsizefactor', 'condition', 'subject', data=data) + ylim(0, 1), newpage=F)

x11(width=11, height=4.5)
print(plot.means.ggplot('rt', 'lag', 'setsizefactor', 'condition', 'subject', data=data[data$corr==1,]), newpage=F)
x11(width=11, height=4.5)
print(plot.means.ggplot('predCorrRT', 'lag', 'setsizefactor', 'condition', 'subject', data=data), newpage=F)

x11(width=11, height=4.5)
print(plot.means.ggplot('rt', 'lag', 'setsizefactor', 'condition', 'subject', data=data[data$corr==0,]), newpage=F)
x11(width=11, height=4.5)
print(plot.means.ggplot('predErrRT', 'lag', 'setsizefactor', 'condition', 'subject', data=data), newpage=F)

