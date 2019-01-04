# CrossValidation_Functions.R
library(MASS)
library(leaps)
library(R.utils) # for doCall()
library(plyr) # for adply()
library(ggplot2)
library(reshape2) # for dcast()
library(xtable)
library(doParallel)

# Note:
# In multiMethodSim(), we had a typo until 6/17/16
# of using mvrnorm() instead of rnorm() to generate epsilon.
# Problematic because mvrnorm takes variance,
# yet we were giving it sd like rnorm() takes.
# It's fixed now, but won't be backwards-compatible.
# So if we wanted to replicate results from an older run,
# we'd have to undo this fix temporarily.
# But honestly, we should just rerun new sims of results we care about,
# instead of repeating these incorrect sims & results.

# Note:
# In code below, beta refers to the vector of non-intercept terms;
# and betaWithInt is the vector augmented with B_0 in the first place.
# Similarly for betahat and betahatWithInt.



#### Simulation wrapper functions ####

# Wrapper *around* the sim-wrapper,
# to run many different combos of settings.
# Give it a *list* of data settings (n, p, etc.),
#   which it will expand into an all-cross-combinations dataframe;
# and a *string vector* of function names to use as procedures,
#   ideally taking just 2 args: functionName(Y, X),
#   with all other args set as defaults;
# and tell it whether to useParallel processing, and if so how many cores.
# 
# It returns the simulated dataset / results,
# and as a side effect it dumps the dataset to a file
# (along with all input settings and the random seed)
runAndSaveSim = function(indir, outdir,
                         simSettingsList, whichProcedures,
                         useParallel = FALSE, nCores = NULL,
                         mySeed = NULL, mySuffix = NULL,
                         simSettingsDf = NULL) {
  
  if(is.null(simSettingsDf)) {
    # Convert list of options into a dataframe with all cross-combinations
    simSettingsDf = do.call(expand.grid0, simSettingsList)
  }
  
  if(is.null(mySeed)) {
    # Make a random seed based on the date, if unspecified
    mySeed = as.integer(format(Sys.time(), "%y%m%d"))
  }
  if(is.null(mySuffix)) {
    # Make a string of the date, to use as a file name suffix
    mySuffix = format(Sys.time(), "%y-%m-%d")
  }

  ### SPECIFY PROCEDURES, CREATE WRAPPER
  # Create a wrapper around doCall(multiSim, ...)
  # so that it runs once per row of the simSettings dataframe;
  # also specify whichProcedures here,
  # so that each replicated dataset runs every procedure
  multiSimPly = function(dfrow, indir, whichProcedures) {
    source(paste0(indir, "CrossValidation_Functions.R"))
    doCall("multiSim",
           whichProcedures = whichProcedures,
           alwaysArgs = as.list(dfrow))
  } 
  
  ### SET UP OPTIONS FOR adply(), INCL. WHETHER PARALLEL OR NOT
  if(useParallel) {
    packageVec = c("MASS", "leaps", "R.utils") # are these all we need to pass on?
    stopifnot(nCores <= detectCores())
    cl = makeCluster(nCores)
    registerDoParallel(cl)
    clusterExport(cl, whichProcedures) # in case we define global fns ouside of ...Functions.R file
    clusterSetRNGStream(cl, mySeed) # parallel equiv of set.seed()
    runSim = function(df) {
      adply(df, .margins = 1, .fun = multiSimPly,
            indir = indir, whichProcedures = whichProcedures,
            .progress = "none",
            .parallel = TRUE,
            # preschedule = TRUE makes reproducible,
            # though no load-balancing so may be slower
            .paropts = list(.options.snow = list(preschedule = TRUE),
                            .packages = packageVec))
    }
  } else { # if not running in parallel
    set.seed(mySeed)
    runSim = function(df) {
      adply(df, .margins = 1, .fun = multiSimPly,
            indir = indir, whichProcedures = whichProcedures,
            .progress = "text",
            .parallel = FALSE)
    }
  }
  
  ### RUN AND SAVE SIMULATION
  # Use adply to run complete simulation
  # (run each row of simSettings and rbind the results)
  print(system.time({
    SimData = runSim(simSettingsDf)
  }))
  
  if(useParallel){
    stopCluster(cl)
  }
  
  # Dump the results, as well as the random seed and all the settings
  # (unsure of best way to dump any extra proc-functions defined for this run)
  dump(c("SimData", "mySeed",
         "indir", "outdir",
         "simSettingsList", "whichProcedures",
         "useParallel", "nCores"),
       file = paste0(outdir, "SimData_", mySuffix, ".R"))

  return(SimData)
}



# Restructured sim wrapper to...
# (1) call named functions for generating beta, X, and epsilon,
# (2) output a list with not just means but also SDs over nrSims
multiSim = function(nrSims = 100, n = 115, p = 6, K = 5,
                    beta.fun = "beta.Linear", betamin = 5, betamax = 9, betarange = NULL,
                    Sigma.fun = "Sigma.Identity", mu = 1/(2*K), fac = 1 - mu,
                    X.fun = "X.Normal",
                    epsilon.fun = "epsilon.Normal",
                    whichProcedures = c("CrossValFS", "OMP"),
                    ...){
  # so that we don't have to deal with NAs for SD
  # and dimension-dropping for rowMeans etc...
  stopifnot(nrSims > 1)
  
  # By specifying betarange, you override betamax
  # (?is there a better way to do this?)
  if(!is.null(betarange)) {
    stopifnot(betarange >= 0)
    betamax = betamin + betarange
  }
  
  # We can use names(list(...)) to see what got passed in via ...
  # and this way we don't need to specify threshold=NULL in default args.
  # 
  # TODO: check if this is the right way to deal with ... arguments?
  # Without attach, it'll check the names but still won't see the values;
  # without on.exit(detach), it'll stay attached even after function exits.
  # Or would attachLocally() be better for any reason?
  # 
  # ...Can't we just deal with this inside the OMP function instead?
  # print(names(list(...)))
  # attach(list(...)); on.exit(detach(list(...)))
  # if("OMP" %in% whichProcedures){
  #   if(! "threshold" %in% names(list(...))) {
  #     threshold = NULL
  #   }
  #   if(is.null(threshold)){
  #     if(! "eta" %in% names(list(...))) {
  #       eta = NULL
  #     }
  #     if(is.null(eta)){
  #       eta = 1
  #     }
  #     # Use threshold with sigma=1, assuming known;
  #     # right now our epsilon.fun functions assume sigma=1 too.
  #     sigma = 1
  #     threshold = sigma*sqrt(2*(1+eta)*log(p)/n)
  #   }
  #   stopifnot(threshold >= 0)
  # } else {
  #   # If OMP is *not* one of the functions, we still want to define threshold,
  #   # so we can pass it below.
  #   # I'd *rather* append this to "..." but not sure how or if that's even possible.
  #   threshold = NULL
  # }
  
  # Assuming for now that we generate Sigma once per call to multiSim(),
  # but might be feeding in different functions for it
  # in the wrapper *around* multiSim()
  Sigma = doCall(Sigma.fun, p = p, mu = mu, fac = fac)
  
  # Names of the outcomes reported by runAndEvalProcedure()
  outcomeNames = c("foundTrueModel", 
                    "nrFalsePos", "nrFalseNegs",
                    "estimationSSE", "elapsed")
  nrOutcomes = length(outcomeNames)
  nrProcedures = length(whichProcedures)
  
  # Use sapply() to mimic replicate() but allow ... to be passed
  simReps = sapply(1:nrSims, function(...) { 
    beta = doCall(beta.fun, betamin = betamin, betamax = betamax, K = K, p = p)
    betaWithInt = c(0, beta)
    X = doCall(X.fun, n = n, p = p, Sigma = Sigma, ...)
    epsilon = doCall(epsilon.fun, n = n, ...)
    Y = X%*%beta + epsilon
    
    # Within each sim replication,
    # we're applying each procedure to the same dataset,
    # *not* generating a whole new dataset for each procedure
    results = matrix(0, nrProcedures * nrOutcomes, 3)
    colnames(results) = c("procedure", "outcome", "value")
    for(ii in 1:nrProcedures) {
      oneResult = runAndEvalProcedure(Y, X, 1:K, betaWithInt,
                                      procedure = whichProcedures[ii],
                                      ...)
      results[((ii-1)*nrOutcomes + 1):(ii*nrOutcomes), ] =
        cbind(ii, 1:nrOutcomes, unname(oneResult))
    }
    return(results)
  }, ..., simplify = "array") # mimic replicate() but allow ... to be passed
  
  valueIndex = which(colnames(simReps) == "value")
  values = simReps[, valueIndex, ]
  simSummary = data.frame(simReps[, -valueIndex, 1],
                          mean = rowMeans(values),
                          sd = rowSDs(values))
  simSummary$procedure = factor(simSummary$procedure,
                                     levels = 1:nrProcedures,
                                     labels = whichProcedures)
  simSummary$outcome = factor(simSummary$outcome,
                                   levels = 1:nrOutcomes,
                                   labels = outcomeNames)
  return(simSummary)
}



# Setup function to run some multi-method simulations:
# works with both FS+CV and OMP+CW8 (or either alone).
# (Idea is that we can run multiple methods on the same dataset,
#  for every simulated dataset.)
# Give it input parameters;
# create one (deterministic) beta vector;
# then repeatedly generate new (random) X matrices and epsilon vectors,
# and fit them using FS+CV and/or OMP+CW8.
multiMethodSim = function(nrSims = 100, n = 115, p = 6, K = 5,
                     betamin = 5, mu = 0.1,
                     sigma = 1, Sigma = diag(p),
                     whichProcedures = c("CrossValFS", "OMP"),
                     nFolds = 5, threshold = NULL,
                     trainProportion = 0.5, nSplits = 1){
  if("OMP" %in% whichProcedures){
    if(is.null(threshold)){
      # Use threshold with eta=1, assuming sigma known
      eta = 1
      threshold = sigma*sqrt(2*(1+eta)*log(p)/n)
    }
    stopifnot(threshold >= 0)
  }
  
  # Let beta decrease by integer values until it reaches betamin,
  # in K'th position, then 0s afterwards
  # (Jing wants this to be random,
  #  but that seems like it'll add noise to sims w/o clarifying much,
  #  esp. since neither our proofs nor C&W's account for random beta...)
  beta = c((betamin+K-1):betamin, rep(0, p-K))
  betaWithInt = c(0, beta)
  
  simReps = replicate(nrSims, {
    # Generate a new draw of random X for each sim;
    # and normalize it BEFORE generating Y
    X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    X = normalize(X)
    # Going back to using sigma^2 / n for the noise variance
    ## OH NO! We've been doing it wrong: mvrnorm uses Sigma matrix i.e. variance, not sd
    # epsilon = mvrnorm(n, 0, sigma/sqrt(n))
    ## Replace with rnorm and rerun to see what happens to our results
    ## (checking this as of 8:30pm, 6/17/16)
    epsilon = rnorm(n, mean = 0, sd = sigma/sqrt(n))
    Y = X%*%beta + epsilon
    results = vector("numeric", 0)
    for(ii in 1:length(whichProcedures)) {
      results = rbind(results,
                      runAndEvalProcedure(Y, X, 1:K, betaWithInt,
                                          procedure = whichProcedures[ii],
                                          threshold = threshold,
                                          nFolds = nFolds,
                                          trainProportion = trainProportion,
                                          nSplits = nSplits))
    }
    return(results)
  })
  if(length(whichProcedures) == 1){
    return(rowMeans(simReps, dim = 2))
  } else {
    return(cbind(whichProcedure = 1:length(whichProcedures),
                 rowMeans(simReps, dim = 2)))
  }
}



# Setup function to run some simple simulations:
# give it input parameters,
# create one (deterministic) beta vector and (random) X matrix,
# then repeatedly generate new Y vectors and fit them
# using FS with K chosen by cross-validation
simpleSim = function(nrSims = 100, n = 115, p = 6, K = 5,
                     betamin = 5, mu = 0.1, sigma = 1){
  
  # Let beta decrease by integer values until it reaches betamin,
  # in K'th position, then 0s afterwards
  beta = c((betamin+K-1):betamin, rep(0, p-K))
  betaWithInt = c(0, beta)
  
  # What's a reasonable covariance matrix for X
  # to have many values around the cutoff for mu, but still pos. def.?
  # Trying a Toeplitz matrix with +mu on first off-diagonal
  # and -mu on second off-diagonal,
  # then 0s beyond that.
  C = toeplitz(c(1, mu, -mu, rep(0, p-3)))
  #X = mvrnorm(n = n, mu = rep(0, p), Sigma = C, empirical = TRUE)
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = C)
  
  simReps = replicate(nrSims, {
    #epsilon = mvrnorm(n, 0, sigma, empirical = TRUE)
    epsilon = mvrnorm(n, 0, sigma)
    Y = X%*%beta + epsilon
    return(runAndEvalProcedure(Y, X, 1:K, betaWithInt))
  })

  return(rowMeans(simReps))
}
# Later we could try as in Yu & Zhao
# and not *set* conditions like mu and gamma,
# just try a wide range of X matrices and *find* these conditions,
# then *plot* them... e.g. PctFoundTrueModel vs mu and vs gamma?
# So our functions would have to report mu and gamma for each rep,
# along with the outcomes, instead of only reporting avg over reps.



# Function to evaluate the output of our FS+CV procedure
# or of our implementation of Cai & Wang's OMP+stopping rule procedure.
# Takes modelOutput = list(betahatWithInt, modelterms),
# trueModel = indices of nonzero coefs (where intercept's index is 0),
# and trueBetaWithInt = c(Intercept, beta).
# Reports whether or not true model was chosen;
# nr of false positives and false negatives;
# and L2 estimation error (MSE) of coefficients
evalProcedure = function(procOutput, trueModel, trueBetaWithInt){
  p = length(trueBetaWithInt) - 1
  foundTrueModel = identical(procOutput$modelterms, trueModel)
  nrFalsePos = sum(!procOutput$modelterms %in% trueModel)
  nrFalseNegs = sum(!trueModel %in% procOutput$modelterms)
  # Whoops, I had been calculating SUM instead of MEAN squared error...
  # Gotta rerun/ignore any previously-saved tables of results...?
  # No, I should just stick with SSE, and rename tables.
  # SSE is more relevant, since the extra variables (larger denom in mean)
  # are *always* all true 0s (in our sims),
  # so a drop in MSE doesn't mean things got better---
  # whereas a drop, vs flat, vs rise in SSE are all meaningful.
  estimationSSE = sum((trueBetaWithInt - procOutput$betahatWithInt)^2)
  return(c(foundTrueModel = foundTrueModel,
           nrFalsePos = nrFalsePos, nrFalseNegs = nrFalseNegs,
           estimationSSE = estimationSSE))
}



# Function to run and evaluate model performance;
# also reporting the runtime.
# If the procedure takes extra arguments
# (like OMP needs threshold, or CrossValFS will take nFolds),
# just pass them along and ... will take care of it;
# but thanks to R.utils::doCall(),
# it won't break if the procedure doesn't need that extra argument
runAndEvalProcedure = function(Y, X, trueModel, trueBetaWithInt,
                               procedure = "CrossValFS", ...){
  timing = system.time({ procOutput = doCall(procedure,
                                             Y = Y, X = X, ...,
                                             .ignoreUnusedArgs = TRUE) })
  return(c(evalProcedure(procOutput, trueModel, trueBetaWithInt),
           timing["elapsed"]))
}






#### FS and cross-validation functions ####

# CrossValFS is now just a wrapper for FS_CV_K(), for backwards compatibility
CrossValFS = function(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1){
  FS_CV_K(Y, X, nFolds, trainProportion, nSplits)
}

# Write separate functions for
# FS+CV over K (choose model size)
# vs
# FS+CV over ResCor (choose threshold for corr b/w residuals and X columns)



# Function to run Forward Stepwise,
# including choosing K (nr of non-intercept terms to include)
#   by cross-validation (or sample-splitting if nFolds = 1).
# Can rerun CV or SS and take majority vote over multiple splits,
#   if nSplits > 1.
# Returns integer vector showing indices of which predictors were chosen,
# and numeric vector same shape as c(Intercept, beta) with estimated coefs.
FS_CV_K = function(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1){
  Khats = replicate(nSplits, CrossValChooseK(Y, X, nFolds, trainProportion))
  Khat = Mode(Khats)
  dataset = data.frame(Y = Y, X)
  full.formula = formula(paste0("Y~", paste0(colnames(dataset)[-1], collapse="+")))
  
  if(Khat == 0) {
    estcoefs = coef(lm(Y ~ 1, dataset))
  } else {
    estcoefs = coef(regsubsets(full.formula, dataset,
                               method = "forward", nvmax = Khat),
                    id = Khat)
  }

  return(tidyEsts(estcoefs, p = ncol(X)))
}



# Same as FS_CV_K(), but instead of using CV to choose K directly,
# use CV to choose a threshold for the ResCorr:
# the max absolute correlation between current model's residual and any X column.
# Choose a threshold by CV, then go back to full model
# and use the first model along path whose ResCorr is below that threshold.
# (Analogous to C&W's threshold for OMP, but
#  (1) chosen by CV instead of theory, and
#  (2) normalized correlation, not just inner product,
#      so it's NOT necessarily monotonic along model path...)
# 
# WAIT: do we really want |Corr(R_t, X_i)| ?
# It should really be |Corr(R(Y|X_t), R(X_i|X_t))| instead.
# Would that be monotonic instead?
FS_CV_ResCorr = function(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1){
  # EVENTUALLY:
  # Run full model path on full dataset
  # Get the (monotonic-ized) fullResCorr path
  # Use SS or CV to train model paths;
  #   find and test the models along fullResCorr path;
  #   and choose the fullResCorr that gave best (mean or modal) testRSS
  # Return to full dataset and get the chosen ResCorr's model estcoefs
  # Pass them to tidyEsts and return
  
  # BUT FOR NOW, until we get that working...
  # just return a hardcoded null estimated model,
  # so we can test whether interface works even if internals don't
  tidyEsts(estcoefs = c("(Intercept)" = 0),
           p = ncol(X))
}



# FS model chosen by algorithm at start of Miller (2002), Ch.4:
# Generate new artificial X variables as random uniform draws.
# (By default, make as many artificial X's as real ones.)
# Run FS model path and stop as soon as an artifical variable is chosen;
# select the largest model with no artificial variables.
# (TODO: also allow "nSplits" to be used here;
#  though instead of really being new splits,
#  it's the number of times to generate a new artificial dataset
#  and re-choose model size.)
# (TODO: handle the possible edge case, if p > n,
#  that NO fake Xs are chosen... What will which.max do then?)
FS_M4 = function(Y, X, pNew = ncol(X), nSplits = 1) {
  n = nrow(X); p = ncol(X)
  fake.formula = formula(paste0("Y~",
                                paste0(paste0("X", 1:(p + pNew)),
                                       collapse="+")))
  real.formula = formula(paste0("Y~",
                                paste0(paste0("X", 1:p),
                                       collapse="+")))
  
  Khats = vector("integer", length = nSplits)
  for(ss in 1:nSplits) {
    Xfake = normalize(matrix(runif(n * pNew), n))
    dataset = data.frame(Y = Y, cbind(X, Xfake))
    leaps.out = regsubsets(fake.formula, dataset, method = "forward",
                           nvmax = min(n, p + 1)) # Needn't go more than p steps
    # Find the first model that accepted a fake X:
    # find the first row of summary(leaps.out)$which
    # that has TRUE in any of the last pNew columns.
    # Then subtract 1 to go back a step
    anyFakes = apply(summary(leaps.out)$which[, tail(1:(1+p+pNew), pNew)],
                     1, any)
    if(any(anyFakes)) {
      Khats[ss] = which.max(anyFakes) - 1
    } else {
      # If we never chose a fake variable (possible if p > n),
      # use the last model in the path
      Khats[ss] = length(anyFakes)
    }
  }
  
  Khat = Mode(Khats)
  # Not necessarily OK to reuse last run of model:
  # the first Khat steps *might* include fake Xs
  # if we used multiple "splits" and the last one wasn't a winner...
  if(nSplits > 1) {
    dataset = data.frame(Y = Y, X)
    leaps.out = regsubsets(real.formula, dataset, method = "forward", nvmax = Khat)
  }
  if(Khat == 0) {
    estcoefs = coef(lm(Y ~ 1, dataset))
  } else {
    estcoefs = coef(leaps.out, id = Khat)
  }
  return(tidyEsts(estcoefs, p = ncol(X)))
}


# Common code for FS_...() functions,
# cleaning up estimated coefficients
# (just the selected terms, coming out of lm or regsubsets)
# and turning them into a list of two items:
# the full coef vector including 0s for unchosen terms,
# and the integer vector of selected model terms' indices
tidyEsts = function(estcoefs, p) {
  # Which predictors were chosen?
  # Return a sorted integer vector of indices,
  # e.g. for (X3, X1) we'd return c(1, 3)
  modeltermNames = names(estcoefs[-1]) # ignore intercept
  modelterms = sort(as.integer(substr(modeltermNames,
                                      start = 2,
                                      stop = nchar(modeltermNames))))
  
  # What was the complete estimated betaWithInt vector?
  # Return a vector same shape as c(Intercept, beta),
  # with 0s for unchosen terms
  betahatWithInt = matrix(rep(0, p+1), 1,
                          dimnames = list(NULL,
                                          c("(Intercept)",
                                            paste0("X", 1:p))))
  betahatWithInt[, names(estcoefs)] = estcoefs
  
  return(list(betahatWithInt = betahatWithInt,
              modelterms = modelterms))
}



# Function to choose K by cross-validation
# (though for now, just implements sample splitting)
# Each time CrossValChooseK() runs,
# the sample-split/fold assignment is randomized;
# so if desired, we can re-run it many times
# to vote on best choice of K after many re-splits.
CrossValChooseK = function(Y, X, nFolds = 1, trainProportion = 0.5){
  n = length(Y)

  # If nFolds = 1, just do simple sample-splitting
  if(nFolds == 1){
    # Split in half, or whatever training proportion is desired
    # (80/20 split is more common,
    #  but keeping 50/50 default for back-compatibility for now)
    nTrain = ceiling(n * trainProportion)
    iis = sample(n)
    iiTrain = iis[1:nTrain]
    iiTest = iis[(nTrain+1):n]
    # Train once
    betahatsWithInt = trainFS(Y[iiTrain], X[iiTrain, ])
    # Test once
    RSSs = testFS(Y[iiTest], X[iiTest, ], betahatsWithInt)
    K = which.min(RSSs) - 1 # subtract 1 b/c 1st entry is intercept-only
  } else {
    # Define the folds
    foldIDs = sample(rep(1:nFolds, ceiling(n/nFolds))[1:n])
    RSSs = vector("numeric")
    for(ff in 1:nFolds){
      # Train on the other (nfolds-1) folds
      betahatsWithInt = trainFS(Y[!foldIDs == ff],
                         X[!foldIDs == ff, ])
      # Test on this fold
      RSSs = rbind(RSSs, testFS(Y[foldIDs == ff],
                                X[foldIDs == ff, ],
                                betahatsWithInt))
    }
    # Average test RSS over all folds,
    # and choose k that minimizes it
    K = which.min(colMeans(RSSs)) - 1 # subtract 1 b/c 1st entry is intercept-only
  }
  return(K)
}



# Take a dataset,
# run forward stepwise up through full model,
# create string vector of the formulas for each step,
# and return list of the coefficient vectors (betahatsWithInt) for each step
# from 1 (intercept only) to p+1 (full model)
trainFS = function(Y, X){
  n = nrow(X); p = ncol(X)
  dataset = data.frame(Y = Y, X)

  # Run forward stepwise, from intercept-only to min(n, p) predictors
  full.formula = formula(paste0("Y~", paste0(colnames(dataset)[-1], collapse="+")))

#   step.out = step(lm(Y ~ 1, data = dataset),
#                   scope = full.formula, direction = "forward", trace = 0,
#                   steps = min(n, p),
#                   k = 0) # set k=0 for no penalty, so we add all possible terms
#   # Return the sequence of model formulas,
#   # using the ordered sequence of terms entered in step.out$anova$Step
#   formulas = sapply(1:length(step.out$anova$Step),
#                     function(x) paste0("Y ~ 1",
#                                        paste0(step.out$anova$Step[1:x],
#                                               collapse = " ")))
#   # Compute betahatWithInt for each model in the sequence
#   betahatsWithInt = lapply(formulas,
#                     function(x) lm(x, dataset)$coefficients)
  
  # Speed up by replacing step() and lm()
  # with regsubsets() from leaps package
  leaps.out = regsubsets(full.formula, dataset,
                         method = "forward", nvmax = min(n, p))
#   idmax = ifelse(any(is.na(leaps.out$rss)),
#                  which.min(!is.na(leaps.out$rss)),
#                  length(leaps.out$rss) - 1)
  
  # Allow the possibility that intercept-only model is best
  betahatsWithInt = list(coef(lm(Y~1)))
  # Append the other coefs from leaps.out
  if(leaps.out$nvmax >= 2) {
    betahatsWithInt[2:leaps.out$nvmax] = coef(leaps.out, id=1:(leaps.out$nvmax - 1))
  }
  return(betahatsWithInt)
}



# Take a dataset and list of vectors of coef ests (betahatsWithInt),
# and return computed test RSS of each model.
# (Downstream, take this output
#  and choose the model size K (nr of non-intercept predictors)
#  that minimizes the average test RSS over all folds or splits.)
testFS = function(Y, X, betahatsWithInt){
  Xp1 = as.matrix(data.frame(Int = 1, X))
  colnames(Xp1)[1] = "(Intercept)"
  
  RSSs = sapply(betahatsWithInt,
                function(x) sum((Y - (Xp1[, names(x), drop = FALSE] %*% x))^2))
  names(RSSs) = 0:(length(RSSs)-1)
  return(RSSs)
}






#### OMP functions ####

# OMP function takes in (Y,X) and a stopping rule:
# right now, stopping rule is just the CW8 rule threshold
# (threshold is on max abs inner product
#  between residual and normalized X columns.)
# Produces betahatWithInt and modelterms, same as our FS functions.
OMP = function(Y, X, threshold = NULL, eta = 1, sigma = 1) {
  n = nrow(X); p = ncol(X)
  
  if(is.null(threshold)) {
    threshold = sigma*sqrt(2*(1+eta)*log(p)/n)
  }
  
  # ? Does OMP require zero-mean and unit-norm Y?
  # If we are not estimating an intercept as part of OMP,
  # seems like we should at least de-mean Y and all the Xs.
  # Probably no need to unit-norm Y (scale shouldn't matter there)
  # but also shouldn't hurt.
  Res = normalize(matrix(Y, ncol = 1))
  # Columns of X definitely must be 0-mean, 1-norm.
  Xnormed = normalize(X)
  colnames(Xnormed) = paste0("X", 1:p)
  Xout = 1:p
  Xin = vector("integer")
  InnerProds = abs(t(Res) %*% Xnormed)
  
  # Could probably be more efficient if
  # --we drop rows of Xnormed before calc'ing InnerProds ?
  # --we use qr() instead of resid(lm()) ?
  while(max(InnerProds) > threshold &
        length(Xin) < min(n, p)) {
    best = which.max(InnerProds[])
    Xin = c(Xin, best)
    Xout = Xout[!Xout == best]
    Res = resid(lm(Res ~ Xnormed[, Xin]))
    InnerProds = abs(t(Res) %*% Xnormed)
  }
  
  # Get betahatWithInt for the chosen vars,
  # but also restructure it so it gives us ALL p+1 vars,
  # with 0s as appropriate
  dataset = data.frame(Y = Y, X)
  betahatWithInt = as.matrix(dataset[1, ])
  colnames(betahatWithInt)[1] = "(Intercept)"
  betahatWithInt[1, ] = 0
  omp.formula = formula(paste0("Y~", ifelse(length(Xin) == 0,
                                            "1",
                                            paste0("X", Xin, collapse="+"))))
  estcoefs = coef(lm(omp.formula, dataset))
  betahatWithInt[, names(estcoefs)] = estcoefs

  modelterms = sort(Xin)
  return(list(betahatWithInt = betahatWithInt,
              modelterms = modelterms))
}






#### Small helper functions ####

# Helper function for multiSim(),
# reporting the standard deviation of each result row/column combo
# over the 3rd dimension (nrSims)
facetSDs = function(x) {
  apply(x, 1:2, sd)
}

# L2() calculates L2 norm of a vector
L2 = function(vec) {
  sqrt(sum(vec^2))
}

# normalize() subtracts mean and rescales to unit L2 norm
# for each column of a matrix
# (so if giving it a vector, you must make it a 1-column matrix first)
normalize = function(mat) {
  mat = sweep(mat, 2, colMeans(mat))
  mat = sweep(mat, 2, apply(mat, 2, L2), FUN = "/")
}

# Mode(), to calculate the univariate mode of integer data,
#   or the smallest modal value if there are ties.
# Assumes the inputs are non-negative integers; ought to write a test...
# (R's built-in mode() is about storage type, not central tendency)
Mode = function(vec) {
  # which.max(tabulate(vec)) # does not work for 0 (or negative integers)
  as.integer(names(which.max(table(vec))))
}

# Analogous to rowMeans for 2D matrix:
# within each row, take sd() across columns
rowSDs = function(mat) {
  apply(mat, 1, sd)
}

# Create strings for printing: "mean (SD)"
# with a fixed nr of digits past the decimal
printMeansSDs = function(mean, sd, digits = 2) {
  paste0(formatC(mean, digits, format="f"),
         " (",
         formatC(sd, digits, format="f"),
         ")")
}

# Compute MOEs: by default uses 2*SE, where SE = SD/sqrt(n),
# but you can replace fac for other kinds of symmetric CI,
# e.g. z_alpha = 1.96 for Normal 90% CI
calcMOEs = function(sd, nrSims, fac = 2) {
  fac * sd / sqrt(nrSims)
}

# expand.grid() but with default set to keep strings as strings
# (name inspired by paste0)
expand.grid0 = function(...) {
  expand.grid(..., stringsAsFactors = FALSE)
}


#### Functions to generate beta, Sigma, epsilon ####

# Note again: beta refers to the vector of non-intercept parameters;
# betaWithInt refers to that vector augmented by an intercept term

# Generate Beta: deterministic, null
beta.Null = function(p = 6){
  stopifnot(0 < p)
  return(rep(0, p))
}

# Generate Beta: deterministic, constant
beta.Constant = function(betamin = 5, K = 5, p = 6){
  stopifnot(0 < K & K <= p)
  stopifnot(0 < betamin)
  beta = c(rep(betamin, K),
           rep(0, p - K))
  return(beta)
}

# Generate Beta: deterministic, linear
# (linear decrease from betamax to betamin, then 0s after)
beta.Linear = function(betamin = 5, betamax = betamin + (K-1), K = 5, p = 6){
  stopifnot(0 < K & K <= p)
  stopifnot(0 < betamin & betamin <= betamax)
  beta = c(seq(betamax, betamin, length.out = K),
           rep(0, p - K))
  return(beta)
}

# Generate Beta: random, uniform
# (linear decrease from betamax to betamin, then 0s after)
beta.Uniform = function(betamin = 5, betamax = betamin + (K-1), K = 5, p = 6){
  stopifnot(0 < K & K <= p)
  stopifnot(0 < betamin & betamin <= betamax)
  beta.nonNull = sort(runif(K), decreasing = TRUE)
  beta.nonNull = beta.nonNull - beta.nonNull[K]
  beta.nonNull = beta.nonNull / max(beta.nonNull) * (betamax - betamin) + betamin
  beta = c(beta.nonNull,
           rep(0, p - K))
  return(beta)
}



# Generate Sigma: identity (no correlations)
Sigma.Identity = function(p = 6){
  stopifnot(0 < p)
  return(diag(p))
}

# Generate Sigma: constant off-diagonal correlation,
# ones on diagonal
Sigma.ConstCorr = function(p = 6, mu = 0.1){
  stopifnot(0 < p)
  stopifnot(0 < mu & mu < 1)
  ones = rep(1, p)
  Sigma = mu * ones %*% t(ones) + (1-mu) * diag(p)
  return(Sigma)
}

# Generate Sigma: Toeplitz with power decay
# (ones on diagonal, then constant off-diagonals decreasing geometrically)
Sigma.Toeplitz = function(p = 6, mu = 0.1, fac = 1 - mu){
  stopifnot(0 < p)
  stopifnot(0 < mu & mu < 1)
  stopifnot(0 < fac & fac < 1)
  topRow = c(1, mu * fac^(0:(p-2)))
  Sigma = toeplitz(topRow)
  return(Sigma)
}



# Generate X matrix: Normal with 0 mean and variance Sigma,
# then each column normalized to zero mean and unit L2 length 
# (unless doNormalize = FALSE)
X.Normal = function(n = 115, p = 6, Sigma = Sigma.Identity(p), doNormalize = TRUE){
  stopifnot(1 <= n & 1 <= p)
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  if(doNormalize) {
    X = normalize(X)
  }
  return(X)
}



# Generate epsilon: Normal with 0 mean and variance 1/n
# (to have approximate unit L2 length on average)
# ...unless rescaleByN = FALSE, then variance is sigma^2 (1 by default)
epsilon.Normal = function(n = 115, sigma = 1, rescaleByN = TRUE){
  stopifnot(1 <= n)
  sd = ifelse(rescaleByN, sigma / sqrt(n), sigma)
  return(rnorm(n, mean = 0, sd = sd))
}

# Generate epsilon: t with 3 df,
# rescaled by variance of n such draws
# (to have approximate unit L2 length on average)
# ...unless rescaleByN = FALSE, then variance is sigma^2 (1 by default)
epsilon.t3 = function(n = 115, sigma = 1, rescaleByN = TRUE){
  stopifnot(1 <= n)
  # When df > 2, variance of t is [df / (df - 2)]
  # so here each entry's variance is 3.
  fac = ifelse(rescaleByN, n * 3 / sigma^2, 3 / sigma^2)
  return(rt(n, df = 3) / sqrt(fac))
}

