# CV_Functions.R
library(MASS)
library(leaps)
library(R.utils) # for doCall()
library(plyr) # for adply()
library(ggplot2)
library(reshape2) # for dcast()
library(xtable)
library(doParallel)
library(Matrix)

# Note:
# Since 6/17/16, we've revised code
# so that multiSim CANNOT pass arguments
# to the CV and OMP functions other than the data Y, X.
# Instead, its whichProcedures argument
# must give the names of wrapper-functions
# that are fully-specified except for Y, X,
# so that Y, X are their only arguments.
# For example, write a new wrapper such as
## myCV = function(Y, X) CrossValFS(Y, X, nF = a, nS = B, tr = C)
# and then pass in
## multiSim(---, whichProcedures = "myCV")

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
    source(paste0(indir, "CV_Functions.R"))
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
            # though no load-balancing so may be slower;
            # seems that *consecutive* rows of df get sent
            # to the same core *together*;
            # so to manually spread out the load,
            # try to have big jobs and small jobs adjacent in df,
            # by having the variable which has most-variable-effect
            # (on timing) as the first one in simSettingsList
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
multiSim = function(nrSims = 100,
                    n = 115, p = 6, K = 5,
                    beta.fun = "beta.Linear", betamin = 5, betamax = 9,
                    Sigma.fun = "Sigma.Identity", mu = 1/(2*K), fac = 1 - mu, type = NULL,
                    X.fun = "X.Normal", doNormalize = FALSE,
                    epsilon.fun = "epsilon.Normal", sigma = 1, rescaleByN = FALSE,
                    whichProcedures = NULL,
                    runGC = FALSE){ # runGC = TRUE improves msmt of "elapsed", but slows sims
  # Require multiple sims,
  # so that we don't have to deal with NAs for SD
  # and dimension-dropping for rowMeans etc.
  stopifnot(nrSims > 1)
  
  # Assuming for now that we generate Sigma once per call to multiSim(),
  # but might be feeding in different functions for it
  # in the wrapper *around* multiSim()
  Sigma = doCall(Sigma.fun, p = p, mu = mu, fac = fac, type = type, K = K)
  
  # Names of the outcomes reported by runAndEvalProcedure()
  outcomeNames = c("foundTrueModel", 
                    "nrFalsePos", "nrFalseNegs",
                    "estimationSSE", "elapsed")
  nrOutcomes = length(outcomeNames)
  nrProcedures = length(whichProcedures)
  
  simReps = replicate(nrSims, { 
    beta = doCall(beta.fun, betamin = betamin, betamax = betamax, K = K, p = p)
    betaWithInt = c(0, beta)
    X = doCall(X.fun, n = n, p = p, Sigma = Sigma, doNormalize = doNormalize)
    epsilon = doCall(epsilon.fun, n = n, sigma = sigma, rescaleByN = rescaleByN)
    Y = X%*%beta + epsilon
    
    # Within each sim replication,
    # we're applying each procedure to the same dataset,
    # *not* generating a whole new dataset for each procedure
    results = matrix(0, nrProcedures * nrOutcomes, 3)
    colnames(results) = c("procedure", "outcome", "value")
    for(ii in 1:nrProcedures) {
      oneResult = runAndEvalProcedure(Y, X, 1:K, betaWithInt,
                                      procedure = whichProcedures[ii],
                                      runGC = runGC)
      results[((ii-1)*nrOutcomes + 1):(ii*nrOutcomes), ] =
        cbind(ii, 1:nrOutcomes, unname(oneResult))
    }
    return(results)
  })
  
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
  # Be aware this is SSE, not MSE
  estimationSSE = sum((trueBetaWithInt - procOutput$betahatWithInt)^2)
  return(c(foundTrueModel = foundTrueModel,
           nrFalsePos = nrFalsePos, nrFalseNegs = nrFalseNegs,
           estimationSSE = estimationSSE))
}



# Function to run and evaluate model performance;
# also reporting the runtime.
# If the procedure takes extra arguments
# (like OMP needs threshold, or CrossValFS will take nFolds),
# DO NOT just pass them along for ... to take care of!
# This used to work, but since 6/17/16 we changed code
# so that you MUST specify a procedure wrapper with all options set,
# so that it ONLY takes Y and X, no other arguments!
runAndEvalProcedure = function(Y, X, trueModel, trueBetaWithInt,
                               procedure = NULL, runGC = FALSE){
  timing = system.time({ procOutput = doCall(procedure,
                                             Y = Y, X = X) },
                       gcFirst = runGC)
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


# Run FS path with known (true) value of K
# (NOT estimating it with CV or otherwise);
# and report results in same format as FS_CV_K and other functions
FS_KnownK = function(Y, X, K){
  # If we run this when true K > p,
  # impossible to get the right model, so don't bother
  # (?Is this OK? Should work fine for Prob(correct selection),
  #  but could throw off results for nr FPs or FNs;
  #  but I suppose we'll delete any such simulation results
  #  from reported results anyway.)
  if(K > ncol(X)){
    return(tidyEsts(c("(Intercept)" = 0), p = ncol(X)))
  }
  
  dataset = data.frame(Y = Y, X)
  full.formula = formula(paste0("Y~", paste0(colnames(dataset)[-1], collapse="+")))
  if(K == 0) {
    estcoefs = coef(lm(Y ~ 1, dataset))
  } else {
    estcoefs = coef(regsubsets(full.formula, dataset,
                               method = "forward", nvmax = K),
                    id = K)
  }
  return(tidyEsts(estcoefs, p = ncol(X)))
}



# Function to run Forward Stepwise,
# including choosing K (nr of non-intercept terms to include)
#   by cross-validation (or sample-splitting if nFolds = 1).
# Can rerun CV or SS and take majority vote over multiple splits,
#   if nSplits > 1.
# Returns integer vector showing indices of which predictors were chosen,
# and numeric vector same shape as c(Intercept, beta) with estimated coefs.
#
# Edit 8/2/16: Now allows to decide whether to vote or average
# across Folds and across Splits,
# by whether voteFolds and voteSplits are TRUE or FALSE.
# Can also combine both (first do Folds, then Splits are reps of the Fold process).
# HOWEVER, cannot (makes no sense to) avg over Splits after voting over Folds.
FS_CV_K = function(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1,
                   voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE){
  # If just one split, assume we're "voting" over it,
  # so we use Mode() instead of rowMeansMinM1()
  stopifnot(!(nSplits == 1 & !voteSplits))
  # Cannot vote over Folds but avg over Splits
  stopifnot(!(nFolds > 1 & nSplits > 1 & (!voteSplits) & voteFolds))
  # If avg'ing (not voting) across Splits,
  # need intermed. fns to return RSS values, not a set of Khats.
  returnRSSs = ifelse(voteSplits, FALSE, TRUE)

  CVout = replicate(nSplits, CrossValChooseK(Y, X, nFolds, trainProportion,
                                             returnRSSs, voteFolds, SeqCV))
  if(voteSplits) {
    Khat = Mode(CVout)
  } else { # then returnRSSs was TRUE, so avg them and find best K
    Khat = rowMeansMinM1(CVout, Seq = SeqCV)
  }
  
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
# Each time CrossValChooseK() runs,
# the sample-split/fold assignment is randomized;
# so if desired, we can re-run it many times
# to vote on best choice of K after many re-splits.
#
# 8/2/16: Edit to allow returning (avg'ed) test RSS vector,
# instead of Khat,
# in case we wish to average RSS over splits rather than vote.
CrossValChooseK = function(Y, X, nFolds = 1, trainProportion = 0.5,
                           returnRSSs = FALSE, voteFolds = FALSE,
                           SeqCV = FALSE){
  stopifnot(!(returnRSSs & voteFolds))
  stopifnot(nFolds > 0)
  stopifnot(is.wholenumber(nFolds) | is.wholenumber(1/nFolds))
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
    if(returnRSSs) {
      return(RSSs)
    }
    K = rowMeansMinM1(RSSs, Seq = SeqCV)
  } else {
    # If nFolds is between 0 and 1, it's "inverted" V-fold:
    # train on one fold, test on the others.
    # E.g. if nFolds = 1/3, then use 3 folds.
    if(nFolds < 1) {
      invert = TRUE
      nFolds = round(1/nFolds)
      print(nFolds)
    } else {
      invert = FALSE
    }
    # Define the folds
    foldIDs = sample(rep(1:nFolds, ceiling(n/nFolds))[1:n])
    RSSs = vector("numeric")
    for(ff in 1:nFolds){
      if(invert) {
        # Train on this fold
        betahatsWithInt = trainFS(Y[foldIDs == ff],
                                  X[foldIDs == ff, ])
        # Test on the other (nfolds-1) folds
        RSSs = cbind(RSSs, testFS(Y[!foldIDs == ff],
                                  X[!foldIDs == ff, ],
                                  betahatsWithInt))
      } else {
        # Train on the other (nfolds-1) folds
        betahatsWithInt = trainFS(Y[!foldIDs == ff],
                                  X[!foldIDs == ff, ])
        # Test on this fold
        RSSs = cbind(RSSs, testFS(Y[foldIDs == ff],
                                  X[foldIDs == ff, ],
                                  betahatsWithInt))
      }
    }
    # Average testRSS over all folds
    if(returnRSSs) {
      avgRSSs = rowMeans(RSSs)
      return(avgRSSs)
    }
    # Or, vote on K across folds
    if(voteFolds){
      Kvec = apply(RSSs, 2, function(x) {rowMeansMinM1(x, Seq = SeqCV)} )
      K = Mode(Kvec)
    } else {
      # Or, choose K that minimizes avg testRSS
      K = rowMeansMinM1(RSSs, Seq = SeqCV)
    }
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
  
  rssFun = function(x) {
    sum((Y - (Xp1[, names(x), drop = FALSE] %*% x))^2)
  }
  
  RSSs = sapply(betahatsWithInt, rssFun)
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
    # Are we right to include 1/sqrt(n) here?
    # C&W didn't do that explicitly;
    # they assumed the NORMALIZED-X model has constant noise variance,
    # but we generate RAW-X model with constant noise variance,
    # so scaling sigma by sqrt(n) seems appropriate here.
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

# Now just run OMP with known K,
# for sake of evaluating how often it succeeds
OMP_KnownK = function(Y, X, K) {
  # Code copied from OMP() above;
  # ought to refactor to avoid code copying, eventually
  n = nrow(X); p = ncol(X)
  dataset = data.frame(Y = Y, X)
  betahatWithInt = as.matrix(dataset[1, ])
  colnames(betahatWithInt)[1] = "(Intercept)"
  betahatWithInt[1, ] = 0
  
  if(K == 0) {
    estcoefs = coef(lm(Y ~ 1, dataset))
  } else {
    Res = normalize(matrix(Y, ncol = 1))
    Xnormed = normalize(X)
    colnames(Xnormed) = paste0("X", 1:p)
    Xout = 1:p
    Xin = vector("integer")
    InnerProds = abs(t(Res) %*% Xnormed)
    while(length(Xin) < K) {
      best = which.max(InnerProds[])
      Xin = c(Xin, best)
      Xout = Xout[!Xout == best]
      Res = resid(lm(Res ~ Xnormed[, Xin]))
      InnerProds = abs(t(Res) %*% Xnormed)
    }
    omp.formula = formula(paste0("Y~", ifelse(length(Xin) == 0,
                                              "1",
                                              paste0("X", Xin, collapse="+"))))
    estcoefs = coef(lm(omp.formula, dataset))
  }
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

# Check if (known-to-be symmetric) matrix is positive definite,
# modified from inside of MASS:mvrnorm
isPosDef = function(Sigma) {
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  return(min(ev) > 0)
}

# Given a matrix of testRSS values
# (each column is a split or fold,
#  each row is a model size starting from 0),
# average within rows and
# report the model size with lowest avg testRSS.
# (If Seq=TRUE, return the first local minimum,
#  i.e. smallest model size beyond which the testRSS rises,
#  not the global minimum.)
# Given a vector (just one split), skip the avg'ing.
rowMeansMinM1 = function(mat, Seq = FALSE){
  stopifnot(is.vector(mat) | is.matrix(mat))
  if(length(dim(mat)) == 2) {
    mat = rowMeans(mat)
  }
  if(Seq) {
    # first local min: smallest model size
    # such that the next model's testRSS is *not smaller*
    return(which.max(diff(mat) >= 0) - 1) # subtract 1 b/c 1st entry is intercept-only
  } else {
    return(which.min(mat) - 1) # subtract 1 b/c 1st entry is intercept-only
  }
}

# Check if a number rounds to itself within tolerance.
# Copied from ?is.integer,
# since is.integer checks storage type
# rather than whether a double is "basically" an integer
is.wholenumber = function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
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

# Generate Sigma: perturbed identity,
# with nonzero correlation in just one pair of variables.
# type = c("withinTrue", "withinFalse", "across")
# indicates whether the nonzero corr is
# between two True Xs (1, 2),
# between two False Xs (p-1, p),
# or across True and False Xs (1, p).
Sigma.OneBigEntry = function(p = 6, mu = 0.1,
                             type = "withinTrue"){
  stopifnot(0 < p)
  stopifnot(0 < mu & mu < 1)
  Sigma = Sigma.Identity(p)
  entry = switch(type,
                 withinTrue = c(1, 2),
                 withinFalse = c(p-1, p),
                 across = c(1, p))
  Sigma[entry[1], entry[2]] = mu
  Sigma[entry[2], entry[1]] = mu
  return(Sigma)
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

# Generate Sigma: constant absolute off-diagonal correlation,
# but negative for first K variables and positive for others;
# ones on diagonal.
# (If SigmaApprox=TRUE, start with this matrix and edit to be pos.definite:
#  Jing suggested that if there are any negative eigenvalues,
#  replace them with value of smallest positive eigenvalue,
#  and reconstruct Sigma from that -- but Jerzy discovered that
#  Matrix::nearPD() has a more sophisticated approach built-in.)
# https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/nearPD.html
# https://www.rdocumentation.org/packages/sfsmisc/versions/1.1-1/topics/posdefify
# http://eprints.ma.man.ac.uk/232/01/covered/MIMS_ep2006_70.pdf
Sigma.NegKPosOthers = function(p = 6, mu = 0.1, K = 5, SigmaApprox = FALSE){
  stopifnot(0 < p & 0 < K)
  stopifnot(0 < mu & mu < 1)
  Sigma = matrix(mu, nrow = p, ncol = p)
  Sigma[1:K, 1:K] = -mu
  diag(Sigma) = 1
  if(SigmaApprox) {
    PD = nearPD(Sigma, corr = TRUE, conv.norm.type = "F")
    stopifnot(PD$converged)
    Sigma = as.matrix(PD$mat)
    ## 9/24/2017: Jerzy coded up something similar manually,
    ## but it turns out to be a less-sophisticated special case
    ## of what Matrix::nearPD() already does:
    # eigs = eigen(Sigma, symmetric = TRUE, only.values = FALSE)
    # if(min(eigs$values) <= 0) {
    #   evals = eigs$values
    #   evals[evals <= 0] = min(evals[evals > 0])
    #   evecs = eigs$vectors
    #   Sigma = evecs %*% diag(evals) %*% t(evecs)
    #   Sigma = cov2cor(Sigma)
    # }
    ## Not obvious whether my version or nearPD() is better generally,
    ## but at least in the special case of
    ## Sigma.NegKPosOthers(p=11, mu=5/(2*10), K=10, Approx = T)
    ## the nearPD() version gives negative corrs for 1:K and pos for last,
    ## while the naive version gives near-zero corrs everywhere.
    ## In other situations, neither naive nor nearPD() gives neg. corrs
    ## for 1:K, and I'm not sure how to enforce them otherwise,
    ## but I suppose this is a decent first thing to try.
  } else {
    stopifnot(isPosDef(Sigma))
  }
  return(Sigma)
}
# Also predefine a version with Approx=TRUE
Sigma.NegKPosOthers.Approx = function(p = 6, mu = 0.1, K = 5) {
  Sigma.NegKPosOthers(p, mu, K, SigmaApprox = TRUE)
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

# As above, but now NOT rescaling to unit variance,
# just like epsilon.t2 function below
epsilon.t3Var3 = function(n = 115, sigma = 1, rescaleByN = TRUE){
  stopifnot(1 <= n)
  fac = ifelse(rescaleByN, n / sigma^2, 1 / sigma^2)
  return(rt(n, df = 3) / sqrt(fac))
}



# Generate epsilon: t with 2 df.
# Not going to try rescaling to unit variance here,
# because variance is undefined for t_2.
# But still allow us to scale it up/down by sigma,
# or by sigma/sqrt(n) is rescaleByN is TRUE.
epsilon.t2 = function(n = 115, sigma = 1, rescaleByN = TRUE){
  stopifnot(1 <= n)
  fac = ifelse(rescaleByN, n / sigma^2, 1 / sigma^2)
  return(rt(n, df = 2) / sqrt(fac))
}

