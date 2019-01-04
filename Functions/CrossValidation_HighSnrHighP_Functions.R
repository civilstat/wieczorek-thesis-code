# Simulating small example datasets
# to show that, with high SNR,
# CV with large p helps account for
# post-selection effect of greedy FS



L2 = function(x) sum(x^2)



genData = function(n = 30, p = 2, B0 = 10) {
  X = data.frame(matrix(rnorm(p*n), nrow = n))
  epsilon = rnorm(n)
  Y = B0 + epsilon
  return(list(X = X, Y = Y))
}



# FS with CV (actually just a single sample split for now),
# randomly splitting data into Train and Test, then
# training and then reporting test RSS for comparing just these 3 models:
# 0 (intercept-only model),
# 1Fixed (intercept and X_1, no matter what the other Xs are), and 
# 1Greedy (intercept and X_i, where X_i is greedily chosen for max |corr| with Y)
FS1CV = function(Y, X, nFolds = 1, trainProportion = 0.8) {
  n = length(Y)
  nTrain = ceiling(n * trainProportion)
  iis = sample(n)
  iiTrain = iis[1:nTrain]
  iiTest = iis[(nTrain+1):n]
  
  FS1CV_OnSplits(Y, X, iiTrain, iiTest)
}



# Compute training fits and test RSSs for a given Train/Test split
FS1CV_OnSplits = function(Y, X, iiTrain, iiTest) {
  dataset = data.frame(Y, X)
  
  # Greedily choose X_i based on training data
  # whichX = which.max(abs(cor(cbind(Y[iiTrain], X[iiTrain, ])))[1, -1])
  whichX = which.max(abs(cor(Y[iiTrain], X[iiTrain, ])))
  
  # Fit all 3 models on training data
  trainFit0 = lm(Y ~ 1, dataset[iiTrain, ])
  trainFit1Fixed = lm(Y ~ ., dataset[iiTrain, c(1, 2)])
  trainFit1Greedy = lm(Y ~ ., dataset[iiTrain, c(1, whichX+1)])
  
  # Compute test error
  testRss0 = L2(Y[iiTest] - predict(trainFit0, newdata = X[iiTest, ]))
  testRss1Fixed = L2(Y[iiTest] - predict(trainFit1Fixed, newdata = X[iiTest, ]))
  testRss1Greedy = L2(Y[iiTest] - predict(trainFit1Greedy, newdata = X[iiTest, ]))
  
  # Return test errors to be compared
  c(testRss1Fixed = testRss1Fixed,
    testRss0 = testRss0,
    testRss1Greedy = testRss1Greedy)
}



# Generate data and run sample split
HighSnrHighP = function(n = 30, p = 2, B0 = 10,
                        verbose = FALSE, trainProportion = 0.8){
  data = genData(n, p, B0)

  if(verbose) {
    pairs(cbind(data$Y, data$X))
    print(cor(cbind(data$Y, data$X)))
    cat("\n")
  }
  
  testRSSs = FS1CV(data$Y, data$X, nFolds = 1, trainProportion = trainProportion)
  return(testRSSs)
}



#### Run it many times to confirm it's a trend: ####
# Very common to have ZeroBeatsGreedy (so CV works well with higher p);
# uncommon but not very rare to have FixedBeatsZero (so CV can fail with low p)
manyReps = function(nrReps = 500, n = 30, p = 2, B0 = 10, trainProportion = 0.8) {
  out = replicate(nrReps, { # 500 reps take ~10 sec at p=2
    oneRep = HighSnrHighP(n, p, B0, trainProportion = trainProportion)
    c(FixedBeatsZero = unname(oneRep["testRss1Fixed"] < oneRep["testRss0"]),
      ZeroBeatsGreedy = unname(oneRep["testRss0"] < oneRep["testRss1Greedy"]))
  })
  rowMeans(out)
}



#### Evaluate multi-split 80/20 sample-splitting: ####

# Rewrite HighSnrHighP() to allow multiple sample-splitting.
# Generate one dataset,
# and a list of M 80/20 splits (iTrain, iTest).
# For each split, evaluate testRss0, testRss1Fixed, and testRss1Greedy.
# Decide whether Fixed would choose K=0 or K=1; same for Greedy.
# Across splits, see whether K=0 or K=1 got more votes;
# if there's a tie, choose lower (K=0) for parsimony.
# (In later generalizations of this code, choose lowest modal K).
#
# Finally, wrap it inside replicate()
# and count how often Fixed chose K=0 correctly,
# vs how often Greedy chose K=0 correctly.
# 
# What do I expect to see?
# Each alg's rate should be higher than it was with a single sample-split and no voting;
# also, gap between Fixed and Greedy should be smaller;
# and I bet the effect of p will be less extreme.



# Generate data and run sample split
HighSnrHighP_MultiSplit = function(n = 30, p = 2, B0 = 10,
                                   nrSplits = 5, trainProportion = 0.8){
  data = genData(n, p, B0)

  out = replicate(nrSplits, {
    oneSplit = FS1CV(data$Y, data$X, nFolds = 1, trainProportion)
    c(KhatFixed = as.integer(unname(oneSplit["testRss0"] > oneSplit["testRss1Fixed"])),
      KhatGreedy = as.integer(unname(oneSplit["testRss0"] > oneSplit["testRss1Greedy"])))
  })
  # Find the modes: table() each row and see whether 0 or 1 was chosen more often
  c(KhatFixed = as.integer(names(which.max(table(out[1,])))),
    KhatGreedy = as.integer(names(which.max(table(out[2,])))))
}



manyReps_MultiSplit = function(nrReps = 500, n = 30, p = 2, B0 = 10,
                               nrSplits = 5, trainProportion = 0.8) {
  out = replicate(nrReps, { # 500 reps take ~40 sec at p=2
    HighSnrHighP_MultiSplit(n, p, B0, nrSplits, trainProportion)
  })
  # rowMeans() gives % of time 1 was chosen (incorrectly),
  # so take 1-rowMeans() to get % of time 0 was chosen (correctly)
  wins = 1 - rowMeans(out)
  # Note that "correct" here means correct model size K=0,
  # NOT necessarily that the correct model terms are chosen
  # (well, here it's the same thing, but could not be if true K>0)
  names(wins) = c("PropCorrect_Fixed", "PropCorrect_Greedy")
  wins
}



# Generate data and run sample split
MultiSplit_Freq0 = function(n = 30, p = 2, B0 = 10,
                                   nrSplits = 5, trainProportion = 0.5){
  data = genData(n, p, B0)
  
  out = replicate(nrSplits, {
    oneSplit = FS1CV(data$Y, data$X, nFolds = 1, trainProportion)
    # Did we correctly choose K=0?
    c(KhatFixed = as.integer(unname(oneSplit["testRss0"] <= oneSplit["testRss1Fixed"])),
      KhatGreedy = as.integer(unname(oneSplit["testRss0"] <= oneSplit["testRss1Greedy"])))
  })
  # Find the means: frequency of correctly choosing K=0 under each procedure
  c(Freq0Fixed = mean(out[1,]),
    Freq0Greedy = mean(out[2,]))
}



# Generate many datasets and do many sample splits on each
manyReps_MultiSplit_Freq0 = function(nrReps = 500, n = 30, p = 2, B0 = 10,
                               nrSplits = 7, trainProportion = 0.5) {
  out = replicate(nrReps, { # 500 reps take ~60 sec at p=2, nrSplits = 7
    MultiSplit_Freq0(n, p, B0, nrSplits, trainProportion)
  })
  # Don't summarize -- we want to plot histograms etc
  out
}

