source("http://pastebin.com/raw.php?i=JVyTrYRD")  # source Devium


data(iris)
tmp.data <- iris[, -5]
tmp.group <- iris[, 5]  # species
tmp.y <- matrix(as.numeric(tmp.group), ncol = 1)  # make numeric matrix



train.test.index.main <- test.train.split(nrow(tmp.data), n = 1, strata = tmp.group, 
                                         split.type = "duplex", data = tmp.data)

train.id <- train.test.index.main == "train"

# partition data to get the trainning set
tmp.data <- tmp.data[train.id, ]
tmp.group <- tmp.group[train.id]
tmp.y <- tmp.y[train.id, ]

# the variables could be scaled now, or done internally in the model for
# each CV split (leave-one-out)
# scaled.data<-data.frame(scale(tmp.data,center=TRUE, scale=TRUE))
scaled.data <- tmp.data



mods <- make.OSC.PLS.model(tmp.y, pls.data = scaled.data, comp = 2, OSC.comp = 1, 
                           validation = "LOO", method = "oscorespls", cv.scale = TRUE, progress = FALSE)
# extract model
final <- get.OSC.model(obj = mods, OSC.comp = 1)
# view out-of-bag error for cross-validation splits
plot.OSC.results(mods, plot = "RMSEP", groups = tmp.group)


plot.OSC.results(mods, plot = "scores", groups = tmp.group)

plot.PLS.results(obj = final, plot = "scores", groups = tmp.group)



train.test.index = test.train.split(nrow(scaled.data), n = 100, strata = as.factor(tmp.y))  # strata controls if the species are sampled from equally
permuted.stats <- permute.OSC.PLS(data = scaled.data, y = as.matrix(tmp.y), 
                                  n = 50, ncomp = 2, osc.comp = 1, progress = FALSE, train.test.index = train.test.index)
# look how our model compares to random chance
OSC.validate.model(model = mods, perm = permuted.stats)


# reset data
scaled.data <- iris[, -5]
tmp.group <- iris[, 5]
tmp.y <- matrix(as.numeric(tmp.group), ncol = 1)

# make predictions for the test set
mods <- make.OSC.PLS.model(tmp.y, pls.data = scaled.data, comp = 2, OSC.comp = 1, 
                           validation = "LOO", method = "oscorespls", cv.scale = TRUE, progress = FALSE, 
                           train.test.index = train.test.index.main)

# get the true (actual) and predicted values round them to integers to
# represent discreet species labels
plot.data = data.frame(predicted = round(mods$predicted.Y[[2]][, 1], 0), actual = mods$test.y)
# note these are numeric but we would prefer to interpret classification
# of species a class
plot.data$predicted <- factor(plot.data$predicted, labels = levels(iris[, 5]), 
                              levels = 1:3)
plot.data$actual <- factor(plot.data$actual, labels = levels(iris[, 5]), levels = 1:3)

table(plot.data)


pairs(iris[, -5], pch = 21, bg = rainbow(nlevels(iris[, 5]), alpha = 0.75)[iris[, 
                                                                                5]], upper.panel = NULL, cex = 2)
par(xpd = TRUE)
legend(0.75, 1, levels(iris[, 5]), fill = rainbow(nlevels(iris[, 5]), alpha = 0.75), 
       bty = "n")
