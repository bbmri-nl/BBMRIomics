.glmnetFit <- function(y, X, trainid, testid, family, cv.opt, type, alpha){

    Xtrain <- X[trainid,]
    ytrain <- y[trainid]

    Xtest <-  X[testid,]
    ytest <- y[testid]

    fit <- cv.glmnet(Xtrain, ytrain, alpha = alpha, family = family)
    test <- predict(fit, Xtest, s = cv.opt, type = type)
    predicted <- predict(fit, X, s = cv.opt, type = type)

    list(test=as.vector(test), predicted=as.vector(predicted))
}

.gbmFit <- function(y, X, trainid, testid, family, type, n.trees, interaction.depth, cv.folds, shrinkage){

    Xtrain <- X[trainid,]
    ytrain <- y[trainid]

    Xtest <-  X[testid,]
    ytest <- y[testid]

    fit <- gbm(y~., data=data.frame(y=ytrain, x=Xtrain), distribution = family, n.trees=n.trees, interaction.depth=interaction.depth, cv.folds=cv.folds, shrinkage=shrinkage)

    test <- predict(fit, data.frame(y=ytest, x=Xtest), type = "response", n.trees = n.trees)
    if(type == "class")
        test <- apply(test, 1, function(x) colnames(test)[which.max(x)])

    predicted <- predict(fit, data.frame(y=y, x=X), type = "response", n.trees=n.trees)
    if(type == "class")
        predicted <- apply(predicted, 1, function(x) colnames(predicted)[which.max(x)])

    list(test=as.vector(test), predicted=as.vector(predicted))
}

.plsFit <- function(y, X, trainid, testid, ncomp=50, cv.opt, type){

    Xtrain <- X[trainid,]
    ytrain <- y[trainid]

    Xtest <-  X[testid,]
    ytest <- y[testid]

    if(type != "class") {
        fit <- plsr(y~., ncomp=min(nrow(X), ncomp), data=data.frame(y=ytrain, x=Xtrain), validation="CV")
        ncomp.onesigma <- selectNcomp(fit, method = cv.opt, plot = FALSE)
        test <- predict(fit, ncomp = ncomp.onesigma, newdata = data.frame(y=ytest, x=Xtest), type = type)
        predicted <- predict(fit, ncomp = ncomp.onesigma, newdata = data.frame(y=y, x=X), type=type)
    }
    else {
        fit <- cppls(y~., ncomp=min(nrow(X), ncomp), data=data.frame(y=ytrain, x=Xtrain), validation="CV")
        ncomp.onesigma <- selectNcomp(fit, method = cv.opt, plot = FALSE)
        test <- predict(fit, ncomp = ncomp.onesigma, newdata = data.frame(y=ytest, x=Xtest), type = "score")
        predicted <- predict(fit, ncomp = ncomp.onesigma, newdata = data.frame(y=y, x=X), type = "score")
    }

    list(test=as.vector(test), predicted=as.vector(predicted))
}

##' build prediction models for phenotype based on given matrix of features
##'
##' all penalized regression models e.g. n << p
##' @title predict phenotypes
##' @param phenotype vector with phenotype either numeric of categorical
##' @param features feature matrix
##' @param train.frac fraction fo the data that is used for training
##' @param methods regression methods used
##' @param verbose TRUE/FALSE
##' @param ... optional arguments
##' @importFrom glmnet cv.glmnet
##' @importFrom pls plsr selectNcomp
##' @importFrom gbm gbm 
##' @return validation and prediction
##' @author mvaniterson
phenotyping <- function(phenotype, features, train.frac=2/3, methods = c("ridge", "elastic-net", "lasso", "gradient-boosting", "partial-least-squares"), verbose, ...){

    if(length(phenotype) != ncol(features))
        stop("Number of phenotypes is not equal to the number of columns of the features!")

    methods <- match.arg(methods,
                         choices =  c("ridge", "elastic-net", "lasso", "gradient-boosting", "partial-least-squares"),
                         several.ok = TRUE)

    if(is.numeric(phenotype)) {
        message("Phenotyping on continuous phenotype ...")

        trainid <- sample(which(!is.na(phenotype)), size=floor(train.frac*sum(!is.na(phenotype))))
        family <- "gaussian"
        type <- "response"
        .performance <- .performanceContinuous

    } else if( nlevels(factor(phenotype)) == 2 ){
        message("Phenotyping on dichotomous phenotype ...")

        sp <- split(which(!is.na(phenotype)), phenotype[!is.na(phenotype)]) ##stratified sampling of test and train samples
        trainid <- unlist(lapply(sp, function(x) sample(x, length(x)*train.frac)), use.names=FALSE)
        family <- "bernoulli"
        type <- "class"
        .performance <- .performanceCategorical
    }
    else {
        message("Phenotyping on categorical phenotype ...")

        sp <- split(which(!is.na(phenotype)), phenotype[!is.na(phenotype)]) ##stratified sampling of test and train samples
        trainid <- unlist(lapply(sp, function(x) sample(x, length(x)*train.frac)), use.names=FALSE)
        family <- "multinomial"
        type <- "class"
        .performance <- .performanceCategorical
    }

    message("Using ", 100*round(train.frac,2), "% of data of train and ", 100*round(1 - train.frac, 2), "% for testing ...")
    testid <- setdiff(1:length(phenotype), c(trainid, which(is.na(phenotype))))

    phenotyped <- lapply(methods, function(method) {
        message("fit ", method, " regression model ...")
        phenotyped <- switch(method,
                             "ridge" = .glmnetFit(phenotype, t(features), trainid, testid,
                                                  family = family, cv.opt = "lambda.min", type = type,
                                                  alpha=0),

                             "elastic-net" = .glmnetFit(phenotype, t(features), trainid, testid,
                                                        family = family, cv.opt = "lambda.min", type = type,
                                                        alpha=1/2),

                             "lasso" = .glmnetFit(phenotype, t(features), trainid, testid,
                                                  family = family, cv.opt = "lambda.min", type = type,
                                                  alpha=1),

                             "gradient-boosting" = .gbmFit(phenotype, t(features), trainid, testid,
                                                           family = family, type = type,
                                                           n.trees=200, interaction.depth=4, cv.folds=5, shrinkage=0.005),

                             "partial-least-squares" = .plsFit(phenotype, t(features), trainid, testid,
                                                               cv.opt = "onesigma", type = type, ncomp=50),

                             message("Method, ", method, ", not implemented!"))

        list(validation = .performance(phenotype[testid], phenotyped$test, verbose=verbose),
             predicted = phenotyped$predicted)
    })

    validation <- do.call("rbind", lapply(phenotyped, function(x) x$validation))
    predicted <- do.call("cbind", lapply(phenotyped, function(x) x$predicted))
    colnames(predicted) <- rownames(validation) <- methods

    invisible(list(predicted=predicted, validation=validation, trainid=trainid, testid=testid))
}

.performanceContinuous <- function(reported, predicted, verbose) {
    correlation <- cor(reported, predicted)
    meddiff <- median(abs(reported - predicted))
    performance <- c(correlation, meddiff)
    names(performance) <- c("correlation (Pearson)", "error (median abs. diff.)")
    performance
}

.performanceCategorical <- function(reported, predicted, verbose) {

    cm <- as.matrix(table(Reported=reported, Predicted=predicted))

    if(verbose){
        message("Confusion matrix:\n")
        print(cm)
        message("\n")
    }

    n <- sum(cm) # number of instances
    nc <- nrow(cm) # number of classes
    diag <- diag(cm) # number of correctly classified instances per class
    rowsums <- apply(cm, 1, sum) # number of instances per class
    colsums <- apply(cm, 2, sum) # number of predictions per class
    p <- rowsums / n # distribution of instances over the classes
    q <- colsums / n # distribution of instances over the predicted classes

    ##accuracy
    accuracy <- sum(diag) / n
    names(accuracy) <- "accuracy (Overall)"

    if(verbose)
        message(names(accuracy), ": ", round(as.numeric(accuracy), 2))

    ##per class prf
    recall <- diag / rowsums
    precision <- diag / colsums
    f1 <- 2 * precision * recall / (precision + recall)

    names(recall) <- paste0("recall (", names(recall), ")")
    names(precision) <- paste0("precision (", names(precision), ")")
    names(f1) <- paste0("f1 (", names(f1), ")")

    performance <- c(accuracy, recall, precision, f1)
    invisible(performance)
}


