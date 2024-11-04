ExtractData<-function (object, data, knots) 
{
    knt <- dat <- list()
    vecMat <- if (!is.list(object$xt) || is.null(object$xt$sumConv)) 
        TRUE
    else object$xt$sumConv
    for (i in 1:length(object$term)) {
        dat[[object$term[i]]] <- get.var(object$term[i], data, 
            vecMat = vecMat)
        knt[[object$term[i]]] <- get.var(object$term[i], knots, 
            vecMat = vecMat)
    }
    names(dat) <- object$term
    m <- length(object$term)
    if (!is.null(attr(dat[[1]], "matrix")) && vecMat) {
        n <- length(dat[[1]])
        X <- data.frame(dat)
        X <- uniquecombs(X)
        if (nrow(X) < n * 0.9) {
            for (i in 1:m) dat[[i]] <- X[, i]
            attr(dat, "index") <- attr(X, "index")
        }
    }
    if (object$by != "NA") {
        by <- get.var(object$by, data)
        if (!is.null(by)) {
            dat[[m + 1]] <- by
            names(dat)[m + 1] <- object$by
        }
    }
    return(list(data = dat, knots = knt))
}
