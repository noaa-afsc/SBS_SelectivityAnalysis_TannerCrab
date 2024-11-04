PredictMat<-function (object, data, n = nrow(data)) 
{
    pm <- Predict.matrix3(object, data)
    qrc <- attr(object, "qrc")
    if (inherits(qrc, "sweepDrop")) {
        deriv <- if (is.null(object$deriv) || object$deriv == 
            0) 
            FALSE
        else TRUE
        if (!deriv && !is.null(object$margin)) 
            for (i in 1:length(object$margin)) if (!is.null(object$margin[[i]]$deriv) && 
                object$margin[[i]]$deriv != 0) 
                deriv <- TRUE
        if (!deriv) 
            pm$X <- pm$X[, -qrc[1], drop = FALSE] - matrix(qrc[-1], 
                nrow(pm$X), ncol(pm$X) - 1, byrow = TRUE)
        else pm$X <- pm$X[, -qrc[1], drop = FALSE]
    }
    if (!is.null(pm$ind) && length(pm$ind) != n) {
        if (is.null(attr(pm$X, "by.done")) && object$by != "NA") {
            by <- get.var(object$by, data)
            if (is.null(by)) 
                stop("Can't find by variable")
        }
        else by <- rep(1, length(pm$ind))
        q <- length(pm$ind)/n
        ind <- 0:(q - 1) * n
        offs <- attr(pm$X, "offset")
        if (!is.null(offs)) 
            offX <- rep(0, n)
        else offX <- NULL
        X <- matrix(0, n, ncol(pm$X))
        for (i in 1:n) {
            ind <- ind + 1
            X[i, ] <- colSums(by[ind] * pm$X[pm$ind[ind], , drop = FALSE])
            if (!is.null(offs)) {
                offX[i] <- sum(offs[pm$ind[ind]] * by[ind])
            }
        }
        offset <- offX
    }
    else {
        offset <- attr(pm$X, "offset")
        if (!is.null(pm$ind)) {
            X <- pm$X[pm$ind, , drop = FALSE]
            if (!is.null(offset)) 
                offset <- offset[pm$ind]
        }
        else X <- pm$X
        if (is.null(attr(pm$X, "by.done"))) {
            if (object$by != "NA") {
                by <- get.var(object$by, data)
                if (is.null(by)) 
                  stop("Can't find by variable")
                if (is.factor(by)) {
                  by.dum <- as.numeric(object$by.level == by)
                  X <- by.dum * X
                  if (!is.null(offset)) 
                    offset <- by.dum * offset
                }
                else {
                  if (length(by) != nrow(X)) 
                    stop("`by' variable must be same dimension as smooth arguments")
                  X <- as.numeric(by) * X
                  if (!is.null(offset)) 
                    offset <- as.numeric(by) * offset
                }
            }
        }
        rm(pm)
        attr(X, "by.done") <- NULL
        if (n != nrow(X)) {
            q <- nrow(X)/n
            ind <- 1:n
            Xs <- X[ind, ]
            if (!is.null(offset)) {
                get.off <- TRUE
                offs <- offset[ind]
            }
            else {
                get.off <- FALSE
                offs <- NULL
            }
            for (i in 2:q) {
                ind <- ind + n
                Xs <- Xs + X[ind, , drop = FALSE]
                if (get.off) 
                  offs <- offs + offset[ind]
            }
            offset <- offs
            X <- Xs
        }
    }
    if (!is.null(qrc)) {
        j <- attr(object, "nCons")
        if (j > 0) {
            k <- ncol(X)
            if (inherits(qrc, "qr")) {
                indi <- attr(object, "indi")
                if (is.null(indi)) {
                  if (sum(is.na(X))) {
                    ind <- !is.na(rowSums(X))
                    X1 <- t(qr.qty(qrc, t(X[ind, , drop = FALSE]))[(j + 
                      1):k, , drop = FALSE])
                    X <- matrix(NA, nrow(X), ncol(X1))
                    X[ind, ] <- X1
                  }
                  else {
                    X <- t(qr.qty(qrc, t(X))[(j + 1):k, , drop = FALSE])
                  }
                }
                else {
                  nx <- length(indi)
                  nc <- j
                  nz <- nx - nc
                  if (sum(is.na(X))) {
                    ind <- !is.na(rowSums(X))
                    if (nz > 0) 
                      X[ind, indi[1:nz]] <- t(qr.qty(qrc, t(X[ind, 
                        indi, drop = FALSE]))[(nc + 1):nx, ])
                    X <- X[, -indi[(nz + 1):nx], drop = FALSE]
                    X[!ind, ] <- NA
                  }
                  else {
                    if (nz > 0) 
                      X[, indi[1:nz]] <- t(qr.qty(qrc, t(X[, 
                        indi, drop = FALSE]))[(nc + 1):nx, , 
                        drop = FALSE])
                    X <- X[, -indi[(nz + 1):nx], drop = FALSE]
                  }
                }
            }
            else if (inherits(qrc, "sweepDrop")) {
            }
            else if (length(qrc) > 1) {
                m <- qrc[-c(1, length(qrc))]
                if (length(m) > 0) 
                  X <- t(XZKr(X, m))
            }
            else if (qrc > 0) {
                X <- X[, -qrc, drop = FALSE]
            }
            else if (qrc < 0) {
                X <- t(diff(t(X)))
            }
        }
    }
    if (!is.null(object$diagRP)) 
        X <- X %*% object$diagRP
    del.index <- attr(object, "del.index")
    if (!is.null(del.index)) 
        X <- X[, -del.index, drop = FALSE]
    attr(X, "offset") <- offset
    X
}
