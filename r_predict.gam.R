predict.gam<-function (object, newdata, type = "link", se.fit = FALSE, terms = NULL, 
    exclude = NULL, block.size = NULL, newdata.guaranteed = FALSE, 
    na.action = na.pass, unconditional = FALSE, iterms.type = NULL, 
    ...) 
{
    if (unconditional) {
        if (is.null(object$Vc)) 
            warning("Smoothness uncertainty corrected covariance not available")
        else object$Vp <- object$Vc
    }
    if (type != "link" && type != "terms" && type != "iterms" && 
        type != "response" && type != "lpmatrix" && type != "newdata") {
        warning("Unknown type, reset to terms.")
        type <- "terms"
    }
    if (!inherits(object, "gam")) 
        stop("predict.gam can only be used to predict from gam objects")
    if (missing(newdata)) 
        na.act <- object$na.action
    else {
        if (is.null(na.action)) 
            na.act <- NULL
        else {
            na.txt <- if (is.character(na.action) || is.function(na.action)) 
                mgcv:::get.na.action(na.action)
            else "na.pass"
            if (na.txt == "na.pass") 
                na.act <- "na.exclude"
            else if (na.txt == "na.exclude") 
                na.act <- "na.omit"
            else na.act <- na.action
        }
    }
    nd.is.mf <- FALSE
    yname <- attr(attr(object$terms, "dataClasses"), "names")[attr(object$terms, 
        "response")]
    if (newdata.guaranteed == FALSE) {
        if (missing(newdata)) {
            newdata <- object$model
            new.data.ok <- FALSE
            nd.is.mf <- TRUE
            response <- newdata[[yname]]
        }
        else {
            new.data.ok <- TRUE
            if (is.data.frame(newdata) && !is.null(attr(newdata, 
                "terms"))) {
                if (sum(!(names(object$model) %in% names(newdata)))) 
                  stop("newdata is a model.frame: it should contain all required variables\n")
                nd.is.mf <- TRUE
                response <- newdata[[yname]]
            }
            else {
                resp <- mgcv::get.var(yname, newdata, FALSE)
                naresp <- FALSE
                if (!is.null(object$family$predict) && !is.null(resp)) {
                  if (!is.null(object$pred.formula)) 
                    object$pred.formula <- attr(object$pred.formula, 
                      "full")
                  response <- TRUE
                  Terms <- terms(object)
                  if (is.matrix(resp)) {
                    if (sum(is.na(rowSums(resp))) > 0) 
                      stop("no NAs allowed in response data for this model")
                  }
                  else {
                    if (sum(is.na(resp)) > 0) {
                      naresp <- TRUE
                      rar <- range(resp, na.rm = TRUE)
                      thresh <- rar[1] * 1.01 - rar[2] * 0.01
                      resp[is.na(resp)] <- thresh
                      newdata[[yname]] <- thresh
                    }
                  }
                }
                else {
                  response <- FALSE
                  Terms <- delete.response(terms(object))
                }
                allNames <- if (is.null(object$pred.formula)) 
                  all.vars(Terms)
                else all.vars(object$pred.formula)
                if (length(allNames) > 0) {
                  ff <- if (is.null(object$pred.formula)) 
                    reformulate(allNames)
                  else object$pred.formula
                  if (sum(!(allNames %in% names(newdata)))) {
                    warning("not all required variables have been supplied in  newdata!\n")
                  }
                  newdata <- eval(model.frame(ff, data = newdata, 
                    na.action = na.act), parent.frame())
                  if (naresp) 
                    newdata[[yname]][newdata[[yname]] <= thresh] <- NA
                }
                na.act <- attr(newdata, "na.action")
                response <- if (response) 
                  get.var(yname, newdata, FALSE)
                else NULL
            }
        }
    }
    else {
        na.act <- NULL
        new.data.ok = TRUE
        if (!is.null(attr(newdata, "terms"))) 
            nd.is.mf <- TRUE
        response <- get.var(yname, newdata, FALSE)
    }
    if (new.data.ok) {
        nn <- names(newdata)
        mn <- colnames(object$model)
        for (i in 1:length(newdata)) if (nn[i] %in% mn && is.factor(object$model[, 
            nn[i]])) {
            levm <- levels(object$model[, nn[i]])
            levn <- if (any(is.na(levm))) 
                levels(factor(newdata[[i]], exclude = NULL))
            else levels(factor(newdata[[i]]))
            if (sum(!levn %in% levm) > 0) {
                msg <- paste("factor levels", paste(levn[!levn %in% 
                  levm], collapse = ", "), "not in original fit", 
                  collapse = "")
                warning(msg)
            }
            if (is.matrix(newdata[[i]])) {
                dum <- factor(newdata[[i]], levels = levm, exclude = NULL)
                dim(dum) <- dim(newdata[[i]])
                newdata[[i]] <- dum
            }
            else newdata[[i]] <- factor(newdata[[i]], levels = levm, 
                exclude = NULL)
        }
        if (type == "newdata") 
            return(newdata)
        if (length(newdata) == 1) 
            newdata[[2]] <- newdata[[1]]
        if (is.null(dim(newdata[[1]]))) 
            np <- length(newdata[[1]])
        else np <- dim(newdata[[1]])[1]
        nb <- length(object$coefficients)
        if (is.null(block.size)) 
            block.size <- 1000
        if (block.size < 1) 
            block.size <- np
    }
    else {
        np <- nrow(object$model)
        nb <- length(object$coefficients)
    }
    if (type == "lpmatrix") 
        block.size <- NULL
    if (is.null(block.size)) {
        n.blocks <- 1
        b.size <- array(np, 1)
    }
    else {
        n.blocks <- np%/%block.size
        b.size <- rep(block.size, n.blocks)
        last.block <- np - sum(b.size)
        if (last.block > 0) {
            n.blocks <- n.blocks + 1
            b.size[n.blocks] <- last.block
        }
    }
    lpi <- if (is.list(object$formula)) 
        attr(object$formula, "lpi")
    else NULL
    nlp <- if (is.null(lpi)) 
        1
    else length(lpi)
    n.smooth <- length(object$smooth)
    if (type == "lpmatrix") {
        H <- matrix(0, np, nb)
    }
    else if (type == "terms" || type == "iterms") {
        term.labels <- attr(object$pterms, "term.labels")
        para.only <- attr(object, "para.only")
        if (is.null(para.only)) 
            para.only <- FALSE
        n.pterms <- length(term.labels)
        fit <- array(0, c(np, n.pterms + as.numeric(para.only == 
            0) * n.smooth))
        if (se.fit) 
            se <- fit
        ColNames <- term.labels
    }
    else {
        fit <- if (nlp > 1) 
            matrix(0, np, nlp)
        else array(0, np)
        if (se.fit) 
            se <- fit
        fit1 <- NULL
    }
    stop <- 0
    if (is.list(object$pterms)) {
        if (type == "iterms") {
            warning("type iterms not available for multiple predictor cases")
            type <- "terms"
        }
        pstart <- attr(object$nsdf, "pstart")
        pind <- rep(0, 0)
        Terms <- list()
        pterms <- object$pterms
        for (i in 1:length(object$nsdf)) {
            Terms[[i]] <- delete.response(object$pterms[[i]])
            if (object$nsdf[i] > 0) 
                pind <- c(pind, pstart[i] - 1 + 1:object$nsdf[i])
        }
    }
    else {
        Terms <- list(delete.response(object$pterms))
        pterms <- list(object$pterms)
        pstart <- 1
        pind <- 1:object$nsdf
    }
    drop.intercept <- object$family$drop.intercept
    if (is.null(drop.intercept)) {
        drop.intercept <- rep(FALSE, length(Terms))
    }
    else {
        for (i in 1:length(Terms)) {
            if (drop.intercept[i] == TRUE) 
                attr(Terms[[i]], "intercept") <- 1
        }
    }
    drop.ind <- attr(object$nsdf, "drop.ind")
    s.offset <- NULL
    any.soff <- FALSE
    if (n.blocks > 0) 
        for (b in 1:n.blocks) {
            start <- stop + 1
            stop <- start + b.size[b] - 1
            if (n.blocks == 1) 
                data <- newdata
            else data <- newdata[start:stop, ]
            X <- matrix(0, b.size[b], nb + length(drop.ind))
            Xoff <- matrix(0, b.size[b], n.smooth)
            offs <- list()
            for (i in 1:length(Terms)) {
                if (new.data.ok) {
                  if (nd.is.mf) 
                    mf <- model.frame(data, xlev = object$xlevels)
                  else {
                    mf <- model.frame(Terms[[i]], data, xlev = object$xlevels)
                    if (!is.null(cl <- attr(pterms[[i]], "dataClasses"))) 
                      .checkMFClasses(cl, mf)
                  }
                  oc <- if (length(object$contrasts) == 0) 
                    object$contrasts
                  else object$contrasts[names(object$contrasts) %in% 
                    attr(Terms[[i]], "term.labels")]
                  Xp <- model.matrix(Terms[[i]], mf, contrasts = oc)
                }
                else {
                  Xp <- model.matrix(Terms[[i]], object$model)
                  mf <- newdata
                }
                if (!is.null(terms) || !is.null(exclude)) {
                  assign <- attr(Xp, "assign")
                  if (min(assign) == 0 && ("(Intercept)" %in% 
                    exclude || (!is.null(terms) && !"(Intercept)" %in% 
                    terms))) 
                    Xp[, which(assign == 0)] <- 0
                  tlab <- attr(Terms[[i]], "term.labels")
                  ii <- which(assign %in% which(tlab %in% exclude))
                  if (length(ii)) 
                    Xp[, ii] <- 0
                  if (!is.null(terms)) {
                    ii <- which(assign %in% which(!tlab %in% 
                      terms))
                    if (length(ii)) 
                      Xp[, ii] <- 0
                  }
                }
                offi <- attr(Terms[[i]], "offset")
                if (is.null(offi)) 
                  offs[[i]] <- 0
                else {
                  offs[[i]] <- mf[[names(attr(Terms[[i]], "dataClasses"))[offi + 
                    1]]]
                }
                if (drop.intercept[i]) {
                  xat <- attributes(Xp)
                  ind <- xat$assign > 0
                  Xp <- Xp[, xat$assign > 0, drop = FALSE]
                  xat$assign <- xat$assign[ind]
                  xat$dimnames[[2]] <- xat$dimnames[[2]][ind]
                  xat$dim[2] <- xat$dim[2] - 1
                  attributes(Xp) <- xat
                }
                if (object$nsdf[i] > 0) 
                  X[, pstart[i] - 1 + 1:object$nsdf[i]] <- Xp
            }
            if (!is.null(drop.ind)) 
                X <- X[, -drop.ind]
            if (n.smooth) 
                for (k in 1:n.smooth) {
                  klab <- object$smooth[[k]]$label
                  if ((is.null(terms) || (klab %in% terms)) && 
                    (is.null(exclude) || !(klab %in% exclude))) {
                    Xfrag <- PredictMat(object$smooth[[k]], data)
                    X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
                    Xfrag.off <- attr(Xfrag, "offset")
                    if (!is.null(Xfrag.off)) {
                      Xoff[, k] <- Xfrag.off
                      any.soff <- TRUE
                    }
                  }
                  if (type == "terms" || type == "iterms") 
                    ColNames[n.pterms + k] <- klab
                }
            if (!is.null(object$Xcentre)) {
                X <- sweep(X, 2, object$Xcentre)
            }
            if (type == "lpmatrix") {
                H[start:stop, ] <- X
                if (any.soff) 
                  s.offset <- rbind(s.offset, Xoff)
            }
            else if (type == "terms" || type == "iterms") {
                lass <- if (is.list(object$assign)) 
                  object$assign
                else list(object$assign)
                k <- 0
                for (j in 1:length(lass)) if (length(lass[[j]])) {
                  ind <- 1:length(lass[[j]])
                  nptj <- max(lass[[j]])
                  if (nptj > 0) 
                    for (i in 1:nptj) {
                      k <- k + 1
                      ii <- ind[lass[[j]] == i] + pstart[j] - 
                        1
                      fit[start:stop, k] <- X[, ii, drop = FALSE] %*% 
                        object$coefficients[ii]
                      if (se.fit) 
                        se[start:stop, k] <- sqrt(pmax(0, rowSums((X[, 
                          ii, drop = FALSE] %*% object$Vp[ii, 
                          ii]) * X[, ii, drop = FALSE])))
                    }
                }
                if (n.smooth && !para.only) {
                  for (k in 1:n.smooth) {
                    first <- object$smooth[[k]]$first.para
                    last <- object$smooth[[k]]$last.para
                    fit[start:stop, n.pterms + k] <- X[, first:last, 
                      drop = FALSE] %*% object$coefficients[first:last] + 
                      Xoff[, k]
                    if (se.fit) {
                      if (type == "iterms" && attr(object$smooth[[k]], 
                        "nCons") > 0) {
                        if (length(object$cmX) < ncol(X)) 
                          object$cmX <- c(object$cmX, rep(0, 
                            ncol(X) - length(object$cmX)))
                        if (!is.null(iterms.type) && iterms.type == 
                          2) 
                          object$cmX[-(1:object$nsdf)] <- 0
                        X1 <- matrix(object$cmX, nrow(X), ncol(X), 
                          byrow = TRUE)
                        meanL1 <- object$smooth[[k]]$meanL1
                        if (!is.null(meanL1)) 
                          X1 <- X1/meanL1
                        X1[, first:last] <- X[, first:last]
                        se[start:stop, n.pterms + k] <- sqrt(pmax(0, 
                          rowSums((X1 %*% object$Vp) * X1)))
                      }
                      else se[start:stop, n.pterms + k] <- sqrt(pmax(0, 
                        rowSums((X[, first:last, drop = FALSE] %*% 
                          object$Vp[first:last, first:last, drop = FALSE]) * 
                          X[, first:last, drop = FALSE])))
                    }
                  }
                  colnames(fit) <- ColNames
                  if (se.fit) 
                    colnames(se) <- ColNames
                }
                else {
                  if (para.only && is.list(object$pterms)) {
                    term.labels <- unlist(lapply(object$pterms, 
                      attr, "term.labels"))
                  }
                  colnames(fit) <- term.labels
                  if (se.fit) 
                    colnames(se) <- term.labels
                  if (para.only) {
                    order <- if (is.list(object$pterms)) 
                      unlist(lapply(object$pterms, attr, "order"))
                    else attr(object$pterms, "order")
                    term.labels <- term.labels[order == 1]
                    fit <- fit[, order == 1, drop = FALSE]
                    colnames(fit) <- term.labels
                    if (se.fit) {
                      se <- se[, order == 1, drop = FALSE]
                      colnames(se) <- term.labels
                    }
                  }
                }
            }
            else {
                fam <- object$family
                if (!is.null(object$family$setInd)) 
                  object$family$setInd(start:stop)
                k <- attr(attr(object$model, "terms"), "offset")
                if (nlp > 1) {
                  if (is.null(fam$predict) || type == "link") {
                    off.ind <- (1:n.smooth)[as.logical(colSums(abs(Xoff)))]
                    for (j in 1:nlp) {
                      ind <- lpi[[j]]
                      fit[start:stop, j] <- X[, ind, drop = FALSE] %*% 
                        object$coefficients[ind] + offs[[j]]
                      if (length(off.ind)) 
                        for (i in off.ind) {
                          if (object$smooth[[i]]$first.para %in% 
                            ind) 
                            fit[start:stop, j] <- fit[start:stop, 
                              j] + Xoff[, i]
                        }
                      if (se.fit) 
                        se[start:stop, j] <- sqrt(pmax(0, rowSums((X[, 
                          ind, drop = FALSE] %*% object$Vp[ind, 
                          ind, drop = FALSE]) * X[, ind, drop = FALSE])))
                      if (type == "response") {
                        linfo <- object$family$linfo[[j]]
                        if (se.fit) 
                          se[start:stop, j] <- se[start:stop, 
                            j] * abs(linfo$mu.eta(fit[start:stop, 
                            j]))
                        fit[start:stop, j] <- linfo$linkinv(fit[start:stop, 
                          j])
                      }
                    }
                  }
                  else {
                    attr(X, "lpi") <- lpi
                    ffv <- fam$predict(fam, se.fit, y = if (is.matrix(response)) 
                      response[start:stop, ]
                    else response[start:stop], X = X, beta = object$coefficients, 
                      off = offs, Vb = object$Vp)
                    if (is.matrix(fit) && !is.matrix(ffv[[1]])) {
                      fit <- fit[, 1]
                      if (se.fit) 
                        se <- se[, 1]
                    }
                    if (is.matrix(ffv[[1]]) && (!is.matrix(fit) || 
                      ncol(ffv[[1]]) != ncol(fit))) {
                      fit <- matrix(0, np, ncol(ffv[[1]]))
                      if (se.fit) 
                        se <- fit
                    }
                    if (is.matrix(fit)) {
                      fit[start:stop, ] <- ffv[[1]]
                      if (se.fit) 
                        se[start:stop, ] <- ffv[[2]]
                    }
                    else {
                      fit[start:stop] <- ffv[[1]]
                      if (se.fit) 
                        se[start:stop] <- ffv[[2]]
                    }
                  }
                }
                else {
                  offs <- if (is.null(k)) 
                    rowSums(Xoff)
                  else rowSums(Xoff) + model.offset(mf)
                  fit[start:stop] <- X %*% object$coefficients + 
                    offs
                  if (se.fit) 
                    se[start:stop] <- sqrt(pmax(0, rowSums((X %*% 
                      object$Vp) * X)))
                  if (type == "response") {
                    linkinv <- fam$linkinv
                    if (is.null(fam$predict)) {
                      dmu.deta <- fam$mu.eta
                      if (se.fit) 
                        se[start:stop] <- se[start:stop] * abs(dmu.deta(fit[start:stop]))
                      fit[start:stop] <- linkinv(fit[start:stop])
                    }
                    else {
                      ffv <- fam$predict(fam, se.fit, y = if (is.matrix(response)) 
                        response[start:stop, ]
                      else response[start:stop], X = X, beta = object$coefficients, 
                        off = offs, Vb = object$Vp)
                      if (is.null(fit1) && is.matrix(ffv[[1]])) {
                        fit1 <- matrix(0, np, ncol(ffv[[1]]))
                        if (se.fit) 
                          se1 <- fit1
                      }
                      if (is.null(fit1)) {
                        fit[start:stop] <- ffv[[1]]
                        if (se.fit) 
                          se[start:stop] <- ffv[[2]]
                      }
                      else {
                        fit1[start:stop, ] <- ffv[[1]]
                        if (se.fit) 
                          se1[start:stop, ] <- ffv[[2]]
                      }
                    }
                  }
                }
            }
            rm(X)
        }
    if (!is.null(object$family$setInd)) 
        object$family$setInd(NULL)
    if ((type == "terms" || type == "iterms") && (!is.null(terms) || 
        !is.null(exclude))) {
        cnames <- colnames(fit)
        if (!is.null(terms)) {
            if (sum(!(terms %in% cnames))) 
                warning("non-existent terms requested - ignoring")
            else {
                fit <- fit[, terms, drop = FALSE]
                if (se.fit) {
                  se <- se[, terms, drop = FALSE]
                }
            }
        }
        if (!is.null(exclude)) {
            if (sum(!(exclude %in% cnames))) 
                warning("non-existent exclude terms requested - ignoring")
            else {
                exclude <- which(cnames %in% exclude)
                fit <- fit[, -exclude, drop = FALSE]
                if (se.fit) {
                  se <- se[, -exclude, drop = FALSE]
                }
            }
        }
    }
    if (type == "response" && !is.null(fit1)) {
        fit <- fit1
        if (se.fit) 
            se <- se1
    }
    rn <- rownames(newdata)
    if (type == "lpmatrix") {
        colnames(H) <- names(object$coefficients)
        rownames(H) <- rn
        if (!is.null(s.offset)) {
            s.offset <- napredict(na.act, s.offset)
            attr(H, "offset") <- s.offset
        }
        if (!is.null(offs)) {
            offs <- offs[1:nlp]
            for (i in 1:nlp) offs[[i]] <- napredict(na.act, offs[[i]])
            attr(H, "model.offset") <- if (nlp == 1) 
                offs[[1]]
            else offs
        }
        H <- napredict(na.act, H)
        if (length(object$nsdf) > 1) {
            attr(H, "lpi") <- lpi
        }
    }
    else {
        if (se.fit) {
            if (is.null(nrow(fit))) {
                names(fit) <- rn
                names(se) <- rn
                fit <- napredict(na.act, fit)
                se <- napredict(na.act, se)
            }
            else {
                rownames(fit) <- rn
                rownames(se) <- rn
                fit <- napredict(na.act, fit)
                se <- napredict(na.act, se)
            }
            H <- list(fit = fit, se.fit = se)
        }
        else {
            H <- fit
            if (is.null(nrow(H))) 
                names(H) <- rn
            else rownames(H) <- rn
            H <- napredict(na.act, H)
        }
    }
    if ((type == "terms" || type == "iterms") && attr(object$terms, 
        "intercept") == 1) 
        attr(H, "constant") <- object$coefficients[1]
    H
}
