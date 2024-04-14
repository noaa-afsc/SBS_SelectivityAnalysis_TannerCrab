myGam<-function (formula, family = gaussian(), data = list(), weights = NULL, 
    subset = NULL, na.action, offset = NULL, method = "GCV.Cp", 
    optimizer = c("outer", "newton"), control = list(), scale = 0, 
    select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, 
    gamma = 1, fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, 
    drop.unused.levels = TRUE, drop.intercept = NULL, nei = NULL, 
    discrete = FALSE, ...) 
{
    control <- do.call("gam.control", control)
    if (is.null(G) && discrete) {
        cl <- match.call()
        cl[[1]] <- quote(bam)
        cl$fit = FALSE
        G <- eval(cl, parent.frame())
    }
    if (is.null(G)) {
        gp <- interpret.gam(formula)
        cl <- match.call()
        mf <- match.call(expand.dots = FALSE)
        mf$formula <- gp$fake.formula
        mf$family <- mf$control <- mf$scale <- mf$knots <- mf$sp <- mf$min.sp <- mf$H <- mf$select <- mf$drop.intercept <- mf$nei <- mf$gamma <- mf$method <- mf$fit <- mf$paraPen <- mf$G <- mf$optimizer <- mf$in.out <- mf$discrete <- mf$... <- NULL
        mf$drop.unused.levels <- drop.unused.levels
        mf[[1]] <- quote(myModel.Frame.Default);#quote(stats::model.frame)
        pmf <- mf
        mf <- eval(mf, parent.frame())
        if (nrow(mf) < 2) 
            stop("Not enough (non-NA) data to do anything meaningful")
        terms <- attr(mf, "terms")
        if (!is.null(nei)) {
            k <- attr(mf, "na.action")
            if (!is.null(k)) {
                nei <- nanei(nei, as.numeric(k))
            }
        }
        vars <- all_vars1(gp$fake.formula[-2])
        inp <- parse(text = paste("list(", paste(vars, collapse = ","), 
            ")"))
        if (!is.list(data) && !is.data.frame(data)) 
            data <- as.data.frame(data)
        dl <- eval(inp, data, parent.frame())
        names(dl) <- vars
        var.summary <- variable.summary(gp$pf, dl, nrow(mf))
        rm(dl)
        if (is.list(formula)) {
            environment(formula) <- environment(formula[[1]])
            pterms <- list()
            tlab <- rep("", 0)
            for (i in 1:length(formula)) {
                pmf$formula <- gp[[i]]$pf
                pterms[[i]] <- attr(eval(pmf, parent.frame()), 
                  "terms")
                tlabi <- attr(pterms[[i]], "term.labels")
                if (i > 1 && length(tlabi) > 0) 
                  tlabi <- paste(tlabi, i - 1, sep = ".")
                tlab <- c(tlab, tlabi)
            }
            attr(pterms, "term.labels") <- tlab
        }
        else {
            pmf$formula <- gp$pf
            pmf <- eval(pmf, parent.frame())
            pterms <- attr(pmf, "terms")
        }
        if (is.character(family)) 
            family <- eval(parse(text = family))
        if (is.function(family)) 
            family <- family()
        if (is.null(family$family)) 
            stop("family not recognized")
        if (family$family[1] == "gaussian" && family$link == 
            "identity") 
            am <- TRUE
        else am <- FALSE
        if (!control$keepData) 
            rm(data)
        if (is.null(family$drop.intercept)) {
            lengthf <- if (is.list(formula)) 
                length(formula)
            else 1
            if (is.null(drop.intercept)) 
                drop.intercept <- rep(FALSE, lengthf)
            else {
                drop.intercept <- rep(drop.intercept, length = lengthf)
                if (sum(drop.intercept)) 
                  family$drop.intercept <- drop.intercept
            }
        }
        else drop.intercept <- as.logical(family$drop.intercept)
        if (inherits(family, "general.family") && !is.null(family$presetup)) 
            eval(family$presetup)
        gsname <- if (is.list(formula)) 
            "gam.setup.list"
        else "gam.setup"
        G <- do.call(gsname, list(formula = gp, pterms = pterms, 
            data = mf, knots = knots, sp = sp, min.sp = min.sp, 
            H = H, absorb.cons = TRUE, sparse.cons = 0, select = select, 
            idLinksBases = control$idLinksBases, scale.penalty = control$scalePenalty, 
            paraPen = paraPen, drop.intercept = drop.intercept))
        G$var.summary <- var.summary
        G$family <- family
        if ((is.list(formula) && (is.null(family$nlp) || family$nlp != 
            gp$nlp)) || (!is.list(formula) && !is.null(family$npl) && 
            (family$npl > 1))) 
            stop("incorrect number of linear predictors for family")
        G$terms <- terms
        G$mf <- mf
        G$cl <- cl
        G$am <- am
        if (is.null(G$offset)) 
            G$offset <- rep(0, G$n)
        G$min.edf <- G$nsdf
        if (G$m) 
            for (i in 1:G$m) G$min.edf <- G$min.edf + G$smooth[[i]]$null.space.dim
        G$formula <- formula
        G$pred.formula <- gp$pred.formula
        environment(G$formula) <- environment(formula)
    }
    else {
        if (!is.null(sp) && any(sp >= 0)) {
            if (is.null(G$L)) 
                G$L <- diag(length(G$sp))
            if (length(sp) != ncol(G$L)) 
                stop("length of sp must be number of free smoothing parameters in original model")
            ind <- sp >= 0
            spind <- log(sp[ind])
            spind[!is.finite(spind)] <- -30
            G$lsp0 <- G$lsp0 + drop(G$L[, ind, drop = FALSE] %*% 
                spind)
            G$L <- G$L[, !ind, drop = FALSE]
            G$sp <- rep(-1, ncol(G$L))
        }
    }
    if (!fit) {
        class(G) <- "gam.prefit"
        return(G)
    }
    G$conv.tol <- control$mgcv.tol
    G$max.half <- control$mgcv.half
    object <- estimate.gam(G, method, optimizer, control, in.out, 
        scale, gamma, nei = nei, ...)
    if (!is.null(G$L)) {
        object$full.sp <- as.numeric(exp(G$L %*% log(object$sp) + 
            G$lsp0))
        names(object$full.sp) <- names(G$lsp0)
    }
    names(object$sp) <- names(G$sp)
    object$paraPen <- G$pP
    object$formula <- G$formula
    if (is.list(object$formula)) 
        attr(object$formula, "lpi") <- attr(G$X, "lpi")
    object$var.summary <- G$var.summary
    object$cmX <- G$cmX
    object$model <- G$mf
    object$na.action <- attr(G$mf, "na.action")
    object$control <- control
    object$terms <- G$terms
    object$pred.formula <- G$pred.formula
    attr(object$pred.formula, "full") <- reformulate(all.vars(object$terms))
    object$pterms <- G$pterms
    object$assign <- G$assign
    object$contrasts <- G$contrasts
    object$xlevels <- G$xlevels
    object$offset <- G$offset
    if (!is.null(G$Xcentre)) 
        object$Xcentre <- G$Xcentre
    if (control$keepData) 
        object$data <- data
    object$df.residual <- nrow(G$X) - sum(object$edf)
    object$min.edf <- G$min.edf
    if (G$am && !(method %in% c("REML", "ML", "P-ML", "P-REML"))) 
        object$optimizer <- "magic"
    else object$optimizer <- optimizer
    object$call <- G$cl
    class(object) <- c("gam", "glm", "lm")
    if (is.null(object$deviance)) 
        object$deviance <- sum(residuals(object, "deviance")^2)
    names(object$gcv.ubre) <- method
    environment(object$formula) <- environment(object$pred.formula) <- environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
    if (!is.null(object$model)) 
        environment(attr(object$model, "terms")) <- .GlobalEnv
    if (!is.null(attr(object$pred.formula, "full"))) 
        environment(attr(object$pred.formula, "full")) <- .GlobalEnv
    object
}
