`smooth_dim.gam.prefit`<-function(obj){
  vapply(obj[["smooth"]],FUN=`[[`,FUN.VALUE=integer(1),"dim");
}

`smooth_terms.gam.prefit`<-function (object, ...) {
    lapply(object[["smooth"]], `[[`, "term")
}

`parametric_terms.gam.prefit`<-function (model, ...) {
    lss_terms <- function(i, terms) {
        labs <- labels(delete.response(terms[[i]]))
        names(labs) <- unlist(lapply(labs, function(lab, i) {
            lab <- if (i > 1L) {
                paste0(lab, ".", i - 1)
            }
            else {
                lab
            }
            lab
        }, i = i))
        labs
    }
    tt <- model$pterms
    if (is.list(tt)) {
        labs <- unlist(lapply(seq_along(tt), lss_terms, terms = tt))
    }
    else {
        if (length(attr(tt, "term.labels") > 0L)) {
            tt <- delete.response(tt)
            labs <- labels(tt)
            names(labs) <- labs
        }
        else {
            labs <- character(0)
        }
    }
    labs
}

`family.gam.prefit`<-function(obj){
  return(obj$family);
}



