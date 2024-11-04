Predict.matrix3<-function (object, data) 
{
    dk <- ExtractData(object, data, NULL)
    X <- Predict.matrix(object, dk$data)
    ind <- attr(dk$data, "index")
    list(X = X, ind = ind)
}
