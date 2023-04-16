\name{predict.lm.from.rf}
\alias{predict.lm.from.rf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for rflasso
}
\description{
Make predictions, generate RMSE from pooled terminal leaves, and obtain lasso coefficients. Supports "lasso" where leaves from the Random Forest model are pooled and a lasso model is fit to the data. Also, supports "Average" where the standard Random Forest prediction is returned for benchmarking purposes. Lasso is tuned by GCV using the "msgps" package.
}
\usage{
predict.lm.from.rf(bagged.model,newdata,model="lasso",predictions.only=FALSE,number.of.cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bagged.model}{Compatible with models generated by the "build.tree" function.}
  \item{newdata}{New data as an R data.frame for which predictions are to be made.}
  \item{model}{Either "lasso" or "Average".}
  \item{predictions.only}{If TRUE, then the function will only return predictions for each sample in "newdata". When FALSE, this function will also return lasso coefficients (if applicable) and RMSE of in sample data.}
  \item{number.of.cores}{Number of cores to be utilized.}
}


\author{
Alexander Wendelborn
apwendelborn@gmail.com
}
\note{
Generating lasso models from the pooled terminal leaves sometimes outperforms the mean predictions generated by the standard Random Forest. Use the function "out.of.bag.error" to see which performs better on out of bag samples.
}


\seealso{
build.tree, out.of.bag.error, auto.optimize
}