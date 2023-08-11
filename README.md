# rflasso
Fits random forests and then makes predictions with lasso based on pooled samples from the terminal nodes.

To download this package, first download the devtools library and load it. Next, run install_github("apwendelborn/rflasso").

The general purpose of this package is to provide a model which improves upon random forest. In a standard random forest, predictions are made by finding the terminal nodes associated with a new data point. The dependent variable for all of these terminal nodes is then averaged to make a prediction. rflasso performs a similar operation except instead of averaging the dependent variable and discarding all the information on the independent variables, a lasso model is fit to the pooled data and predictions are made from this linear model. This provides a significant improvement in performance over random forest. rflasso rivals "boosting" in terms of performance while providing much greater model interpretability.

This implementation works best with datasets of less than 10,000 observations.
