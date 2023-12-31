% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob_metrics.R
\name{prob_metrics_sf}
\alias{prob_metrics_sf}
\alias{average_precision.sf}
\alias{brier_class.sf}
\alias{classification_cost.sf}
\alias{gain_capture.sf}
\alias{mn_log_loss.sf}
\alias{pr_auc.sf}
\alias{roc_auc.sf}
\alias{roc_aunp.sf}
\alias{roc_aunu.sf}
\title{Probability metrics for \code{sf} objects}
\usage{
\method{average_precision}{sf}(data, ...)

\method{brier_class}{sf}(data, ...)

\method{classification_cost}{sf}(data, ...)

\method{gain_capture}{sf}(data, ...)

\method{mn_log_loss}{sf}(data, ...)

\method{pr_auc}{sf}(data, ...)

\method{roc_auc}{sf}(data, ...)

\method{roc_aunp}{sf}(data, ...)

\method{roc_aunu}{sf}(data, ...)
}
\arguments{
\item{data}{an \link[sf:sf]{sf::sf} object}

\item{...}{any other parameters to pass to the \code{data.frame} version of
the metric.}
}
\value{
A tibble with columns \code{.metric}, \code{.estimator}, and \code{.estimate}
and 1 row of values.
}
\description{
\code{tidysdm} provides methods to handle \link[sf:sf]{sf::sf} objects for the following
\link{yardstick} metrics:

\code{\link[yardstick:average_precision]{yardstick::average_precision()}}

\code{\link[yardstick:brier_class]{yardstick::brier_class()}}

\code{\link[yardstick:classification_cost]{yardstick::classification_cost()}}

\code{\link[yardstick:gain_capture]{yardstick::gain_capture()}}

\code{\link[yardstick:mn_log_loss]{yardstick::mn_log_loss()}}

\code{\link[yardstick:pr_auc]{yardstick::pr_auc()}}

\code{\link[yardstick:roc_auc]{yardstick::roc_auc()}}

\code{\link[yardstick:roc_aunp]{yardstick::roc_aunp()}}

\code{\link[yardstick:roc_aunu]{yardstick::roc_aunu()}}
}
