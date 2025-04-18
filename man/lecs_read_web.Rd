% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/web.R
\name{lecs_read_web}
\alias{lecs_read_web}
\title{Read LECS data from website}
\usage{
lecs_read_web(start_date = NULL, base_url = NULL)
}
\arguments{
\item{start_date}{Date to get date from. Date object or string coerceable}

\item{base_url}{URL to query. Defaults to WHOI LECS data
to a date}
}
\value{
a single column data frame
}
\description{
Reads HTML table from LECS_data into a single column data frame
}
