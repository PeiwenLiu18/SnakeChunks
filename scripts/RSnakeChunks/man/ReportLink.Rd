% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_link.R
\name{ReportLink}
\alias{ReportLink}
\title{IncludeLink}
\usage{
ReportLink(source, target, displayPath = TRUE, report.text = NULL)
}
\arguments{
\item{source}{path of the report file, where the link will be inserted}

\item{target}{path of the file to link}

\item{displayPath=TRUE}{in the link text, display either the path or the basename of the target}

\item{report.text=NULL}{if not null, append the link to it}
}
\description{
IncludeLink
}
\examples{

## Display the (relative) path  in the link
ReportLink(
  source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
  target = "RNA-seq/results/diffexpr/a_vs_b/DEG/deg_table.tsv",
  displayPath = TRUE)

## Display only the file name, not its path
ReportLink(
  source = "RNA-seq/results/diffexpr/reports/myreport.Rmd", 
  target = "RNA-seq/results/diffexpr/a_vs_b/DEG/deg_table.tsv",
  displayPath = FALSE)

}
\author{
Jacques van Helden
}
