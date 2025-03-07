% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidybulk_functions.R
\name{make_ciber_volcano_table}
\alias{make_ciber_volcano_table}
\title{return dataframe of cell proportions log2FC of cibersort plot of cell type proportions.}
\usage{
make_ciber_volcano_table(
  ciber_res_paired_df,
  factor_one,
  factor_two,
  wilcox_df
)
}
\arguments{
\item{ciber_res_paired_df}{is dataframe with cell proportions of two conditions}

\item{factor_one}{is the name of the column in ciber_res_paired representing condition one}

\item{factor_two}{is the name of the column in ciber_res_paired representing condition one. Will calculate foldchange of factor_two/factor_one}

\item{wilcox_df}{is dataframe of listing statistic, p.value, method, alternative (hypothesis) of cell_type proportions from cibersort results}
}
\description{
return dataframe of cell proportions log2FC of cibersort plot of cell type proportions.
}
