\name{Perform PPR and LDA on a dataframe}
\alias{PPR_LDA}
\title{
  PPR_LDA
}
\description{
  Performs projection pursuit regression (using ppr from stats) and linear discriminant analysis (using lda from MASS) on a dataframe after including or excluding certain columns (PPR requires all columns to be numerical). The PP1, PP2, LD1, and LD2 directions are then dotted with each row to produce PP1, PP2, LD1, and LD2 columns in the new dataframe. The function returns a list of the new dataframe, the PPR and LDA function outputs, and the PPR and LDA direction matrices/lists.
}
