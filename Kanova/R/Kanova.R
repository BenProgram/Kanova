# Title     : TODO
# Objective : TODO
# Created by: bkaso
# Created on: 23-Sep-20


#' Title
#'
#' @param Ymatrix vector of your response variable
#' @param Xmatrix matrix of your independent variables
#' @param alfa power of your test
#' @param nms name list of your independent variables
#' @param style how should your table be formatted
#' there are four styles: latex output, pander format, kable format, ggpubr  output
#'
#' @return Anova table
#' @export
#'
#' @examples
#' nms <- c("A", "B", "C")
#' Kanova(Y,X,alf,nms,style ="native")
#' Kanova(Y,X,alf,style ="native")
#' Kanova(Y,X,alf,style = "latex")
#' Kanova(Y,X,alf,nms,"pander")
#' Kanova(Y,X,alf,style = "kable")
#' Kanova(Y,X,alf,nms,"table")

library(ggpubr)
library(xtable)
library(knitr)
library(pander)

Kanova <- function(Ymatrix, Xmatrix, alfa = 0.05, nms = NULL, style = NULL) {

    n <- nrow(Ymatrix)
  p <- ncol(Xmatrix)

  #Sum of Squares

  SST <- SST(Ymatrix)
  SSE <- SSE(Ymatrix, Xmatrix)
  SSR <- SSR(Ymatrix, Xmatrix)

  #Mean Squares, F and R-squared

  MSE <- SSE / (n - p)
  MSR <- SSR / (p - 1)
  MST <- SST / (n - 1)

  #F-value and R-squares
  F_Value <- MSR / MSE
  R_Squared <- SSR / SST
  Adj_R_Squared <- 1 - MSE / MST
  pval <- PV(F_Value, n, p)

  Regression <- c(p - 1, round(SSR, 3), round(MSR, 3), round(F_Value, 3))
  Error <- c(n - p, round(SSE, 3), round(MSE, 3), "")
  Total <- c(n - 1, round(SST, 3), "", "")

  Correlation <- round(R_Squared, 4)
  Adjusted <- round(Adj_R_Squared, 4)

  if (alfa == 0) {
    alfa = 0.05
    F_critical <- round(qf(1 - alfa, p - 1, n - p), 4)
  }else {
    F_critical <- round(qf(1 - alfa, p - 1, n - p), 4)
  }

  anv <- (rbind(Regression, Error, Total))
  expl <- cbind(Correlation, Adjusted, F_critical, round(pval, 4))
  anv <- rbind(anv, c("R-sqr", "R-Adj", "F-Crit", "Pr(>F"), expl)
  anv <- as.table(anv)
  colnames(anv) <- (c("df", "Sum sQ", "Mean Sq", "F Value"))
  names(dimnames(anv)) <- list("", "\nAnalysis of Variance Table\n")

  if (style == "native")
  {
    print(anv)
    print(Coefs(Ymatrix, Xmatrix, nms))
  }
  else if (style == "table")
  {
    if (F_Value >= qf(1 - alfa, p - 1, n - p)) {
      main.title <- "Analysis of Variance Table"
      tab <- ggtexttable((anv), theme = ttheme("light", base_size = 13))
      tab %>%
        tab_add_title(text = main.title, face = "bold", padding = unit(5, "line")) %>%
        tab_add_footnote(text = paste("*We reject H0 at the ", alf, " level of significance"), size = 13, face = "italic")
    }
    else
    {
      main.title <- "Analysis of Variance Table"
      tab <- ggtexttable((anv), theme = ttheme("light", base_size = 13))
      tab %>%
        tab_add_title(text = main.title, face = "bold", padding = unit(5, "line")) %>%
        tab_add_footnote(text = paste("*We fail to reject H0 at the ", alf, " level of significance"), size = 13, face = "italic")
    }
  }
  else if (style == "kable") {

    print(kable(anv))
    kable(Coefs(Ymatrix, Xmatrix, nms))
  }
  else if (style == "pander") {

    print(pander(anv))
    pander(Coefs(Ymatrix, Xmatrix, nms))
  }
  else if (style == "latex")
  {

    print(xtable(anv))
    xtable(Coefs(Ymatrix, Xmatrix, nms))

  }
  else {

    print(anv)

  }
}


SST <- function(Y) {
  sst <- t(Y) %*% Y - (1 / n) %*% t(Y) %*% J %*% Y
  return(sst)
}


SSR <- function(Y, X) {

  ssr <- t(Y) %*%
    X %*%
    solve(t(X) %*% X) %*%
    t(X) %*%
    Y - (1 / n) %*% t(Y) %*% J %*% Y
  return(ssr)
}

SSE <- function(Y, X) {
  SSE <- t(Y) %*% Y - t(Y) %*%
    X %*%
    solve(t(X) %*% X) %*%
    t(X) %*%
    Y
}

PV <- function(Fval, n, p) {
  pv <- 1 - pf(Fval, p - 1, n - p)
  return(pv)
}

Coefs <- function(Y, X, nms) {
  if (is.null(nms)) {
    fit <- lm(Y ~ X[, -1])
    smr <- summary(fit)$coef
    smr <- as.table(smr)
    rownames(smr) <- c("Intercept", LETTERS[1:ncol(X) - 1])
    names(dimnames(smr)) <- list("", "\nCoefficients Test of Significance\n")
  }else {
    fit <- lm(Y ~ X[, -1])
    smr <- summary(fit)$coef
    smr <- as.table(smr)
    names(dimnames(smr)) <- list("", "\n Coefficients Test of Significance \n") #make sure your response is last
    rownames(smr) <- c("Intercept", nms)
  }
  return(smr)
}
