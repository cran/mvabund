## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(knitr)

## ----include=FALSE------------------------------------------------------------
worms <- data.frame(matrix(nrow = 3, ncol = 10) )
rownames(worms) <- c("Treatment", "Count", "# pitfalls")
names(worms) <- NULL

worms[1,] <- c("C", "R", "R", "R", "C", "R", "R", "R", "R", "R")
worms[2,] <- c(0,3,1,3,1,2,12,1,18,0)
worms[3,] <- c(5,5,5,5,5,5,5,4,5,5)

## ----worm counts--------------------------------------------------------------
knitr::kable(worms, "simple")

## ----setup--------------------------------------------------------------------
library(mvabund)
library(ecostats)

data(reveg) 
attach(reveg)

skimr::skim(abund$Haplotaxida) # Great function to get an overview of the data

## ----fit model----------------------------------------------------------------
worms_offset <-  manyglm(Haplotaxida~treatment+offset(log(pitfalls)), family="negative.binomial", data=abund)

anova(worms_offset)

