## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align='center'
  )

## ----setup--------------------------------------------------------------------
library(mvabund)

data(Tasmania) 
attach(Tasmania)

skimr::skim(copepods) # Great function to get an overview of the data

## -----------------------------------------------------------------------------
copepod_abund <- mvabund(copepods)

## ----fig.height=6, fig.width= 7-----------------------------------------------
plot(copepod_abund~treatment, col = block)

## -----------------------------------------------------------------------------
cope.nb <- manyglm(copepods ~ treatment*block, family =  "negative.binomial")

## ----fig.height=5, fig.width= 6-----------------------------------------------
plot(cope.nb)

## ----fig.height=5, fig.width= 6-----------------------------------------------
meanvar.plot(copepods~tr.block, col = treatment)

## ----c------------------------------------------------------------------------
anova(cope.nb, p.uni = "adjusted")

## -----------------------------------------------------------------------------
summary(cope.nb) 

## -----------------------------------------------------------------------------
predict(cope.nb, type = "response") 

