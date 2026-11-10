library(tidyverse)
library(rlang)
library(zoo)

source("functions/auc_poly.R")
source("functions/drug_concentration_app.R")

test_that("24-hour AUC/MIC uses trapezoidal integration", {
  triangle <- data.frame(time = c(0, 1, 2), concentration = c(0, 2, 0))
  result <- auc_poly(triangle, MIC = 2, dosing_interval = 2, n_interval = 1)

  # Triangle area = 2 mg*h/L; / MIC 2 mg/L; scaled from 2 h to 24 h.
  expect_equal(result$AUC, 12)
  expect_length(result$x_poly, length(result$y_poly))
})

test_that("AUC is independent of the supplied sampling resolution", {
  coarse <- data.frame(time = seq(0, 8, by = 1))
  coarse$concentration <- coarse$time
  fine <- data.frame(time = seq(0, 8, by = 0.25))
  fine$concentration <- fine$time

  expect_equal(
    auc_poly(coarse, MIC = 1, dosing_interval = 8, n_interval = 1)$AUC,
    auc_poly(fine, MIC = 1, dosing_interval = 8, n_interval = 1)$AUC
  )
})

test_that("time above MIC handles clear limiting cases", {
  below <- drug_concentration_app(dose = 0, MIC = 1)
  above <- drug_concentration_app(dose = 1000, MIC = 1e-9)

  expect_equal(below$time_above_mic, 0)
  expect_equal(above$time_above_mic, 1, tolerance = 0.01)
})

test_that("AUC remains available when concentration never exceeds MIC", {
  result <- drug_concentration_app(dose = 0, MIC = 1, pk_pd_index = "24-AUC")

  expect_equal(result$auc_poly$AUC, 0)
})

test_that("invalid clinical and numerical inputs fail early", {
  expect_error(drug_concentration_app(t.5 = 0), "greater than zero")
  expect_error(drug_concentration_app(protein_binding = 40), "fraction")
  expect_error(drug_concentration_app(MIC = 0), "greater than zero")
  expect_error(drug_concentration_app(n_intervals = 1.5), "whole number")
  expect_error(
    drug_concentration_app(infusion_time = 9, dosing_interval = 8),
    "cannot exceed"
  )
})
