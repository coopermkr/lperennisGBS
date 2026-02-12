#'''''''''''''''''''''''''''''''''''''
#' Vortex Parameters
#' @author Cooper Kimball-Rhines
#' @date 2026-01-13
#'''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)
library(vortexR)

# Load in translocation vortex runs
ALsup <- collate_dat(project = "translocations_AlbanySupplement", runs = 500,
                        dir_in = "5.translocation/VOutput/translocations/VOutput/",
                        dir_out = "6.vortexAnalysis/")

SAsup <- collate_dat(project = "translocations_SaratogaSupplement", runs = 250,
                     dir_in = "5.translocation/VOutput/translocations/VOutput/",
                     dir_out = "6.vortexAnalysis/")

MOsup <- collate_dat(project = "translocations_AlbanySupplement", runs = 250,
                     dir_in = "5.translocation/VOutput/translocations/VOutput/",
                     dir_out = "6.vortexAnalysis/")

effAL <- Ne(data = ALsup, gen = 4, yr0 = 1, yrt = 160, dir_out = "6.vortexAnalysis/")
