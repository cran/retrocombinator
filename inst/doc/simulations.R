## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(retrocombinator)

## -----------------------------------------------------------------------------
simulateEvolution()

## ---- eval = FALSE------------------------------------------------------------
#  # --- NOT RUN (this is an alternative) ---
#  recombParams <- RecombParams(recombMean = 1.0)
#  activityParams <- ActivityParams(lengthCriticalRegion = 20,
#                                   probInactiveWhenMutated = 0.1)
#  simulateEvolution(recombParams = recombParams, activityParams = activityParams)
#  # ----------------------------------------

## ---- eval = FALSE------------------------------------------------------------
#  # --- NOT RUN (this is an alternative) ---
#  
#  # Obtain your sequence in your favourite way
#  library(Biostrings)
#  fastaInput <- readDNAStringSet('path/to/your/FASTA/file')
#  yourSequence <- toString(fastaInput$yourSequence)
#  
#  # Override the sequence parameters and pass it to the simulation
#  sequenceParams <- SequenceParams(initialSequence = yourSequence)
#  simulateEvolution(sequenceParams = sequenceParams)
#  # ----------------------------------------

## -----------------------------------------------------------------------------
data <- parseSimulationOutput('simulationOutput.out')


## ---- eval = FALSE------------------------------------------------------------
#  # --- NOT RUN (this is an alternative) ---
#  
#  # Alternatively, using the pipe operator
#  library(magrittr) # For %>%
#  data <- simulateEvolution() %>%
#    parseSimulationOutput()
#  
#  # ----------------------------------------

## -----------------------------------------------------------------------------
print(names(data$params))
# Prints out all the parameters the simulation was run with

print(colnames(data$sequences))
# step - timestep in the simulation
# realTime - time since the start of the simulation (in millions of years)
# sequenceId - the unique ID of the sequence (to track it over time); initial
#              sequence has sequenceId 0 to (numInitialCopies-1)
# parentMain - the unique ID of the sequence this burst from; 
#              (-1 if nothing)
# parentOther - the unique ID of the sequence its parent recombined with; 
#               (-1 if nothing)
# distanceToInitial - the distance to the initial sequence
# isActive - whether or not the sequence is capable of transposition

print(colnames(data$pairwise))
# step - timestep in the simulation
# realTime - time since the start of the simulation (in millions of years)
# sequenceId1 - an ID of a sequence present at the time
# sequenceId2 - an ID of a different sequence present at the time; not all pairs
#               are given - that is, for sequences a and b, either (a, b) or (b, a) 
#               is present as a row but not both
# distancePairwise - the distance between the two sequences

print(colnames(data$familyRepresentatives))
# step - timestep in the simulation
# realTime - time since the start of the simulation (in millions of years)
# familyId - the unique ID of the family representative (to track it over time)
# creationTime - time this family representative was created (in millions of years)
# sequenceId - a sequence belonging to that family (sequences belonging to a
#              family are listed as separate rows, and all sequences belonging
#              to that family are listed)

print(colnames(data$familyPairwise))
# step - timestep in the simulation
# realTime - time since the start of the simulation (in millions of years)
# familyId1 - an ID of a family representative present at the time
# familyId2 - an ID of a different family representative present at the time;
#             not all pairs are given - that is, for sequences a and b, either (a, b) or (b,
#             a) is present as a row but not both
# distancePairwise - the distance between the two family representatives

## ---- eval = FALSE------------------------------------------------------------
#  summariseEvolution(data) # TODO

## -----------------------------------------------------------------------------
plotEvolution(data, "initial")         # Distance to initial sequence
plotEvolution(data, "pairwise")        # Pairwise distances between sequences
plotEvolution(data, "families")        # Family sizes

## -----------------------------------------------------------------------------
sequenceParams <- SequenceParams(initialSequence = "TCAGTCAGTCAGTCAGTGTG",
                                 numInitialCopies = 10)
activityParams <- ActivityParams(lengthCriticalRegion = 2,
                                 probInactiveWhenMutated = 0.1)
mutationParams <- MutationParams(model = "F81")
burstParams <- BurstParams(burstProbability = 0.2,
                           burstMean = 2,
                           maxTotalCopies = 40)
recombParams <- RecombParams(recombMean = 2.0,
                             recombSimilarity = 0.85)
selectionParams <- SelectionParams(selectionThreshold = 0.25)
familyParams <- FamilyParams(familyCoherence = 0.65,
                             maxFamilyRepresentatives = 15)
simulationParams <- SimulationParams(numSteps = 35,
                                     timePerStep = 1.5)
outputParams <- OutputParams(outputFilename = 'simulationExampleParameters.out',
                             outputNumInitialDistance = 5,
                             outputNumPairwiseDistance = 5,
                             outputNumFamilyLabels = 5,
                             outputNumFamilyMatrix = 5,
                             outputMinSimilarity = 0.1)
seedParams <- SeedParams(toSeed = TRUE, seedForRNG = 1)
simulateEvolution(sequenceParams = sequenceParams,
                  activityParams = activityParams,
                  mutationParams = mutationParams,
                  burstParams = burstParams,
                  recombParams = recombParams,
                  selectionParams = selectionParams,
                  familyParams = familyParams,
                  simulationParams = simulationParams,
                  outputParams = outputParams)

