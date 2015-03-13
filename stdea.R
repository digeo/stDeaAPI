#' Spatio-Temporal DEA efficiency
#' 
#' Estimates a ST-DEA frontier and calculates efficiency measures a la Farrell.
#' @param X Inputs of firms to be evaluated, a K x m matrix of observations of K firms with m inputs (firm x input). In case TRANSPOSE=TRUE the input matrix is transposed to input x firm.
#' @param Y Outputs of firms to be evaluated, a K x n matrix of observations of K firms with n outputs (firm x input). In case TRANSPOSE=TRUE the output matrix is transposed to output x firm.
#' @param RTS Text string or a number defining the underlying DEA technology / returns to scale assumption.
#' @param ORIENTATION Input efficiency "in" (1), output efficiency "out" (2), and graph efficiency "graph" (3). For use with DIRECT , an additional option is "in-out" (0).
#' @return The results are returned in a Farrell object with the following components. The last three components in the list are only part of the object when SLACK=TRUE.
#' 
#' @export
stdea <- function(X, Y, RTS = "vrs", ORIENTATION = "out"){
  #proc.time
  ptm0 <- proc.time()
  
  rts <- c("fdh", "vrs", "drs", "crs", "irs", "irs2", "add", "fdh+", "fdh++", "fdh0")
  
  if (missing(RTS)) 
    RTS <- "vrs"
  if (is.numeric(RTS)) {
    RTStemp <- rts[1 + RTS]
    RTS <- RTStemp
  }
  RTS <- tolower(RTS)
  if (!(RTS %in% rts)) 
    stop(paste("Unknown scale of returns:", RTS))
  
  orientation <- c("in", "out")
  
  if (is.numeric(ORIENTATION)) {
    ORIENTATION_ <- orientation[ORIENTATION + 1]
    ORIENTATION <- ORIENTATION_
  }
  
  ORIENTATION <- tolower(ORIENTATION)
  
  if (!(ORIENTATION %in% orientation)) {
    stop(paste("Unknown value for ORIENTATION:", ORIENTATION))
  }
  
  # Define the M
  M <- 100000
  
  # Step
  stp <<- 0.01
  
  repeats <- 1 / stp
  
  numberOfDMUs <<- dim(imported.data)[1]
  numberOfOutputs <<- dim(Y)[2]
  numberOfInputs <<- dim(X)[2]
  
  minimumDMUs <<- max(numberOfInputs * numberOfOutputs, 3 * (numberOfInputs + numberOfOutputs))
  stopDMU <- minimumDMUs
  
  if(minimumDMUs > numberOfDMUs){
    stop("Not enough Decision Making Units, you should have at least : ", minimumDMUs)
  }
  
  alphaMatrix <- matrix(data = 0, nrow = numberOfDMUs, ncol = numberOfDMUs)
  phiDEA <- matrix(data = 0, nrow = numberOfDMUs - stopDMU + 1, ncol = 2)
  
  # Run DEA and calculate lambda(stored at alphaMatrix) for each DMU
  
  while(stopDMU <= numberOfDMUs){
    tempX <- X[1:stopDMU,]
    tempY <- Y[1:stopDMU,]
    e <- dea(tempX, tempY, RTS, ORIENTATION)
    phiDEA[stopDMU - minimumDMUs + 1, 1] <- stopDMU
    phiDEA[stopDMU - minimumDMUs + 1, 2] <- eff(e)[stopDMU]
    l <- lambda(e)
    
    columnNames <- colnames(l)
    dataFrameLambdas <- data.frame(l)
    
    for (i in 1:length(columnNames)){
      alphaMatrix[stopDMU, as.integer(substr(columnNames[i], 2, nchar(columnNames[i])))] <- dataFrameLambdas[nrow(dataFrameLambdas), i]
    }
    
    stopDMU <- stopDMU + 1
  }
  
  alphaMax <- matrix(data = 0, nrow = numberOfDMUs, ncol = 1)
  # alphaMaxCol <- matrix(data = 0, nrow = nrow(alphaMatrix), ncol = 1)
  
  # Find the max lambda for each DMU
  for (i in minimumDMUs:numberOfDMUs){
    alphaMax[i] <- max(alphaMatrix[i,])
  }
  
  deltaMatrix <- matrix(data = M, nrow = numberOfDMUs, ncol = numberOfDMUs)
  
  # Calculating the Delta matrix
  for (i in minimumDMUs:numberOfDMUs){
    for (j in 1:numberOfDMUs){
      if (alphaMatrix[i,j] > 0){
        deltaMatrix[i,j] <- (i - j)
      }
    }
  }
  
  deltaMin <- matrix(data = 0, nrow = numberOfDMUs, ncol = 1)
  
  # Find the DeltaMin
  for (i in minimumDMUs:numberOfDMUs){
    deltaMin[i] <- 0
    
    for (j in 1:i){
      if((deltaMatrix[i, j] != M) && (deltaMatrix[i, j] > deltaMin[i])){
        deltaMin[i] <- deltaMatrix[i, j]
      }
    }
    if(deltaMin[i] == 0){
      deltaMin[i] <- 1
    }
  }
  
  # Define the constraintTypes
  constraintTypes <- matrix(data = ">=", nrow = (numberOfOutputs + 1), ncol = 1)
  
  # Define the rhs
  rightHandSide <- matrix(data = 0, nrow = numberOfOutputs, ncol = 1)
  
  # STDEAOutput <- array(dim = c((numberOfDMUs - minimumDMUs + 1), (1 / stp + 1), (numberOfDMUs + 3)))
  STDEAOutput <- matrix(data = NA, nrow = ((numberOfDMUs - minimumDMUs + 1) * (repeats + 1)), ncol = (8 + 2 * numberOfOutputs))
  
  technical.efficiency <- matrix(data = NA, nrow = ((numberOfDMUs - minimumDMUs + 1) * (repeats + 1)), ncol = numberOfOutputs)
  
  repCounter <- 1
  
  # ST-DEA
  for (i in minimumDMUs:numberOfDMUs){
    objectiveFunction <- matrix(data = 0, nrow = i, ncol = 1)
    
    
    for(Wsp in seq(0, 1, stp)){
      Wt <- 1 - Wsp
      
      lpmodel <- make.lp((numberOfOutputs + 2), (i + 1))
      
      # The 1st column, is the column for the phi
      set.column(lpmodel, 1, c((-Y[i,]), 1, 0))
      
      for (j in 1:i){
        objectiveFunction[j] <- (((Wsp / alphaMax[i]) * alphaMatrix[i, j]) - ((Wt / deltaMin[i]) * deltaMatrix[i, j]))
        set.column(lpmodel, (j + 1), c((Y[j,]), 0, 1))
      }
      
      set.objfn(lpmodel, c(0, objectiveFunction[1:i]))
      set.constr.type(lpmodel, c(constraintTypes[1:(numberOfOutputs + 1)], "="))
      set.rhs(lpmodel, c(rightHandSide[1:numberOfOutputs], 1, 1))
      
      set.type(lpmodel, 2:(i + 1), "binary")
      
      # Set sense for the Linear Problem
      lp.control(lpmodel, sense = "max")
      
      solve(lpmodel)
      #       write.lp(lpmodel, paste("lpfilename", Wsp * 100, ".lp", sep = ""), "lp")
      
      STDEAOutput[repCounter, 1] <- i
      STDEAOutput[repCounter, 2] <- Wsp
      STDEAOutput[repCounter, 3] <- Wt
      STDEAOutput[repCounter, 4] <- (which.max(get.variables(lpmodel)[2:i]))
      
      for(ii in 1:numberOfOutputs){
        technical.efficiency[repCounter, ii] <- Y[(which.max(get.variables(lpmodel)[2:i])), ii] / Y[i, ii]
      }
      
      STDEAOutput[repCounter, 5] <- phiDEA[i - minimumDMUs + 1, 2]
      STDEAOutput[repCounter, 6] <- min(technical.efficiency[repCounter,])
      STDEAOutput[repCounter, 7] <- 1 / STDEAOutput[repCounter, 5]
      STDEAOutput[repCounter, 8] <- 1 / STDEAOutput[repCounter, 6]
      
      for(ii in 1:numberOfOutputs){
        STDEAOutput[repCounter, 8 + ii] <- Y[(which.max(get.variables(lpmodel)[2:i])), ii] - STDEAOutput[repCounter, 6]  * Y[numberOfDMUs, ii]
        STDEAOutput[repCounter, 8 + numberOfOutputs + ii] <- STDEAOutput[repCounter, 6] * Y[i, ii] +  STDEAOutput[repCounter, 6 + ii]
      }
      
      #repCounter <- (i - minimumDMUs) * (repeats + 1) + Wsp * repeats + 1
      #(i - minimumDMUs) * (repeats + 1) + Wsp * repeats + 1
      
      repCounter <- repCounter + 1
      
      rm(lpmodel)
    }
  }
  
  output.column.names <- matrix(nrow = 1, ncol = (8 + 2 * numberOfOutputs))
  
  output.column.names[1] <- "DMU"
  output.column.names[2] <- "Wsp"
  output.column.names[3] <- "Wt"
  output.column.names[4] <- "Peer DMU"
  output.column.names[5] <- "Phi (DEA)"
  output.column.names[6] <- "Phi (S-T DEA)"
  output.column.names[7] <- "Technical Efficiency (DEA)"
  output.column.names[8] <- "Technical Efficiency (S-T DEA)"
  
  for(i in 1:numberOfOutputs){
    output.column.names[8 + i] <- paste("Slack", as.character(i), sep = " ")
    output.column.names[8 + numberOfOutputs + i] <- paste("Projected Value", as.character(i), sep = " ")
  }
  
  colnames(STDEAOutput) <- output.column.names
  
  return(STDEAOutput)
}