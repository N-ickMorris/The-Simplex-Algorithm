# --------------------------------------------------------------------------------------------------------------------
# dual simplex algorithm for linear programs of canonical form:
    # min or max: Cx
    # s.t. Ax = b
    # s.t. x >= 0 
# --------------------------------------------------------------------------------------------------------------------

# REFER TO: http://www.egwald.ca/operationsresearch/lpdualsimplex.php

# --- Inputs for the function: dsimplx ----

# type is a string indicating the type of optimization to perform: "min" or "max"
# C is a row vector of objective function coefficients
# A is a matrix of contraint coeffients
# b is a column vector of right hand side scalars
# slack is a row vector of 1's or 0's identifying which variables are slack variables 
    # s.t. 1 indicates an entry in row vector: x, is currently a slack variable
    # s.t. 0 indicates an entry in row vector: x, is currently a non-slack variable
# artificial is a row vector of 1's or 0's identifying which variables are artificial variables
    # s.t. 1 indicates an entry in row vector: x, is currently an artificial variable
    # s.t. 0 indicates an entry in row vector: x, is currently a non-artificial variable

dsimplx = function(type = "min", C, A, b, slack, artificial)
{
    # make sure the user enters the optimization direction correctly
    
    if(type != "min" && type != "max")
    {
        stop("You must enter 'min' or 'max' for input 'type'")
    }
    
    # ---- setup ----
    
    # mat is a matrix representing the initial tableau
    
    mat = cbind(as.vector(c(0,b)), rbind(as.vector(C),A))
    
    # this program only solves maximization problems so multiply C (ie. mat[1,]) by -1 if type = "max"
        # C gets multiplied by -1 before the program starts, regardless of what type is initially
        # to switch from "min" to "max" requires C to be multiplied by -1
        # therefore it is reduntant to multiply C by -1 if type = "min" initially
    
    if(type == "max")
    {
        mat[1,] = mat[1,] * -1
    }
    
    # basis is a row vector of 1's and 0's representing which variables are and aren't in the initial basis
    
    basis = slack + artificial
    
    # it is a single value that will act as a counter representing the current iteration of the program 
	# tab is a list of matrices where each matrix is a tableau corresponding to an iteration of the program
	# store tableau data
    
    it = 0
    tab = list(mat)
    
    # ---- Phase 0: Drive all artificial variables out of the basis ----
    
    # if any artificial variables are in the model then the following code will move each one out of the basis 
    
    if(any(artificial == 1))
    {
        # x is a row vector indicating which columns of the matrix 'mat' correspond to artificial variables
        
        x = which(artificial == 1) + 1
        
		# the following for-loop evaluates each artificial variable identified in the row vector 'x'
		
        for(k in 1:length(x))
        {
            # i is a single value corresponding to the row position of the current non-basic variable that will become basic
            # j is a single value corresponding to the column position of the current artificial variable that will become non-basic
            # new is a row vector representing the row that will be in the next tableau, this row corresponds to the current non-basic variable taht will become basic
            
            i = which(as.matrix(mat[,x])[,k] == 1)
            j = which(mat[i,-c(1, which(basis == 1))] != 0)[1] + 1
            new = mat[i,] / mat[i,j]
            
            # the following for-loop performs scalar multiplication and vector addition on each row of the current tableau 'mat' 
                        
            for(r in 1:NROW(mat))
            {
                mat[r,] = mat[r,] - (mat[r,j] * new)
            }
            mat[i,] = new
            
			# after each row in matrix 'mat' has been updated, there has been a change of basis and therefore the row vector 'basis' must be updated accordingly
			
			basis[(x[k] - 1)] = 0
            basis[(j - 1)] = 1
            
            # store tableau data
            
            it = it + 1
            tab[(it + 1)] = list(mat)
        }
    }
    
    # ---- Phase 1: Obtain non-negative entries in the first row of the tableu for the primal and slack variables ----
    
    # the following while-loop iterates until all entries for slack and primal variables in the first row of the current tableau, are non-negative
    
    complete = FALSE
    while(complete == FALSE)
    {
        if(mat[1,which.min(mat[1,which(artificial == 0) + 1])[1] + 1] < 0)
        {
            # j is a row vector corresponding to which slack/primal columns in row 1 of the current tableau 'mat' have a negative coefficients
            
            j = which(mat[1,which(artificial == 0) + 1] < 0) + 1
            
			# the following for-loop iterates until a positive entry has been found in any of the column vectors identified by j
			# i is a single value corresponding the row position of the desired positive value
			# j becomes a single value corresponding to the column position of the desired positive value
			
            for(k in 1:length(j))
            {
                i = which(mat[,j[k]] > 0)[1]
                if(length(i) > 0)
                {
                    j = j[k]
                    break
                }
            }
            
			# if no value has been assigned to i then the desired positive value cannot be found because the model is dual infeasible
			
            if(length(i) == 0)
			{
			    stop("Your model is dual infeasible")
			}
			            
            i = which(mat[,j] > 0)[1]
            
			# ratio is a vector corresponding to the results of the minimum ratio test
			# if a ratio is non-positive the it is discarded
			# if a ratio corresponds to an artificial varibale then it is discarded
			# the first entry of the ratio is discarded becuase this column corresponds to the value of the objective function and basic variables
            
            ratio = -mat[1,] / mat[i,]
            ratio = ifelse(ratio <= 0, NA, ratio)
            ratio = ifelse(artificial == 1, NA, ratio)
            ratio[1] = NA
            
            # determine which ratio is the minimum
			# j becomes a single value corresponding to the column position of the next variable to enter the basis
			# new is a row vector representing the row that will be in the next tableau, this row corresponds to the current non-basic variable that will become basic
            
            j = which.min(ratio)
            new = mat[i,] / mat[i,j]
            
            # the following for-loop iterates through the row identified by 'i' until it finds the column corresponding to the current variable in the basis that will leave the basis
			# x is a single value corresponding to the column position of the next variable to leave the basis
			
            x = 0
            for(k in 1:length(which(mat[i,] == 1)))
            {
                if(length(which(mat[-1,which(mat[i,] == 1)[k]] == 0)) == (NROW(mat) - 2))
                {
                    x = which(mat[i,] == 1)[k]
                }
                
                if(x > 0) break
            }
            
            # the following for-loop performs scalar multiplication and vector addition on each row of the current tableau 'mat'
            
            for(r in 1:NROW(mat))
            {
                mat[r,] = mat[r,] - (mat[r,j] * new)
            }
            mat[i,] = new
			
			# after each row in matrix 'mat' has been updated, there has been a change of basis and therefore the row vector 'basis' must be updated accordingly
			
            basis[(x - 1)] = 0
            basis[(j - 1)] = 1
            
            # store tableau data
            
            it = it + 1
            tab[(it + 1)] = list(mat)
            
        } else
        {
            complete = TRUE
        }
        
        if(complete == TRUE) break
    }
    
    # ---- Phase 2: Obtain non-negative entries in the first column of the tableu for the righthand side vector: b ----
    
    # the following while-loop iterates until all entries in the first column of matrix 'mat', excluding the first entry, are non-negative
    
    complete = FALSE
    while(complete == FALSE)
    {
        if(mat[which.min(mat[-1,1]) + 1,1] < 0)
        {
		    # i is a single value corresponding the row position of a negative entry in column 1 that will become non-negative
		
            i = which.min(mat[-1,1]) + 1
			
			# ratio is a vector corresponding to the results of the minimum ratio test
			# if a ratio is non-positive the it is discarded
			# if a ratio corresponds to an artificial varibale then it is discarded
			# the first entry of the ratio is discarded becuase this column corresponds to the value of the objective function and basic variables            
            
			ratio = -mat[1,] / mat[i,]
            ratio = ifelse(ratio <= 0, NA, ratio)
            ratio = ifelse(artificial == 1, NA, ratio)
            ratio[1] = NA
            
			# determine which ratio is the minimum
			# j is a single value corresponding to the column position of the next variable to enter the basis
			# new is a row vector representing the row that will be in the next tableau, this row corresponds to the current non-basic variable that will become basic
			
            j = which.min(ratio)
            new = mat[i,] / mat[i,j]
            
			# the following for-loop iterates through the row identified by 'i' until it finds the column corresponding to the current variable in the basis that will leave the basis
			# x is a single value corresponding to the column position of the next variable to leave the basis
			
            x = 0
            for(k in 1:length(which(mat[i,] == 1)))
            {
                if(length(which(mat[-1,which(mat[i,] == 1)[k]] == 0)) == (NROW(mat) - 2))
                {
                    x = which(mat[i,] == 1)[k]
                }
                
                if(x > 0) break
            }
            
            # the following for-loop performs scalar multiplication and vector addition on each row of the current tableau 'mat'
            
            for(r in 1:NROW(mat))
            {
                mat[r,] = mat[r,] - (mat[r,j] * new)
            }
            mat[i,] = new
            
			# after each row in matrix 'mat' has been updated, there has been a change of basis and therefore the row vector 'basis' must be updated accordingly
			
			basis[(x - 1)] = 0
            basis[(j - 1)] = 1
            
            # store tableau data
            
            it = it + 1
            tab[(it + 1)] = list(mat)
            
        } else
        {
            complete = TRUE
        }
        
        if(complete == TRUE) break
    }
    
    # ---- results ----
    
    # solution is a column vector containing the value of the objective function and the values of the basic variables at optimality
    
    solution = vector(length = NCOL(mat)) * 0
    solution[1] = mat[1,1]
    
    for(i in 2:NCOL(mat))
    {
        solution[i] = sum(c(0,basis)[i] * mat[,1] * mat[,i])
    }
    
    # vars is a column vector containing the identity of each value in the column vector 'solution'
    
    vars = vector(length = NCOL(mat))
    vars[1] = "obj"
    
    for(i in 2:length(vars))
    {
        vars[i] = paste0("x",(i - 1))
    }
    
    solution = data.frame(t(t(solution)))
    rownames(solution) = t(vars)
    colnames(solution) = "value"
    
    # the identity of each tableau in the list 'tab' is labeled according to it's iteration number
        
    names(tab) = 0:it
    
    # the results that are the output of this function, is a list of the total iterations required, the optimal solution, and the dual simplex tableaus
    
    results = list("iterations" = as.numeric(it), "solution" = solution, "tableaus" = tab)
    return(results)
}

