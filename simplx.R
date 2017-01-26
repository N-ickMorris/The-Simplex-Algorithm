# --------------------------------------------------------------------------------------------------------------------
# primal simplex algorithm for linear programs of standard form:
    # min or max: Cx
    # s.t. Ax = b
    # s.t. x >= 0 
# --------------------------------------------------------------------------------------------------------------------

# --- Inputs for the function: simplx ----

# type is a string indicating the type of optimization to perform: "min" or "max"
# C is a row vector of objective function coefficients
# A is a matrix of contraint coeffients
# b is a column vector of right hand side scalars
# x is a row vector of values representing a feasible solution
# basis is a row vector of 1's or 0's identifying which variables are currently basic or non-basic
    # s.t. 1 indicates an entry in row vector: x, is currently a basic variable
    # s.t. 0 indicates an entry in row vector: x, is currently a non-basic variable

simplx = function(type = "min", C, A, b, x, basis)
{
    # make sure the user enters the optimization direction correctly
    
    if(type != "min" && type != "max")
    {
        stop("You must enter 'min' or 'max' for input 'type'")
    }
    
    # ---- setup ----
    
    optimal = FALSE
    it = 0
    
    # n is a scalar of the total number of variables
    # tab is going to become a table of data regarding each iteration of the simplex algorithm:
        # iteration number
        # objective value
        # step size
        # direction of improvement
        # variable values
    
    n = as.numeric(length(basis))
    tab = cbind(it, as.numeric(C%*%x), 0, t(rep(0,n)), t(x))
    
    # ---- iterate through the simplex algorithm until it reaches optimality or unboundedness ----
    
    while(optimal == FALSE)
    {
        # Cb is a row vector of the objective function coefficients of the current basic variables
        # Cn is a row vector of the objective function coefficients of the current non-basic variables
        
        Cb = C[basis == 1]
        Cn = C[basis == 0]
        
        # B is a matrix of the constraint coefficients of the current basic variables
        # N is a matrix of the constraint coefficients of the current non-basic variables
        
        B = A[,basis == 1]
        N = A[,basis == 0]
        
        # xb is a row vector of the values of the current basic variables
        # xn is a row vector of the values of the current non-basic variables
        
        xb = x[basis == 1]
        xn = x[basis == 0]
        
        # Crc is a row vector of the coefficients of reduced cost of the current non-basic variables
        # if the objective function is to be maximized then:
            # if all entries in row vector: Crc, are non-positive then optimality has occured
        # if the objective function is to be minimized then:
            # if all entries in row vector: Crc, are non-negative then optimality has occured
        
        Crc = Cn - Cb%*%solve(B)%*%N
        
        if(type == "max")
        {
            if(all(Crc <= 0))
            {
                optimal = TRUE
            }	
        } else
        {
            if(all(Crc >= 0))
            {
                optimal = TRUE
            }
        }		
        
        # if optimality has occured, then stop iterating though the simplex algorithm
        
        if(optimal == TRUE) break
        
        # db is a column vector of the direction of improvement of the current basic variables
                
        if(type == "max")
        {
            db = -solve(B)%*%N[,Crc == max(Crc)]	
        } else
        {
            db = -solve(B)%*%N[,Crc == min(Crc)]
        }
        
		# if there was a tie between multiple non-basic variable for the most favorable reduced cost coeffiecient, then just use the first non-basic variable
		
		if(NCOL(db) > 1)
		{
		    db = matrix(db[,1])
		}
		
		# if all entries in column vector: db, are non-negative then unboundedness has occured
		
        if(all(db >= 0))
        {
            stop("Your model is Unbounded")
        }	
        
        # ratio is a column vector of the ratios of the current basic variables for the minimum ratio test
        # lam is the chosen step size in the direction of improvement
        # lam is a result of the minimum positive value in the column vector: ratio
        
        ratio = -xb/ifelse(db >= 0, NA, db)
        lam = min(ratio, na.rm = TRUE)
        
        # the following logic computes the change of basis, while finalizing the direction of improvement
        # the non-basic variable identified in the row vector: Crc, will become a basic variable
        # the basic variable identified in the column vector: ratio, will become a non-basic variable
        
        # bv and nbv are counters to identify the current basic and non-basic variables respectively
        # d is a row vector representing the current direction of improvement
        # into is a logical variable that indicates if the non-basic variable identified in the row vector: Crc, has already become a basic variable
        # out is a logical variable that indicates if the basic variable identified in the column vector: ratio, has already become a non-basic variable
        
        bv = 0
        nbv = 0
        d = vector(length = n)
        into = 0
        out = 0
        
        # iterate through each entry in the row vector: basis
        
        for(i in 1:n)
        {
            # if an entry in the current row vector: basis, has a value of zero then increment nbv
            # otherwise the entry has a value of one so increment bv
            
            if(basis[i] == 0)
            {
                nbv = nbv + 1
            } else
            {
                bv = bv + 1
                
                # assign the current entry in row vector: d, its corresponding value from column vector: db
                
                d[i] = db[bv,]
            }
                        
            # if nbv takes on a value corresponding to the non-basic variable identified in row vector: Crc, and it hasn't already been placed into the basis, then change this entry in row vector: basis, to a value of one to represent the correct non-basic variable becoming basic  
            
            if(into == 0)
            {
                if(type == "max")
                {
                    if(nbv == which(Crc == max(Crc))[1])
                    {
                        basis[i] = 1
                        into = 1
                        
                        # assign the current entry in row vector: d, a value of one because this corresponds to the non-basic variable becoming a basic variable
                        
                        d[i] = 1
                    }
                } else if(nbv == which(Crc == min(Crc))[1])
                {
                    basis[i] = 1
                    into = 1
                    
                    # assign the current entry in row vector: d, a value of one because this corresponds to the non-basic variable becoming a basic variable
                    
                    d[i] = 1
                }
            }			
            
            # if bv takes on a value corresponding to the basic variable identified in column vector: ratio, and it hasn't already been taken out of the basis, then change this entry in row vector: basis, to a value of zero to represent the correct basic variable becoming non-basic
            
            if(out == 0)
            {
                if(bv == which(ratio == lam)[1])
                {
                    basis[i] = 0
                    out = 1
                }
            }
        }
        
        # new extreme point
        
        x = x + lam*d
        
        # increment 'it', to track how many simplex iterations have occured
        
        it = it + 1
        
        # store data for this iteration
        
        tab = rbind(tab, cbind(it, as.numeric(C%*%x), lam, t(d), t(x)))
    }
    
    # ---- results ----
    # finalize the iteration table
    
    tab = data.frame(tab)
    dirs = vector(length = n)
    vars = vector(length = n)
    
    for(i in 1:n)
    {
        dirs[i] = paste0("d",i)
        vars[i] = paste0("x",i)
    }
    
    # defines the header names of the iteration table
    
    colnames(tab) = c("it", "obj", "step", t(dirs), t(vars))
    
    # returns the total number of iterations, the optimal objective value, the optimal solution, and an iteration table
    
    objval = as.numeric(C%*%x)
    BFS = data.frame(t(x))
    colnames(BFS) = t(vars)
    
    results = list("iterations" = as.numeric(it), "objective" = objval, "solution" = BFS, "table" = tab)
    return(results)
}

