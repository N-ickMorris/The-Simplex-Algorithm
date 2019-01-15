# --------------------------------------------------------------------------------------------------------------------
# continuous improvement search algorithm for linear programs
# --------------------------------------------------------------------------------------------------------------------

impS = function(type = "min", obj, A, equality, b, BFS, size = 1)
{
    # ---- setup ----

    require(pracma)
    optimal = FALSE
    it = 0
    tab = cbind(it, obj(BFS), 0, t(rep(0,length(BFS))), t(BFS))

    # ---- iterate through the Continuous Improvement Search Algorithm until it reaches a local optimal solution ----

    while(optimal == FALSE)
    {
        # ---- define the direction of improvement ----

        if(type == "min")
        {
		    # if the dot product of, the gradient and it's negative, is less than zero then there is an improvement direction
			
            if(sum(round((-grad(obj, BFS))*grad(obj, BFS),4)) < 0)
			{
			    d = round(-grad(obj, BFS),4)
			} else
			{
			    d = rep(0,length(BFS))
			}
        } else
        {
		    # if the dot product of, the gradient and itself, is greater than zero then there is an improvement direction
			
		    if(sum(round(grad(obj, BFS)*grad(obj, BFS),4)) > 0)
			{
			    d = round(grad(obj, BFS),4)
			} else
			{
			    d = rep(0,length(BFS))
			}
        }

        # ---- if there is no direction of improvment then we've reached a local optimal ----

        if(round(sum(abs(d)),4) == 0)
        {
            optimal = TRUE
        } else
        {
            # ---- define the largest step size for every constraint, as a vector (ie. lam) ----

            lam = vector(length = NROW(A))

            for(i in 1:NROW(A))
            {
                # if lam is left out of the constraint then the step size is Infinite

                if(round(sum(d*A[i,]),4) == 0)
                {
                    lam[i] = Inf  
                } else
                {
                    # if lam is included in the constraint then compute the step size

                    val = (b[i] - sum(BFS*A[i,])) / abs(sum(d*A[i,]))

                    if(equality[i] == "<=")
                    {
                        # if the constraint is a <= inequality then compute step size this way

                        if(sum(d*A[i,]) < 0)
                        {
                            lam[i] = Inf
                        } else if(sum(d*A[i,]) > 0 && round(val,4) > 0)
                        {
                            lam[i] = val
                        } else
                        {
                            lam[i] = 0
                        }
                    } else if(equality[i] == ">=")
                    {
                        # if the constraint is a >= inequality then compute step size this way

                        if(sum(d*A[i,]) < 0 && round(val,4) < 0)
                        {
                            lam[i] = abs(val)
                        } else if(sum(d*A[i,]) < 0)
                        {
                            lam[i] = 0
                        } else
                        {
                            lam[i] = Inf
                        }
					} else if(equality[i] == "<")
					{
                        # if the constraint is a < inequality then compute step size this way

                        if(sum(d*A[i,]) < 0)
                        {
                            lam[i] = Inf
                        } else if(sum(d*A[i,]) > 0 && round(val,4) > 0)
                        {
                            lam[i] = val - 0.1
                        } else
                        {
                            lam[i] = 0
                        }
				    } else if(equality[i] == ">")
					{
                        # if the constraint is a > inequality then compute step size this way

                        if(sum(d*A[i,]) < 0 && round(val,4) < 0)
                        {
                            lam[i] = abs(val) - 0.1
                        } else if(sum(d*A[i,]) < 0)
                        {
                            lam[i] = 0
                        } else
                        {
                            lam[i] = Inf
                        }
					} else
                    {
                        # the constraint is an equality, so compute step size this way

                        if(sum(d*A[i,]) < 0 && round(val,4) < 0)
                        {
                            lam[i] = abs(val)
                        } else if(sum(d*A[i,]) > 0 && round(val,4) > 0)
                        {
                            lam[i] = val
                        } else
                        {
                             lam[i] = 0
                        }
                    }
                }
            }

            # define the step size as the minimum value of lam multiplied by a scalar 'size'

            lam = round(min(lam) * size, 4)
            it = it + 1
			
			# ---- Check for Unboundedness ----
			
            if(lam == Inf)
		    {
				stop("Your model is Unbounded")
			}

            # ---- if there is a zero step size then we've reached a local optimal ----

            if(lam == 0)
            {
                optimal = TRUE

                # ---- store data for this iteration ----

                tab = rbind(tab, cbind(it, obj(BFS), lam, t(d), t(BFS)))

            } else
            {
                # ---- improve the BFS ----

				BFS = BFS + (lam*d)

                # ---- store data for this iteration ----

                tab = rbind(tab, cbind(it, obj(BFS), lam, t(d), t(BFS)))
            }
        }
    }

    # ---- results ----
    # finalize the iteration table

    tab = data.frame(tab)
    dirs = vector(length = length(BFS))
    vars = vector(length = length(BFS))

    for(i in 1:length(BFS))
    {
        dirs[i] = paste0("d",i)
        vars[i] = paste0("x",i)
    }

	# defines the header names of the iteration table
	
    colnames(tab) = c("it", "obj", "step", t(dirs), t(vars))

    # returns the total number of iterations, the local optimal objective value, the local optimal solution, and an iteration table

    objval = obj(BFS)
    BFS = data.frame(t(BFS))
    colnames(BFS) = t(vars)

    results = list("iterations" = it, "objective" = objval, "solution" = BFS, "table" = tab)
    return(results)
}

# --------------------------------------------------------------------------------------------------------------------
# HW 2 Solution
# --------------------------------------------------------------------------------------------------------------------

# ---- setup ----

obj_ = function(BFS)
{
    sin(BFS[2])*cos(BFS[1])*sin(BFS[1] + (BFS[2]^2))
}
A_ = rbind(c(2,1),c(-1,1),c(1,0),c(0,1))
sign_ = c("<=",">=",">=",">=")
b_ = c(25,-10,-2,0)
BFS_ = c(4.5, 13)
size_ = 1/3

# ---- solve ----

x = impS(type = "max", obj = obj_, A = A_, equality = sign_, b = b_, BFS = BFS_, size = size_)

# ---- results ----

x$iterations
x$objective
x$solution

# --------------------------------------------------------------------------------------------------------------------
# HW 2 Iteration Table
# --------------------------------------------------------------------------------------------------------------------

# This is a table with 9292 rows, R will not print out the entire table because of it's large size

x$table

# This table can be found in the excel spreedsheet provided in the tab 'Solution'
# OR you can create a blank notepad txt file and save it to your desktop, then use the following code...

write.table(x$table,file.choose())

# ...and click on the txt file that you created. This will export the table into your txt file.