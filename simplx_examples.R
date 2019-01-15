# --------------------------------------------------------------------------------------------------------------------
# examples
# --------------------------------------------------------------------------------------------------------------------

# ---- maximization problem ----

# input information:

C = c(60,30,20,0,0,0,0)
A = rbind(c(8,6,1,1,0,0,0),c(4,2,1.5,0,1,0,0),c(2,1.5,0.5,0,0,1,0),c(0,1,0,0,0,0,1))
b = c(48,20,8,5)
x = c(0,0,0,48,20,8,5)
basis = c(0,0,0,1,1,1,1)

# solve:

simplx(type = "max", C = C, A = A, b = b, x = x, basis = basis)

# ---- minimization problem ----

# input information:

C = c(1,5,-2,0,0,0,0)
A = rbind(c(1,1,1,1,0,0,0),c(1,0,0,0,1,0,0),c(0,0,1,0,0,1,0),c(0,3,1,0,0,0,1))
b = c(4,2,3,6)
x = c(2,0,2,0,0,1,4)
basis = c(1,0,1,0,0,1,1)

# solve:

simplx(type = "min", C = C, A = A, b = b, x = x, basis = basis)

# ---- maximization problem (unbounded) ----

# input information:

C = c(3,2,1,0,0)
A = rbind(c(2,-3,2,1,0),c(-1,1,1,0,1))
b = c(3,5)
x = c(0,0,0,3,5)
basis = c(0,0,0,1,1)

# solve:

simplx(type = "max", C = C, A = A, b = b, x = x, basis = basis)

# ---- maximization problem ----

# input information:

C = c(1,-2,1,0,0,0)
A = rbind(c(1,-2,1,1,0,0),c(2,1,-1,0,1,0),c(-1,3,0,0,0,1))
b = c(12,6,9)
x = c(0,0,0,12,6,9)
basis = c(0,0,0,1,1,1)

# solve:

simplx(type = "max", C = C, A = A, b = b, x = x, basis = basis)

# --------------------------------------------------------------------------------------------------------------------
# step through
# --------------------------------------------------------------------------------------------------------------------

type = "max"
C = c(1,-2,1,0,0,0)
A = rbind(c(1,-2,1,1,0,0),c(2,1,-1,0,1,0),c(-1,3,0,0,0,1))
b = c(12,6,9)
x = c(0,0,0,12,6,9)
basis = c(0,0,0,1,1,1)
optimal = FALSE
it = 0
n = as.numeric(length(basis))
m = as.numeric(NROW(A))
tab = cbind(it, as.numeric(C%*%x), 0, t(rep(0,n)), t(x))

Cb = C[basis == 1]
Cn = C[basis == 0]
B = A[,basis == 1]
N = A[,basis == 0]
xb = x[basis == 1]
xn = x[basis == 0]
Crc = Cn - Cb%*%solve(B)%*%N

if(type == "max")
{
    db = -solve(B)%*%N[,Crc == max(Crc)]	
} else
{
    db = -solve(B)%*%N[,Crc == min(Crc)]
}

if(NCOL(db) > 1)
{
    db = matrix(db[,1])
}

ratio = -xb/ifelse(db >= 0, NA, db)
lam = min(ratio, na.rm = TRUE)

bv = 0
nbv = 0
d = vector(length = n)
into = 0
out = 0

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

x = x + lam*d
it = it + 1
tab = rbind(tab, cbind(it, as.numeric(C%*%x), lam, t(d), t(x)))
