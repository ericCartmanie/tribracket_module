from sympy import *
Z = Symbol('Z')

#####################################
######## Knot and Link Code #########
#####################################
gknot = {
    (0):[[-1,1]],
(3,1):[[-1,2,-3,1,-2,3]],
(4,1):[[-1.5,2.5,-3,4,-2.5,1.5,-4,3]],
(5,1):[[-1,2,-3,4,-5,1,-2,3,-4,5]],
(5,2):[[-1,2,-3,4,-5,3,-2,1,-4,5]],
(6,1):[[-1,2,-3.5,4.5,-5.5,6.5,-2,1,-6.5,5.5,-4.5,3.5]],
(6,2):[[-1,2.5,-3.5,4.5,-5.5,1,-6,3.5,-4.5,5.5,-2.5,6]],
(6,3):[[1,-2.5,3.5,-4.5,2.5,-5,6,-1,5,-3.5,4.5,-6]],
#(6,3):[[-1,2,-3.5,4.5,-5.5,1,-6.5,3.5,-4.5,6.5,-2,5.5]],
(6,4):[[-1,2,-3,1,-2,3,-4,5,-6,4,-5,6]],
(6,5):[[-1,2,-3,1,-2,3,4.5,-5.5,6.5,-4.5,5.5,-6.5]],
(7,1):[[-1,2,-3,4,-5,6,-7,1,-2,3,-4,5,-6,7]],
(7,2):[[-1.5,6.5,-7.5,1.5,-2.5,3.5,-4.5,5.5,-6.5,7.5,-5.5,4.5,-3.5,2.5]],
(7,3):[[-1,2,-3,4,-7,6,-5,1,-2,3,-4,5,-6,7]],
(7,4):[[-1,2,-3,4,-5,6,-7,1,-4,3,-2,7,-6,5]],
(7,5):[[-1.5,2.5,-3.5,4.5,-5.5,1.5,-2.5,3.5,-6.5,7.5,-4.5,5.5,-7.5,6.5]],
(7,6):[[-1.5,2.5,-3,4,-5.5,6.5,-4,3,-7.5,1.5,-6.5,5.5,-2.5,7.5]],
(7,7):[[-1.5,2.5,-3,4,-2.5,5.5,-6,7,-5.5,1.5,-4,3,-7,6]],
(8,1):[[-1.5,2.5,-3,4,-2.5,1.5,-5.5,6.5,-7.5,8.5,-4,3,-8.5,7.5,-6.5,5.5]],
(8,2):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-3.5,4.5,-5.5,6.5,-7.5]],
(8,3):[[-1,2,-3,4,-5.5,6.5,-7.5,8.5,-4,3,-2,1,-8.5,7.5,-6.5,5.5]],
(8,4):[[-1.5,2,-3,4,-5,1.5,-6.5,7.5,-8.5,5,-4,3,-2,6.5,-7.5,8.5]],
(8,5):[[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,1.5,-6,7,-8,3,-4,5]],
(8,6):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-7.5,6.5,-3.5,4.5,-5.5]],
(8,7):[[-1.5,2.5,-3.5,1.5,-4,5,-6,7,-8,4,-2.5,3.5,-5,6,-7,8]],
(8,8):[[1,-2,3,-4.5,5.5,-6.5,4.5,-7,8,-5.5,6.5,-3,2,-1,7,-8]],
(8,9):[[-1.5,2.5,-3.5,4.5,-5,6,-7,8,-2.5,3.5,-4.5,1.5,-8,5,-6,7]],
#(8,10):[[-1,2,-3,4,-5,6,-7,1,-2,8.5,-6,7,-8.5,3,-4,5]],
(8,10):[[-1,2,-3,4,-5,1,-2,6.5,-7.5,3,-4,8.5,-6.5,7.5,-8.5,5]],
(8,11):[[1,-2,3.5,-4.5,5.5,-6.5,7.5,-8.5,2,-1,8.5,-5.5,4.5,-3.5,6.5,-7.5]],
(8,12):[[-1,2,-3,4,-5.5,6.5,-4,3,-7.5,8.5,-2,1,-8.5,7.5,-6.5,5.5]],
(8,13):[[1.5,-2.5,3,-4,5,-6,7,-3,8.5,-1.5,2.5,-8.5,4,-7,6,-5]],
(8,14):[[-1,2,-3.5,4.5,-5.5,6.5,-7.5,5.5,-4.5,8.5,-2,1,-6.5,7.5,-8.5,3.5]],
(8,15):[[1.5,-2.5,3.5,-4.5,5.5,-3.5,6.5,-7.5,8.5,-6.5,2.5,-1.5,7.5,-8.5,4.5,-5.5]],
(8,16):[[1.5,-2,3,-4.5,5.5,-6,2,-7.5,4.5,-5.5,8.5,-1.5,7.5,-3,6,-8.5]],
#(8,17):[[1.5,-2.5,3,-4,5.5,-1.5,6.5,-3,4,-7.5,8.5,-5.5,2.5,-6.5,7.5,-8.5]],
(8,17):[[-1.5,2.5,-3,4,-5,6,-2.5,7.5,-4,5,-8.5,1.5,-6,3,-7.5,8.5]],
(8,18):[[-1,2.5,-3.5,4,-5,6.5,-2.5,7,-4,8.5,-6.5,1,-7,3.5,-8.5,5]],
(8,19):[[-1,2,-3,-4,5,1,-2,-6,4,7,-8,-5,6,3,-7,8]],
(8,20):[[-1,2.5,3,-4,-5.5,1,6.5,-3,4,-7.5,8.5,5.5,-2.5,-6.5,7.5,-8.5]],
(8,21):[[-1,2.5,-3.5,-4.5,5.5,-6,-7.5,1,4.5,-5.5,8.5,7.5,-2.5,3.5,6,-8.5]],
(9,2):[[-1.5,2.5,-3.5,4.5,-5.5,6.5,-7.5,8.5,-9.5,1.5,-2.5,9.5,-8.5,7.5,-6.5,5.5,-4.5,3.5]],
(9,24):[[-1.5,2.5,-3,4.5,-5.5,6,-7,8,-2.5,1.5,-6,7,-8,9.5,-4.5,5.5,-9.5,3]],
(9,32):[[-1.5,2.5,-3,4,-5,6,-7,8,-2.5,9.5,-4,7,-8,3,-9.5,1.5,-6,5]],
(10,132):[[1.5,-2.5,3.5,-1.5,-4,5.5,6,-7.5,-8,4,9.5,-6,-10.5,8,2.5,-3.5,7.5,10.5,-5.5,-9.5]],
#(11,34):[[-1,2,-3,4.5,-5.5,6.5,-7.5,1,-8.5,9.5,-2,3,-10.5,8.5,-9.5,10.5,-11,7.5,-4.5,5.5,-6.5,11]],
(11,1):[[-1,2.5,-3.5,-4,5,6.5,-2.5,7,-8,3.5,-6.5,1,-7,8,9.5,-10.5,11.5,-5,4,-9.5,10.5,-11.5]],
(11,2):[[-1,2.5,-3.5,-4.5,5.5,-6.5,7,-8,4.5,-5.5,6.5,9.5,-2.5,10,-11,3.5,-9.5,1,-10,11,8,-7]],
#(11,34):[[-1,2,-3,4.5,5.5,-6.5,7.5,-5.5,8,-9,6.5,-7.5,-10.5,3,-11,-8,9,1,-2,10.5,-4.5,11]],
(11,42):[[-1,2,3.5,-4,5,-6,7,-3.5,8.5,-5,6,9.5,-10.5,11.5,-9.5,1,-2,10.5,-11.5,-7,4,-8.5]],
(-2,3,7):[[-1.5,-3,4,-5,6,-7,8,-9,10,-11,12,1.5,-2.5,-6,7,-8,9,-10,11,-12,3,-4,5,2.5]]
}

gknotlist=((0),(3,1),(4,1),(5,1),(5,2),(6,1),(6,2),(6,3),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7),(8,1),(8,2),(8,3),(8,4),(8,5),(8,6),(8,7),(8,8),(8,9),(8,10),(8,11),(8,12),(8,13),(8,14),(8,15),(8,16),(8,17),(8,18),(8,19),(8,20),(8,21))

glinklist=((2,0,1),(4,0,1),(5,0,1),(6,0,1),(6,0,2),(6,0,3),(6,0,4),(6,0,5),(6,1,1),(7,0,1),(7,0,2),(7,0,3),(7,0,4),(7,0,5),(7,0,6),(7,0,7),(7,1,1),(7,1,2))

glink = {
(0,2):[[-1,1],[-2,2]],
(0,3):[[-1,1],[-2,2],[-3,3]],
(2,0,1):[[-1,2],[-2,1]],
(4,0,1):[[-1,2,-3,4],[-4,3,-2,1]],
(5,0,1):[[-1,2.5,-3.5,4],[1,-4,5.5,-2.5,3.5,-5.5]],
(6,0,1):[[-1,2,-3,4],[1,-5.5,6.5,-4,3,-6.5,5.5,-2]],
#(6,0,2):[[-1,2,-5.5,4,-3.5,6],[1,-2,3.5,-4,5.5,-6]],
(6,0,2):[[-1,2,-3,4,-5,6],[1,-2,3,-6,5,-4]],
(6,0,3):[[-1,2,-3,4,-5,6],[1,-2,3,-4,5,-6]],
(6,0,4):[[-1,2,-3,4],[-5,1,-6,3],[-2,6,-4,5]],
#(6,0,5):[[-1,2.5,-3.5,4],[1,-5,6,-2.5],[3.5,-6,5,-4]],
(6,0,5):[[-1,2,-3,4],[1,-5,6,-2],[3,-6,5,-4]],
(6,1,1):[[-1,-2,3,4],[1,-5,-3,6],[2,5,-4,-6]],
(7,0,1):[[-1.5,2.5,-3,4,-5.5,1.5,-6,3,-7.5,5.5],[-2.5,7.5,-4,6]],
(7,0,2):[[-1,2,-3,4,-5,1,-2,5,-6,7],[3,-7,6,-4]],
(7,0,3):[[-1.5,2.5,-3.5,4.5,-5.5,1.5,-2.5,3.5,-6,7],[-4.5,5.5,-7,6]],
(7,0,4):[[-1.5,2.5,-3.5,4.5,-5.5,3.5,-2.5,1.5,-6,7],[-4.5,5.5,-7,6]],
(7,0,5):[[-1,2.5,-3.5,1,-4,5,-6,7],[-2.5,4,-7,6,-5,3.5]],
(7,0,6):[[-1,2.5,-3.5,4.5,-2.5,5,-6,7],[1,-7,6,-5,3.5,-4.5]],
(7,0,7):[[-1,2,-3,4],[1,-5.5,6.5,-2],[3,-7.5,5.5,-6.5,7.5,-4]],
(7,1,1):[[-1,2.5,3,-4.5,-5,1,6.5,-3,-7.5,5],[-2.5,-6.5,4.5,7.5]],
(7,1,2):[[-1,2.5,3,-4.5,-5,1,-6,-3,7,5],[-2.5,6,4.5,-7]]
#2.4:[[-1,4,-5,6,-3,2],[1,-2,3,-4,5,-6]],
#2.5:[[-1,2,-3,4,-5,6,-7,8],[1,-2,3,-8,7,-4,5,-6]],
#3.1:[[-1,2],[-2,1],[-3,3]],
#3.2:[[-1,2],[-2,3,-4,1],[-3,4]],
}

#####################################
## Helper Functions (tri-brackets) ##
#####################################
def tm(M):
### change the given matrix M into a tuple matrix
### this is implemented to tackle a quirk feature of python :)
# input: a matrix M
# output: M changed into a tuple matrix
    out  = []
    for x in M: out.append(tuple(x))
    return tuple(out)

def tmm(M):
### change the given matrix M into a tuple matrix
### this is implemented to tackle a quirk feature of python :)
# input: a matrix M
# output: M changed into a tuple matrix (deeper conversion than tm(M))

    out  = []
    for x in M: out.append(tm(x))
    return tuple(out)

def lm(M):
### change the given matrix M into a list matrix
### this is implemented to tackle a quirk feature of python :)
# input: a matrix M
# output: M changed into a list matrix
    out = []
    for x in M: out.append(list(x))
    return list(out)

def lmm(M):
### change the given matrix M into a list matrix
### this is implemented to tackle a quirk feature of python :)
# input: a matrix M
# output: M changed into a list matrix (deeper conversion than lm(M))

    out = []
    for x in M: out.append(lm(x))
    return list(out)

def zmatrix(n):
### return n by n Z matrix###
### this table is used as the initial tri-bracket table that is waited to be filled in
# input: x, y, n
# output: a n n* n matrix, with [b][a][c] = (ax - bxy + cy) %n
    out = []
    for i in range(0,n):
        r = []
        for j in range(0,n): r.append(Z)
        out.append(r)
    return out

############################
### alexander Tribracket ###
############################
def alexTriBracket(x,y,n):
# input: x, y and n
# output: a table of alexander tribracket matrix, in which [b,a,c] = ax - bxy + cy
    output = []
    for j in range(0,n):
        output.append(tm(zmatrix(n)))
    output = lmm(output)
    for b in range(1,n+1):
        for a in range(1,n+1):
            for c in range(1,n+1):
                output[b-1][a-1][c-1] = (x*a - x*y*b + c*y)%n
                if output[b-1][a-1][c-1]== 0:
                    output[b-1][a-1][c-1]=n
    return tmm(output)


######################################
### Helper Functions (tri-Modules) ###
######################################
def CTempty(n):
### CTempty refers to an empty co-efficient table all containing characters 'Z'
### this is the initial co-efficient table to be filled for either X_ _ _ or Y_ _ _
# input: modulo n
# output: n n*n matrices, with all entries being 'Z'
    firstIndex = []
    thirdIndex = []
    for x in range(0,n):
        thirdIndex.append(Z)
    for y in range(0,n):
        secondIndex = []
        for z in range(0,n): secondIndex.append(thirdIndex)
        firstIndex.append(secondIndex)
    return firstIndex

def findz(CT):
### findz takes in an empty coefficientTable, finds 1st instance of 'Z' in the CT, and returns its location
# input: Coefficient Table CT
# output: a tuple (x,y,z,l), with CT[x][y][z][l] being the location of the 1st instance of 'Z' in CT

    n = len(CT)
    for i in range(0,len(CT[0])):
        for j in range(0,len(CT[0][0])):
            for k in range(0,len(CT[0][0][0])):
                if CT[0][i][j][k] == Z:
                    return (0,i,j,k)
                if CT[1][i][j][k] == Z:
                    return (1,i,j,k)
    return False

def invs(x,n):
### find multiplicative inverse of x in modulo n;
### if such multiplicative inverse does not exist, a False is returned
# input: integer x and n
# output: multiplicative inverse of x in mod n, False if none exists
##    if gcd(x,n) != 1:
##        return False
#    else:
    for i in range(1,n):
        if (i*x) % n == 1:
            return i
    return False

##########################
### tri-bracket module ###
##########################
def triModFind(T,n):
### find all tri-bracket modules given a tribracket T and a system Z mod n
# input: tri-bracket T, Z mod n
# output: tribracket module based on T and n
    m = len(T)
    #creating a working list, which stores an empty co-effcient table for both X___ and Y___
    # an output list is created, and by default to be empty
    working, output = [[tmm(CTempty(m)), tmm(CTempty(m))]], []
    while len(working)!= 0:
        temporary = working[0]
        working[0:1] = []
        firstZ = findz(temporary)
        #if there is still unfilled entry left in the co-efficient tables
        if firstZ:
            for k in range(1,n):
                # obtain a list matrix version of original list working
                working2 = [lmm(temporary[0]), lmm(temporary[1])]
                # If k is invertible, then fill it, else continue the loop
                if invs(k,n):
                    #fill in k in firstZ in the co-efficient table
                    working2[firstZ[0]][firstZ[1]][firstZ[2]][firstZ[3]] = k
                    #fill in tables by checking all tri-bracket module axioms
                    working3 = triModFill2(working2,T,n)
                    #if the current table passes all axiom checking
                    if working3:
                        working.append(tuple([tmm(working3[0]), tmm(working3[1])]))
        #if all co-efficients have been filled up
        else:
            #attach them to the output
            output.append(tmm(temporary))
    return output


def triModFill(M,T,n):
### helper function for triModFind(T,n);
# input: a working co-efficients table M, a tri-bracket module T, and Z mod n
# output: same working table M, with axioms checked; if axioms fail to pass, a FALSE is returned
    m = len(T)
    X,Y = lmm(M[0]),lmm(M[1])
    for a in range(0,m):
        for b in range(0,m):
            for c in range(0,m):
                for d in range(0,m):
                    #here we start to check axioms accordingly;
# 2 axioms associated with u
                    # equation(1) associated with u: Xe*Xbad = Xf*Xbad - Xf*Yf + Yf*Xbac
                    # Equivalently, Yf*Xbac = Xe*Xbad - Xf*Xbad + Xf*Yf
                    # Equivalently, Xbac = (Yf)^-1 * (Xe*Xbad - Xf*Xbad + Xf*Yf)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                    if Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][a][d]!=Z and X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z:
                        RHS = (X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]) - (X[a][T[b][a][d]-1][T[b][a][c]-1] * X[b][a][d]) +(X[a][T[b][a][d]-1][T[b][a][c]-1]*Y[a][T[b][a][d]-1][T[b][a][c]-1])
                        if X[b][a][c] !=Z:
                            if ((Y[a][T[b][a][d]-1][T[b][a][c]-1] * X[b][a][c])%n) != (RHS%n):
                                return False
                        if X[b][a][c] == Z:
                            X[b][a][c] = (invs(Y[a][T[b][a][d]-1][T[b][a][c]-1],n) * RHS) % n

                    # equation(2) associated with u: Xe*Xbad  = Xg*Xbac
                    # Equivalently, Xbad = (Xe)^-1 * Xg * Xbac
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                    if X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][a][c]!=Z:
                        RHS = X[c][T[b][a][c]-1][T[b][d][c]-1]* X[b][a][c]
                        if X[b][a][d]!=Z:
                            if ((X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d])%n) != (RHS%n):
                                return False
                        if X[b][a][d] == Z:
                            X[b][a][d] = (invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n)* RHS) % n
# 2 axioms associated with v
                    # equation (1): - Xe*Xbad*Ybad - Ye*Xbdc*Ybdc = -Xf * Xbad*Ybad - Yf* Xbac * Ybac
                    # Equivalently, - Ye* Xbad * Ybdc = -Xf * Xbad*Ybad - Yf* Xbac * Ybac + Xe*Xbad*Ybad
                    # Equivalently, Ye* Xbad * Ybdc = Xf * Xbad*Ybad + Yf* Xbac * Ybac - Xe*Xbad*Ybad
                    # Equivalently, Ybdc = (Xbad)^-1 * (Ye)^-1 * (Xf * Xbad*Ybad + Yf* Xbac * Ybac - Xe*Xbad*Ybad)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                    if X[b][a][d]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and Y[b][a][d]!=Z and Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and X[b][a][c]!=Z and Y[b][a][c]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z:
                        RHS = (X[a][T[b][a][d]-1][T[b][a][c]-1]*X[b][a][d]*Y[b][a][d]) + (Y[a][T[b][a][d]-1][T[b][a][c]-1]*X[b][a][c]*Y[b][a][c]) - (X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]*Y[b][a][d])
                        if Y[b][d][c] != Z:
                            if (Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]*Y[b][d][c]%n) != (RHS%n):
                                return False
                        if Y[b][d][c] == Z:
                            Y[b][d][c] = (invs(X[b][a][d],n)* invs(Y[d][T[b][a][d]-1][T[b][d][c]-1],n)* RHS) % n

                    # equation (2): -Xe*Xbad*Ybad - Ye*Xbdc*Ybdc = - Xg*Xbac*Ybac - Yg*Xbdc*Ybdc
                    # Equivalently, -Xe*Xbad*Ybad = - Xg*Xbac*Ybac - Yg*Xbdc*Ybdc + Ye*Xbdc*Ybdc
                    # Equivalently, Xe*Xbad*Ybad = Xg*Xbac*Ybac + Yg*Xbdc*Ybdc - Ye*Xbdc*Ybdc
                    # Equivalently, Ybad = (Xbad)^-1 * (Xe)^-1 * (Xg*Xbac*Ybac + Yg*Xbdc*Ybdc - Ye*Xbdc*Ybdc)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                    # Yg = Y[c][T[b][a][c]][T[b][d][c]]
                    if X[b][a][d]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][a][c]!=Z and Y[b][a][c]!=Z and Y[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z and Y[b][d][c]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z:
                        RHS = (X[c][T[b][a][c]-1][T[b][d][c]-1] * X[b][a][c] * Y[b][a][c]) + (Y[c][T[b][a][c]-1][T[b][d][c]-1]*X[b][d][c]*Y[b][d][c]) - (Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c]*Y[b][d][c])
                        if Y[b][a][d]!=Z:
                            if ((X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]*Y[b][a][d])%n) != (RHS%n):
                                return False
                        if Y[b][a][d] == Z:
                            Y[b][a][d] = (invs(X[b][a][d],n)*invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n) * RHS) % n

# 2 axioms associated with w{
                    # Equation (1): Ye * Ybdc = Xg*Ybac - Xg * Yg + Yg*Ybdc
                    # Equivalently, Xg*Ybac = Ye * Ybdc + Xg * Yg - Yg*Ybdc
                    # Equivalently, Ybac = (Xg)^-1 * (Ye * Ybdc + Xg * Yg - Yg*Ybdc)
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                    # Yg = Y[c][T[b][a][c]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    if X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and Y[b][d][c]!=Z and Y[c][T[b][a][c]-1][T[b][d][c]-1]!=Z:
                        RHS = (Y[d][T[b][a][d]-1][T[b][d][c]-1]*Y[b][d][c]) + (X[c][T[b][a][c]-1][T[b][d][c]-1]*Y[c][T[b][a][c]-1][T[b][d][c]-1]) - (Y[c][T[b][a][c]-1][T[b][d][c]-1]*Y[b][d][c])
                        if Y[b][a][c]!=Z:
                            if ((X[c][T[b][a][c]-1][T[b][d][c]-1]*Y[b][d][c])%n) != (RHS % n):
                                return False
                        if Y[b][a][c] == Z:
                            Y[b][a][c] = (invs(X[c][T[b][a][c]-1][T[b][d][c]-1],n)* RHS) % n

                    # Equation (2): Ye*Ybdc = Yf*Ybac
                    # Equivalently, Ye = (Ybdc)^-1 * Yf * Ybac
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                    if Y[b][d][c]!=Z and Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and Y[b][a][c]!=Z:
                        RHS = Y[a][T[b][a][d]-1][T[b][a][c]-1] * Y[b][a][c]
                        if Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z:
                            if ((Y[d][T[b][a][d]-1][T[b][d][c]-1] * Y[b][d][c])%n) != (RHS%n):
                                return False
                        if Y[d][T[b][a][d]-1][T[b][d][c]-1]==Z:
                            Y[d][T[b][a][d]-1][T[b][d][c]-1] = (invs(Y[b][d][c],n)* RHS) % n

# 2 axioms associated with z{
                    # Equation (1): Xf*Ybad = Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Equivalently, Xf = (Ybad)^-1 * Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                    if Y[b][a][d]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z:
                        RHS = (X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[b][a][d]) - (X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[d][T[b][a][d]-1][T[b][d][c]-1]) + (Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c])
                        if X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z:
                            if ((X[a][T[b][a][d]-1][T[b][a][c]-1] * Y[b][a][d])%n) != (RHS%n):
                                return False
                        if X[a][T[b][a][d]-1][T[b][a][c]-1] == Z:
                            X[a][T[b][a][d]-1][T[b][a][c]-1] = (invs(Y[b][a][d],n)* RHS) % n
    return tuple([tmm(X),tmm(Y)])




def triModFill2(M,T,n):
### helper function for triModFind(T,n);
# input: a working co-efficients table M, a tri-bracket module T, and Z mod n
# output: same working table M, with axioms checked; if axioms fail to pass, a FALSE is returned
    m = len(T)
    X,Y = lmm(M[0]),lmm(M[1])
    keepgoing=True
    while keepgoing:
        keepgoing=False
        for a in range(0,m):
            for b in range(0,m):
                for c in range(0,m):
                    for d in range(0,m):
                    #here we start to check axioms accordingly;
# 2 axioms associated with u
                    # equation(1) associated with u: Xe*Xbad = Xf*Xbad - Xf*Yf + Yf*Xbac
                    # Equivalently, Yf*Xbac = Xe*Xbad - Xf*Xbad + Xf*Yf
                    # Equivalently, Xbac = (Yf)^-1 * (Xe*Xbad - Xf*Xbad + Xf*Yf)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                        if Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][a][d]!=Z and X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z:
                            RHS = (X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d] - X[a][T[b][a][d]-1][T[b][a][c]-1] * X[b][a][d] + X[a][T[b][a][d]-1][T[b][a][c]-1]*Y[a][T[b][a][d]-1][T[b][a][c]-1])%n
                            if not invs(RHS,n): return False #or not invs(Y[a][T[b][a][d]-1][T[b][a][c]-1],n): return False
                            if X[b][a][c] !=Z:
                                if ((Y[a][T[b][a][d]-1][T[b][a][c]-1] * X[b][a][c])%n) != (RHS%n):
                                    return False
                            if X[b][a][c] == Z:
                                X[b][a][c] = (invs(Y[a][T[b][a][d]-1][T[b][a][c]-1],n) * RHS) % n
                                keepgoing = True
                    # equation(2) associated with u: Xe*Xbad  = Xg*Xbac
                    # Equivalently, Xbad = (Xe)^-1 * Xg * Xbac
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                        if X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][a][c]!=Z:
                            RHS = (X[c][T[b][a][c]-1][T[b][d][c]-1]* X[b][a][c])%n
                            if not invs(RHS,n): return False # or not invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n): return False
                            if X[b][a][d]!=Z:
                                if ((X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d])%n) != (RHS%n):
                                    return False
                            if X[b][a][d] == Z:
                                X[b][a][d] = (invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n)* RHS) % n
                                keepgoing = True
# 2 axioms associated with v
                    # equation (1): - Xe*Xbad*Ybad - Ye*Xbdc*Ybdc = -Xf * Xbad*Ybad - Yf* Xbac * Ybac
                    # Equivalently, - Ye* Xbad * Ybdc = -Xf * Xbad*Ybad - Yf* Xbac * Ybac + Xe*Xbad*Ybad
                    # Equivalently, Ye* Xbad * Ybdc = Xf * Xbad*Ybad + Yf* Xbac * Ybac - Xe*Xbad*Ybad
                    # Equivalently, Ybdc = (Xbad)^-1 * (Ye)^-1 * (Xf * Xbad*Ybad + Yf* Xbac * Ybac - Xe*Xbad*Ybad)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                        if X[b][a][d]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and Y[b][a][d]!=Z and Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and X[b][a][c]!=Z and Y[b][a][c]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z:
                            RHS = (X[a][T[b][a][d]-1][T[b][a][c]-1]*X[b][a][d]*Y[b][a][d] + Y[a][T[b][a][d]-1][T[b][a][c]-1]*X[b][a][c]*Y[b][a][c] - X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]*Y[b][a][d])%n
                            if not invs(RHS,n): return False # or not invs(X[b][a][d],n) or not invs(Y[d][T[b][a][d]-1][T[b][d][c]-1],n): return False
                            if Y[b][d][c] != Z:
                                if (Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c]*Y[b][d][c]%n) != (RHS%n):
                                    return False
                            if Y[b][d][c] == Z:
                                Y[b][d][c] = (invs(X[b][d][c],n)* invs(Y[d][T[b][a][d]-1][T[b][d][c]-1],n)* RHS) % n
                                keepgoing = True
                    # equation (2): -Xe*Xbad*Ybad - Ye*Xbdc*Ybdc = - Xg*Xbac*Ybac - Yg*Xbdc*Ybdc
                    # Equivalently, -Xe*Xbad*Ybad = - Xg*Xbac*Ybac - Yg*Xbdc*Ybdc + Ye*Xbdc*Ybdc
                    # Equivalently, Xe*Xbad*Ybad = Xg*Xbac*Ybac + Yg*Xbdc*Ybdc - Ye*Xbdc*Ybdc
                    # Equivalently, Ybad = (Xbad)^-1 * (Xe)^-1 * (Xg*Xbac*Ybac + Yg*Xbdc*Ybdc - Ye*Xbdc*Ybdc)
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                    # Yg = Y[c][T[b][a][c]][T[b][d][c]]
                        if X[b][a][d]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][a][c]!=Z and Y[b][a][c]!=Z and Y[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z and Y[b][d][c]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z:
                            RHS = (X[c][T[b][a][c]-1][T[b][d][c]-1] * X[b][a][c] * Y[b][a][c] + Y[c][T[b][a][c]-1][T[b][d][c]-1]*X[b][d][c]*Y[b][d][c] - Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c]*Y[b][d][c])%n
                            if not invs(RHS,n): return False # or not invs(X[b][a][d],n) or not invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n): return False
                            if Y[b][a][d]!=Z:
                                if ((X[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][a][d]*Y[b][a][d])%n) != (RHS%n):
                                    return False
                            if Y[b][a][d] == Z:
                                Y[b][a][d] = (invs(X[b][a][d],n)*invs(X[d][T[b][a][d]-1][T[b][d][c]-1],n) * RHS) % n
                                keepgoing = True
# 2 axioms associated with w{
                    # Equation (1): Ye * Ybdc = Xg*Ybac - Xg * Yg + Yg*Ybdc
                    # Equivalently, Xg*Ybac = Ye * Ybdc + Xg * Yg - Yg*Ybdc
                    # Equivalently, Ybac = (Xg)^-1 * (Ye * Ybdc + Xg * Yg - Yg*Ybdc)
                    # Xg = X[c][T[b][a][c]][T[b][d][c]]
                    # Yg = Y[c][T[b][a][c]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                        if X[c][T[b][a][c]-1][T[b][d][c]-1]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and Y[b][d][c]!=Z and Y[c][T[b][a][c]-1][T[b][d][c]-1]!=Z:
                            RHS = (Y[d][T[b][a][d]-1][T[b][d][c]-1]*Y[b][d][c] + X[c][T[b][a][c]-1][T[b][d][c]-1]*Y[c][T[b][a][c]-1][T[b][d][c]-1] - Y[c][T[b][a][c]-1][T[b][d][c]-1]*Y[b][d][c])%n
                            if not invs(RHS,n): return False # or not invs(X[c][T[b][a][c]-1][T[b][d][c]-1],n): return False
                            if Y[b][a][c]!=Z:
                                if (X[c][T[b][a][c]-1][T[b][d][c]-1]*Y[b][a][c])%n != (RHS % n):
                                    return False
                            if Y[b][a][c] == Z:
                                Y[b][a][c] = (invs(X[c][T[b][a][c]-1][T[b][d][c]-1],n)* RHS) % n
                                keepgoing = True
                    # Equation (2): Ye*Ybdc = Yf*Ybac
                    # Equivalently, Ye = (Ybdc)^-1 * Yf * Ybac
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Yf = Y[a][T[b][a][d]][T[b][a][c]]
                        if Y[b][d][c]!=Z and Y[a][T[b][a][d]-1][T[b][a][c]-1]!=Z and Y[b][a][c]!=Z:
                            RHS = (Y[a][T[b][a][d]-1][T[b][a][c]-1] * Y[b][a][c])%n
                            if not invs(RHS,n): return False # or not invs(Y[b][d][c],n): return False
                            if Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z:
                                if ((Y[d][T[b][a][d]-1][T[b][d][c]-1] * Y[b][d][c])%n) != (RHS%n):
                                    return False
                            if Y[d][T[b][a][d]-1][T[b][d][c]-1]==Z:
                                Y[d][T[b][a][d]-1][T[b][d][c]-1] = (invs(Y[b][d][c],n)* RHS) % n
                                keepgoing = True
# 2 axioms associated with z{
                    # Equation (1): Xf*Ybad = Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Equivalently, Xf = (Ybad)^-1 * Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Xf = X[a][T[b][a][d]][T[b][a][c]]
                        if Y[b][a][d]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z:
                            RHS = (X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[b][a][d] - X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[d][T[b][a][d]-1][T[b][d][c]-1] + Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c])%n
                            if not invs(RHS,n): return False # or not invs(Y[b][a][d],n): return False
                            if X[a][T[b][a][d]-1][T[b][a][c]-1]!=Z:
                                if ((X[a][T[b][a][d]-1][T[b][a][c]-1] * Y[b][a][d])%n) != (RHS%n):
                                    return False
                            if X[a][T[b][a][d]-1][T[b][a][c]-1] == Z:
                                X[a][T[b][a][d]-1][T[b][a][c]-1] = (invs(Y[b][a][d],n)* RHS) % n
                                keepgoing = True
                    # Equation (2): Yg*Xbdc = Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Equivalently, Yg = (Xbdd)^-1 * Xe*Ybad - Xe*Ye +Ye*Xbdc
                    # Xe = X[d][T[b][a][d]][T[b][d][c]]
                    # Ye = Y[d][T[b][a][d]][T[b][d][c]]
                    # Yg = Y[c][T[b][a][c]][T[b][d][c]]
                        if Y[b][a][d]!=Z and X[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and Y[d][T[b][a][d]-1][T[b][d][c]-1]!=Z and X[b][d][c]!=Z:
                            RHS = (X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[b][a][d] - X[d][T[b][a][d]-1][T[b][d][c]-1]*Y[d][T[b][a][d]-1][T[b][d][c]-1] + Y[d][T[b][a][d]-1][T[b][d][c]-1]*X[b][d][c])%n
                            if not invs(RHS,n): return False # or not invs(X[b][d][c],n): return False
                            if Y[c][T[b][a][c]-1][T[b][d][c]-1]!=Z:
                                if ((Y[c][T[b][a][c]-1][T[b][d][c]-1] * X[b][d][c])%n) != (RHS%n):
                                    return False
                            if  Y[c][T[b][a][c]-1][T[b][d][c]-1]== Z:
                                Y[c][T[b][a][c]-1][T[b][d][c]-1] = (invs(X[b][d][c],n)* RHS) % n
                                keepgoing = True
    return tuple([tmm(X),tmm(Y)])



#############################################################
### color-counting enhancment based on Tribracket Modules ###
#############################################################
def gauss2pd(G):
#input: Gauss code G
#output: planar diagram of G
    out = []
    semiarccount = 0
    complen = []
    for i in range(0,len(G)):
        semiarccount = semiarccount + len(G[i])
        complen.append(len(G[i]))
    nx, cr = [], []
    for k in range(0,semiarccount):
        nx.append(0)
        cr.append(0)
    current = 1
    i = 0
    for x in range(1,len(G)+1):
        for y in range(1,len(G[x-1])+1):
            i = i+1
            if y == complen[x-1]:
                nx[i-1] = current
            else:
                nx[i-1] = i+1
            j = 0
            for z in range(1,len(G)+1):
                for w in range(1,len(G[z-1])+1):
                    j = j +1
                    if G[z-1][w-1] + G[x-1][y-1] == 0:
                        cr[i-1] = j
        current = current + complen[x-1]
    i = 0
    for x in range(1,len(G)+1):
        for y in range(1,len(G[x-1])+1):
            i = i+1
            if G[x-1][y-1] < 0:
                if G[x-1][y-1] % 1 == 0:
                    out.append([ 1,i,nx[cr[i-1]-1],nx[i-1],cr[i-1] ])
                else:
                    out.append([-1,i,cr[i-1],nx[i-1],nx[cr[i-1]-1] ])
    return out

def tpdhfindzero(x):
# helper function for countColors(G,X)
# input: a given Left,Right tuple-list x
# output: index of 1st encounterance of zero in that list, False if no zero detected
    for j in range(0,2):
        for k in range(0,len(x[0])):
            if x[j][k]==0:
                return([j,k])
    return False

##########################
###  primary coloring  ###
##########################
def regionColors(G,X):
# regionColor(G,X):
# input: a gauss code G and a tribracket X
# output: list of all valid colorings of G based on X

    # obtain planar diagram of given gauss code G
    D = gauss2pd(G)
    working, out, x = [],[],[[],[]]
    for j in range(0,2*len(D)):
        x[0].append(0)
        x[1].append(0)
    working = [tm(x)]
    while working != []:
        # take out first entry in the working list
        w = working[0]
        working[0:1] = []
        f = tpdhfindzero(w)
        # if there is no zero in the current working list
        if not f:
            # if all conditions are satisfied
            if colorFill(X,D,w):
                # append w to the output list
                out.append(w)
        # if there is still zero in the current working list
        else:
            # loop thorugh k with range(0,len(X[0]))
            for k in range(1,len(X[0])+1):
                w2 = lm(w)
                w2[f[0]][f[1]] = k
                w3 = colorFill(X,D,w2)
                if w3:
                    # if everything is good so far,append w3 to working
                    working.append(tm(w3))
    return out


def colorFill(X,D,w1):
# Helper function for count Colors, works with a specific list w
# input: Planar Diagram D, tribracket X and a given working lsit w
# output: all axioms checked and filled up accordingly
    Left,Right = list(w1[0]), list(w1[1])
    keepgoing = True
    while keepgoing:
        # keepgoing = False
        keepgoing = False
        # for each crossing:
        for x in D:
            # if we are dealing with a positive crossing:
            if x[0] == 1 or x[0] == 0:
                #check left region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Left[x[1]-1] == 0 and Left[x[2]-1]!=0:
                    Left[x[1]-1] = Left[x[2]-1]
                    keepgoing = True
                if Left[x[1]-1] !=0 and Left[x[2]-1]==0:
                    Left[x[2]-1] = Left[x[1]-1]
                    keepgoing = True
                # else if these 2 are not equal, return false
                if Left[x[1]-1] != Left[x[2]-1]:
                    return False

                # check right region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Right[x[3]-1] == 0 and Right[x[4]-1] !=0:
                    Right[x[3]-1] = Right[x[4]-1]
                    keepgoing = True
                if Right[x[3]-1]!=0 and Right[x[4]-1]==0:
                    Right[x[4]-1] = Right[x[3]-1]
                    keepgoing = True
                # else if these 2 are not equal, return false
                if Right[x[3]-1]!=Right[x[4]-1]:
                    return False

                # check Upper Region Match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Right[x[1]-1] == 0 and Left[x[4]-1] !=0:
                    Right[x[1]-1] = Left[x[4]-1]
                    keepgoing = True
                if Right[x[1]-1]!=0 and Left[x[4]-1]==0:
                    Left[x[4]-1] = Right[x[1]-1]
                    keepgoing = True
                # else if these 2 are not equal, return false
                if Right[x[1]-1]!=Left[x[4]-1]:
                    return False

                # check Lower Region Match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Right[x[2]-1] == 0 and Left[x[3]-1]!=0:
                    Right[x[2]-1] = Left[x[3]-1]
                    keepgoing = True
                if Right[x[2]-1] != 0 and Left[x[3]-1] == 0:
                    Left[x[3]-1] = Right[x[2]-1]
                    keepgoing = True
                # else if these 2 are not equal, return false
                if Right[x[2]-1] != Left[x[3]-1]:
                    return False

            # if we are dealing with a negative crossing:
            if x[0] == -1:
                #check left region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Left[x[2]-1] == 0 and Left[x[3]-1]!=0:
                    Left[x[2]-1] = Left[x[3]-1]
                    keepgoing = True
                if Left[x[2]-1]!=0 and Left[x[3]-1]==0:
                    Left[x[3]-1] = Left[x[2]-1]
                    keepgoing = True
                # else if these 2 are not equal, return False
                if Left[x[2]-1]!=Left[x[3]-1]:
                    return False

                # check right region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Right[x[1]-1] == 0 and Right[x[4]-1] !=0:
                    Right[x[1]-1] = Right[x[4]-1]
                    keepgoing = True
                if Right[x[1]-1] !=0 and Right[x[4]-1] == 0:
                    Right[x[4]-1] = Right[x[1]-1]
                    keepgoing = True
                # else if these 2 are not equal, return False
                if Right[x[1]-1] != Right[x[4]-1]:
                    return False

                # check for upper region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Left[x[1]-1] == 0 and Right[x[2]-1]!=0:
                    Left[x[1]-1] = Right[x[2]-1]
                    keepgoing = True
                if Left[x[1]-1]!=0 and Right[x[2]-1]==0:
                    Right[x[2]-1]=Left[x[1]-1]
                    keepgoing = True
                # else if these 2 are not equal, return False
                if Left[x[1]-1]!=Right[x[2]-1]:
                    return False

                # check for lower region match:
                # if one of 2 is empty, fill the empty one and set keepgoing as true
                if Left[x[4]-1]==0 and Right[x[3]-1]!=0:
                    Left[x[4]-1]=Right[x[3]-1]
                    keepgoing = True
                if Left[x[4]-1]!=0 and Right[x[3]-1]==0:
                    Right[x[3]-1] = Left[x[4]-1]
                    keepgoing = True
                # else if these 2 are not equal, return False
                if Left[x[4]-1]!=Right[x[3]-1]:
                    return False


# the next 2 blocks of code check in if the tri-bracket relationship has been satisfied

            # if we are dealing with a positive crossing
            if x[0] == 1:
                # if Left, Upper and lower regions have all been filled up
                b = Left[x[1]-1]
                a = Right[x[1]-1]
                c = Right[x[2]-1]
                if b !=0 and a !=0 and c != 0:
                    temp = X[b-1][a-1][c-1]
                    # if the right region has not been filled, filled with tribracket table
                    if Right[x[4]-1] == 0:
                        Right[x[4]-1] = temp
                        keepgoing = True
                    # if the right region has already been filled, check if tri-bracket relation gets satisfied
                    if Right[x[4]-1] != temp:
                        return False
            # if we are dealing with a negative crossing
            if x[0] == -1:
                # if Left, Upper and Lower regions have all been filled up
                b = Left[x[3]-1]
                a = Right[x[3]-1]
                c = Right[x[2]-1]
                if b!=0 and a!= 0 and c!=0:
                    temp = X[b-1][a-1][c-1]
                    # if the right region has not been filled, filled with tri-bracket table
                    if Right[x[4]-1] == 0:
                        Right[x[4]-1] = temp
                        keepgoing = True
                    # if the right region has already been filled, check if tri-bracket relation gets satisfied
                    if Right[x[4]-1] != temp:
                        return False
    return tuple([Left,Right])


##########################
### secondary coloring ###
##########################
def stickerColors(L,X,T,V,n):
# stickerColor(L,X,T,V,n):
# input: oriented-link L(in Gauss Code), alexander tribracket X, one X-coloring T, X-tribracket module V, scalar n
# output: compute sticker colorings of primarily-colored-G X by V

    # obtain planar diagram of given gauss code L
    D = gauss2pd(L)
    working, out, x = [],[],[[],[]]
    for j in range(0,2*len(D)):
        x[0].append(0)
        x[1].append(0)
    working = [tm(x)]
    while working != []:
        # take out 1st entry in the working list
        w = working[0]
        working[0:1] = []
        f = tpdhfindzero(w)
        
        # if there is no zero in the current working list
        if not f:
            # if all conditions are satisfied
            if stickerFill(w,D,X,T,V,n):
                # append w to the output list
                out.append(w)
        # if there is still zero in the current working list
        else:
            for k in range(1,n+1):
                w2 = lm(w)
                w2[f[0]][f[1]] = k
                w3 = stickerFill(w2,D,X,T,V,n)
                # if everything is good so far, append w3 to working
                if w3:
                    working.append(tm(w3))
    return out

def stickerFill(w,PD,X,T,V,n):
# Helper function for stickerColor
# input: 1 working list w, Planar Diagram D, alexander tribracket X, X-coloring T, X-tribracket module V, scalar n
# output: All axioms checked and filled up accordingly; False if w fails the check
    x,y = V[0],V[1]
    Left,Right = list(w[0]), list(w[1])
    keepgoing = True
    while keepgoing:
        # Set keepgoing False
        keepgoing = False
        # for each crossing
        for crossing in PD:
            # if we are dealing with a positive crossing
            if crossing[0] == 1 or crossing[0] == 0:
                # We firstly check region match
                # check left region match
                if Left[crossing[1]-1] == 0 and Left[crossing[2]-1] != 0:
                    Left[crossing[1]-1] = Left[crossing[2]-1]
                    keepgoing = True
                if Left[crossing[1] -1] !=0 and Left[crossing[2]-1] == 0:
                    Left[crossing[2]-1] = Left[crossing[1]-1]
                    keepgoing = True
                if Left[crossing[1]-1] != Left[crossing[2]-1]:
                    return False
                # check right region match
                if Right[crossing[4]-1] == 0 and Right[crossing[3]-1] !=0:
                    Right[crossing[4]-1] = Right[crossing[3]-1]
                    keepgoing = True
                if Right[crossing[4]-1] != 0 and Right[crossing[3]-1] == 0:
                    Right[crossing[3]-1] = Right[crossing[4]-1]
                    keepgoing = True
                if Right[crossing[4]-1] != Right[crossing[3]-1]:
                    return False
                # check upper region match
                if Right[crossing[1]-1] != 0 and Left[crossing[4]-1] == 0:
                    Left[crossing[4]-1] = Right[crossing[1]-1]
                    keepgoing = True
                if Right[crossing[1]-1] == 0 and Left[crossing[4]-1] != 0:
                    Right[crossing[1]-1] = Left[crossing[4]-1]
                    keepgoing = True
                if Right[crossing[1]-1] != Left[crossing[4]-1]:
                    return False
                # check lower region match
                if Right[crossing[2]-1]==0 and Left[crossing[3]-1] != 0:
                    Right[crossing[2]-1] = Left[crossing[3]-1]
                    keepgoing = True
                if Right[crossing[2]-1] != 0 and Left[crossing[3]-1] == 0:
                    Left[crossing[3]-1] = Right[crossing[2]-1]
                    keepgoing = True
                if Right[crossing[2]-1] != Left[crossing[3]-1]:
                    return False
            # if we are dealing with a negative crossing
            if crossing[0] == -1:
                # check Left Region match
                if Left[crossing[2]-1] == 0 and Left[crossing[3]-1]!=0:
                    Left[crossing[2]-1] = Left[crossing[3]-1]
                    keepgoing = True
                if Left[crossing[2]-1] != 0 and Left[crossing[3]-1] == 0:
                    Left[crossing[3]-1] = Left[crossing[2]-1]
                    keepgoing = True
                if Left[crossing[2]-1] != Left[crossing[3]-1]:
                    return False
                # check right Region match
                if Right[crossing[1]-1] == 0 and Right[crossing[4]-1]!=0:
                    Right[crossing[1]-1] = Right[crossing[4]-1]
                    keepgoing = True
                if Right[crossing[1]-1]!=0 and Right[crossing[4]-1]==0:
                    Right[crossing[4]-1] = Right[crossing[1]-1]
                    keepgoing = True
                if Right[crossing[1]-1]!=Right[crossing[4]-1]:
                    return False
                # check for upper region match
                if Left[crossing[1]-1] == 0 and Right[crossing[2]-1] != 0:
                    Left[crossing[1]-1] = Right[crossing[2]-1]
                    keepgoing = True
                if Left[crossing[1]-1]!=0 and Right[crossing[2]-1] == 0:
                    Right[crossing[2]-1] = Left[crossing[1]-1]
                    keepgoing = True
                if Left[crossing[1]-1]!=Right[crossing[2]-1]:
                    return False
                # check for lower region match
                if Left[crossing[4]-1] == 0 and Right[crossing[3]-1] != 0:
                    Left[crossing[4]-1] = Right[crossing[3]-1]
                    keepgoing = True
                if Left[crossing[4]-1] != 0 and Right[crossing[3]-1] == 0:
                    Right[crossing[3]-1] = Left[crossing[4]-1]
                    keepgoing = True
                if Left[crossing[4]-1] != Right[crossing[3]-1]:
                    return False
# The next 2 blocks of code check if the sticker coloring rules have been satisfied
            if crossing[0] == 1 or crossing[0]== 0:
                # load up primary region colors
                # T[0] is the left-region-label-list in primary coloring
                # T[1] is the right-region-label-list in primary coloring
                a = T[0][crossing[1]-1]
                b = T[0][crossing[4]-1]
                c = T[1][crossing[2]-1]
                # load up secondary colors
                v = Left[crossing[1]-1]
                u = Left[crossing[4]-1]
                w = Right[crossing[2]-1]
                z = Right[crossing[3]-1]
                if v !=0 and u!=0 and w!=0:
                    # load up the value of temp
                    temp = ((x[a-1][b-1][c-1] * u) - (x[a-1][b-1][c-1] * y[a-1][b-1][c-1] * v) + (y[a-1][b-1][c-1] * w)) % n
                    if temp == 0:
                        temp = n
                    # if z does not have a value yet, fill it with the axiom
                    if z == 0:
                        Right[crossing[3]-1] = temp
                        keepgoing = True
                    if Right[crossing[3]-1] != temp:
                        return False
        
            if crossing[0] == -1:
#                print("This is a negative crossing")
                # load up primary region colors
                # T[0] is the left-region-label-list in primary coloring
                # T[1] is the right-region-label-list in primary coloring
                c = T[1][crossing[2]-1]
                a = T[0][crossing[2]-1]
                b = T[1][crossing[3]-1]
                # load up secondary colors
                w = Right[crossing[2]-1]
                v = Left[crossing[2]-1]
                u = Right[crossing[3]-1]
                z = Right[crossing[4]-1]
                if v !=0 and u!=0 and w!=0:
                    # load up the value of temp
                    temp = ((x[a-1][b-1][c-1] * u) - (x[a-1][b-1][c-1] * y[a-1][b-1][c-1] * v) + (y[a-1][b-1][c-1] * w)) % n
                    if temp == 0:
                        temp = n
                    # if z does not have a value yet, fill it with the axiom
                    if z == 0:
                        Right[crossing[4]-1] = temp
                        keepgoing = True
                    if Right[crossing[4]-1] != temp:
                        return False
    return tuple([Left,Right])



##################
# input: oriented-link G(in Gauss Code), alexander tribracket X, X-tribracket module V, scalar n
def tbmodinv(G,X,V,n):
    ### compute tribracket module enhancement
    out,RC,u=0,regionColors(G,X),Symbol('u')
    for f in RC:
        out=out+u**len(stickerColors(G,X,f,V,n))
    return out
    
# input: oriented-link G(in Gauss Code), alexander tribracket X, X-tribracket module V, scalar n
def tbmodinv2(G,X,V,n):
    ### compute tribracket module enhancement
    out,RC,u=0,regionColors(G,X),Symbol('u')
    sum = len(RC)
    for f in RC:
        out=out+u**len(stickerColors(G,X,f,V,n))
    
    return (out, sum)

def zeromatrix(n):
    ### nxn matrix of zeroes
    out=[]
    for j in range (0,n):
        r=[]
        for k in range(0,n):
            r.append(0)
        out.append(r)
    return out


def findzero(N):
    ### find zero in tribracket
    for a in range(1,len(N)+1):
        for b in range(1,len(N)+1):
            for c in range(1,len(N)+1):
                if N[a-1][b-1][c-1]==0:
                    return [a,b,c]
    return False




def tensorrow(N,t,a,b):
    ### return tensor row of type t
    out=[]
    if t==0: # row
        for j in range(0,len(N)):
            out.append(N[a-1][b-1][j])
    if t==1: # column
        for j in range(0,len(N)):
            out.append(N[a-1][j][b-1])
    if t==2: # vertical column
        for j in range(0,len(N)):
            out.append(N[j][a-1][b-1])
    return out



def tbavail(N,a,b,c): 
    ### list available entries for N[a-1][b-1][c-1]
    out=list(range(1,len(N)+1))
    r1,r2,r3=tensorrow(N,0,a,b), tensorrow(N,1,a,c),tensorrow(N,2,b,c)
    for x in r1:
        if x in out:
            out.remove(x)
    for x in r2:
        if x in out:
            out.remove(x)
    for x in r3:
        if x in out:
            out.remove(x)
    return tuple(out)


def reptest(p): 
    ### Test whether p has nonzero repetitions
    for j in range(0,len(p)):
        for k in range(0,j): 
            if p[j]==p[k] and p[j]!=0:
                return True
    return False



def ntbcheck(N):
    ### test tribracket conditions
    for a in range(1,len(N)+1):
        for b in range(1,len(N)+1): 
            if reptest(tensorrow(N,0,a,b)): 
                return False
            if reptest(tensorrow(N,1,a,b)): 
                return False
            if reptest(tensorrow(N,2,a,b)): 
                return False
            for c in range(1,len(N)+1):
                for d in range(1,len(N)+1):
                    if N[a-1][b-1][N[b-1][c-1][d-1]-1]!=N[a-1][N[a-1][b-1][c-1]-1][N[N[a-1][b-1][c-1]-1][c-1][d-1]-1]:
                        return False
                    if N[N[a-1][b-1][c-1]-1][c-1][d-1]!=N[N[a-1][b-1][N[b-1][c-1][d-1]-1]-1][N[b-1][c-1][d-1]-1][d-1]:
                        return False
    return True


def ntbfill(N1):
    ### fill with horizontal tribracket conditions
    N=lmm(N1)
    keepgoing=True
    while keepgoing:
        keepgoing=False
        for a in range(1,len(N)+1):
            for b in range(1,len(N)+1): 
                if reptest(tensorrow(N,0,a,b)): 
                    return False
                if reptest(tensorrow(N,1,a,b)): 
                    return False
                if reptest(tensorrow(N,2,a,b)): 
                    return False
                for c in range(1,len(N)+1):
                    if N[a-1][b-1][c-1]==0:
                        x=tbavail(N,a,b,c)
                        if len(x)==1:
                            N[a-1][b-1][c-1]=x[0]
                            keepgoing=True
                    for d in range(1,len(N)+1):
                        if N[a-1][b-1][c-1]!=0 and N[a-1][b-1][d-1]!=0 and N[a-1][c-1][d-1]!=0:
                            if N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]==0 and N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]!=0:
                                N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]=N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]
                                keepgoing=True
                            if N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]!=0 and N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]==0:
                                N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]=N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]
                                keepgoing=True
                            if N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]==0 and N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]!=0:
                                N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]=N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]
                                keepgoing=True
                            if N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]!=0 and N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]==0:
                                N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]=N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]
                                keepgoing=True
                            if N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]==0 and N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]!=0:
                                N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]=N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]
                                keepgoing=True
                            if N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]!=0 and N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]==0:
                                N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1]=N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]
                                keepgoing=True
                            if N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]!=N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1] or N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]!=N[d-1][N[a-1][b-1][d-1]-1][N[a-1][c-1][d-1]-1] or N[c-1][N[a-1][b-1][c-1]-1][N[a-1][c-1][d-1]-1]!=N[b-1][N[a-1][b-1][c-1]-1][N[a-1][b-1][d-1]-1]:
                                return False
 #                       if N[b-1][c-1][d-1]!=0 and N[a-1][b-1][c-1]!=0:
 #                           if N[N[a-1][b-1][c-1]-1][c-1][d-1]!=0 and N[a-1][b-1][N[b-1][c-1][d-1]-1]!=0:
 #                               if N[a-1][N[a-1][b-1][c-1]-1][N[N[a-1][b-1][c-1]-1][c-1][d-1]-1]==0:
 #                                   N[a-1][N[a-1][b-1][c-1]-1][N[N[a-1][b-1][c-1]-1][c-1][d-1]-1]=N[a-1][b-1][N[b-1][c-1][d-1]-1]
 #                                   keepgoing=True
 #                               if N[a-1][b-1][N[b-1][c-1][d-1]-1]!=N[a-1][N[a-1][b-1][c-1]-1][N[N[a-1][b-1][c-1]-1][c-1][d-1]-1]:
 #                                   return False
 #                               if N[N[a-1][b-1][N[b-1][c-1][d-1]-1]-1][N[b-1][c-1][d-1]-1][d-1]==0:
 #                                   N[N[a-1][b-1][N[b-1][c-1][d-1]-1]-1][N[b-1][c-1][d-1]-1][d-1]=N[N[a-1][b-1][c-1]-1][c-1][d-1]
 #                                   keepgoing=True
 #                               if N[N[a-1][b-1][c-1]-1][c-1][d-1]!=N[N[a-1][b-1][N[b-1][c-1][d-1]-1]-1][N[b-1][c-1][d-1]-1][d-1]:
 #                                   return False
    return tmm(N)


def ntblist(n):
    ### list n-tribrackets
    N=[]
    for j in range(0,n):
        N.append(tm(zeromatrix(n)))
    working,out=[tmm(N)],[]
    while working!=[]:
        X=working[0]
        working[0:1]=[]
        f=findzero(X)
        if not f:
            if ntbcheck(X):
                out.append(tmm(X))
        if f:
            X1=lmm(X)
            L=tbavail(X,f[0],f[1],f[2])
            for j in L:
                X1[f[0]-1][f[1]-1][f[2]-1]=j
                X2=ntbfill(X1)
                if X2:
                    working.append(tmm(X2))
    return out

def ntblist3(r):
    ### list n-tribrackets
    n=len(r)
    N=[]
    for j in range(0,n):
        N.append(tm(zeromatrix(n)))
    N[0][0]=r
    working,out=[tmm(N)],[]
    while working!=[]:
        X=working[0]
        working[0:1]=[]
        f=findzero(X)
        if not f:
            if ntbcheck(X):
                out.append(tmm(X))
        if f:
            X1=lmm(X)
            L=tbavail(X,f[0],f[1],f[2])
            for j in L:
                X1[f[0]-1][f[1]-1][f[2]-1]=j
                X2=ntbfill(X1)
                if X2:
                    working.append(tmm(X2))
    return out

#########################
#########################
# input: A given tribracket module M
# output: tbmodinv(G,X,V,n) ran for every link L
def testing(X,M,n):
    print("_______________________________________________________________")
    print("\n")
    print("Test on Primary Color 1: ")
    g1s1 = [(2,0,1),(4,0,1),(5,0,1),(6,0,2),(7,0,4),(7,0,6)]
    g1s2 = [(7,0,2),(7,0,3),(7,1,1),(7,1,2)]
    g1s3 = [(6,0,1),(6,0,3),(7,0,1),(7,0,5)]
    print("\n")
    for L in g1s1:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))

    print("_______________________________________________________________")
    for L in g1s2:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))

    print("_______________________________________________________________")
    for L in g1s3:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))
    print("_______________________________________________________________")
    print("\n")
    print("Test on Primary Color 2: ")
    g2s1 = [(0,2)]
    g2s2 = [(6,0,5)]
    g2s3 = [(6,1,1),(7,0,7)]
    print("\n")
    for L in g2s1:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))
        
    print("_______________________________________________________________")
    for L in g2s2:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))
        
    print("_______________________________________________________________")
    for L in g2s3:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))
    print("_______________________________________________________________")
    print("\n")
    print("Test on Primary Color 3: ")
    g3s1 = [(6,0,4)]
    for L in g3s1:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))

    print("_______________________________________________________________")
    print("\n")
    print("Test on Primary Color 4: ")
    g4s1 = [(0,3)]
    for L in g4s1:
        poly = tbmodinv2(glink[L],X,M,n)[0]
        sum = tbmodinv2(glink[L],X,M,n)[1]
        print(str(L) +" is: "+">>>" + str(poly) + "; region colors: "+str(sum))




def tbmodtest(T,M,n):
    ### independently verify that M is a T-module mod n
    X,Y=M[0],M[1]
    for a in range(1,len(X)+1):
        for b in range(1,len(X)+1):
            for c in range(1,len(X)+1):
                for d in range(1,len(X)+1):
                    if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1])%n:
                        return 1#False
                    if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1])%n!=(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]+Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][d-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1])%n:
                        return 2#False
                    if (Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1])%n!=(Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][d-1])%n:
                        return 3#False
                    if (Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][d-1]+Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]-X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1])%n:
                        return 4#False
                    if (X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1])%n!=(Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1])%n:
                        return 5#False
                    if (X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1])%n!=(X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]-X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1])%n:
                        return 6#False
                    if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n!=(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1])%n:
                        return 7#False
                    if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1]+Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n:
                        return 8#False
    return True


def invtest(x,n):
    if x==Z: 
        return True
    for y in range(1,n):
        if (x*y)%n == 1: 
            return True
    return False
    


def tmfill(T,M,n):
    ### independently fill  T-module mod n
    X,Y=lmm(M[0]),lmm(M[1])
    keepgoing=True
    while keepgoing:
        keepgoing=False
        for a in range(1,len(X)+1):
            for b in range(1,len(X)+1):
                for c in range(1,len(X)+1):
                    if not invtest(X[a-1][b-1][c-1],n): 
                        return False
                    if not invtest(Y[a-1][b-1][c-1],n): 
                        return False
                    for d in range(1,len(X)+1):
                        if X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][b-1][c-1]!=Z  and X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z:
                            if X[a-1][b-1][d-1]!=Z:
                                if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1])%n:
                                    return False
                            if X[a-1][b-1][d-1]==Z:
                                X[a-1][b-1][d-1] =(X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*invs(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1],n))%n
                                keepgoing=True
# 2
                        if X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][b-1][c-1]!=Z and X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1] and X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1] !=Z and Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z:
                            if X[a-1][b-1][d-1]!=Z:
                                if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1])%n!=(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]+Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][d-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1])%n:
                                    return False
                            if X[a-1][b-1][d-1]==Z:
                                if not invtest((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]+X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1])%n,n):
                                    return False
                                X[a-1][b-1][d-1]=((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]+X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1])*invs(Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1],n))%n
                                keepgoing=True

# 3
                        if Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and Y[a-1][b-1][d-1]!=Z:
                            if Y[a-1][c-1][d-1]!=Z:
                                if (Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1])%n!=(Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][d-1])%n:
                                    return False
                            if Y[a-1][c-1][d-1]==Z:
                                Y[a-1][c-1][d-1]=(Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][d-1]*invs(Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1],n))%n
                                keepgoing=True
# 4

                        if Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[a-1][c-1][d-1]!=Z and X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[a-1][c-1][d-1]!=Z:
                            if Y[a-1][b-1][d-1]!=Z:
                                if (Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][d-1]+Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]-X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1])%n:
                                    return False
                            if Y[a-1][b-1][d-1]==Z:
                                if not invtest(Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]-Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]+X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1],n):
                                    return False
                                Y[a-1][b-1][d-1]=((Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]-Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][c-1][d-1]+X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1])*invs(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1],n))%n
#
                        if X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and Y[a-1][b-1][c-1]!=Z and Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z:
                            if X[a-1][c-1][d-1]!=Z:
                                if (X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1])%n!=(Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1])%n:
                                    return False
                            if X[a-1][c-1][d-1]==Z:
                                X[a-1][c-1][d-1]=(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1]*invs(Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1],n))%n
                                keepgoing=True
# 
                        if X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and Y[a-1][b-1][c-1]!=Z and X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[a-1][b-1][c-1]!=Z and Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z:
                            if X[a-1][c-1][d-1]!=Z:
                                if (X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1])%n!=(X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]-X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1])%n:
                                    return False
                            if X[a-1][c-1][d-1]==Z:
                                if not invtest(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1]-X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][c-1]
+X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1],n):
                                    return False
                                X[a-1][c-1][d-1]=((X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*Y[a-1][b-1][c-1]-X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[a-1][b-1][c-1]
+X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1])*invs(Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1],n))%n
                                keepgoing=True
#
                        if X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][b-1][c-1]!=Z and Y[a-1][b-1][c-1]!=Z and Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][c-1][d-1]!=Z and Y[a-1][c-1][d-1]!=Z and X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and X[a-1][b-1][c-1]!=Z and Y[a-1][b-1][c-1]!=Z and Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]!=Z and X[a-1][b-1][d-1]!=Z:
                            if Y[a-1][b-1][d-1]!=Z:
                                if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n!=(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1])%n:
                                    return False
                                if Y[a-1][b-1][d-1]==Z:
                                    if not invtest((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1])%n,n):
                                        return False
                                    Y[a-1][b-1][d-1]=((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1]-X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1])%n*(invx(X[b-1][T[a-1][b-1][c-1]-1][T[a-1][b-1][d-1]-1]*X[a-1][b-1][c-1],n)))%n
                                    keepgoing=True
#
                        if X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][b-1][c-1]!=Z and Y[a-1][b-1][c-1]!=Z and Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][c-1][d-1]!=Z and Y[a-1][c-1][d-1]!=Z and X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][b-1][d-1]!=Z and Y[a-1][b-1][d-1]!=Z and Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]!=Z and X[a-1][c-1][d-1]!=Z:
                            if Y[a-1][c-1][d-1]!=Z:
                                if (X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n!=(X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1]+Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1])%n:
                                    return False
                            if Y[a-1][c-1][d-1]==Z:
                                if not invtest((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1]-X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1])%n,n):
                                    return False
                                Y[a-1][c-1][d-1]=((X[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][c-1]*Y[a-1][b-1][c-1]+Y[c-1][T[a-1][b-1][c-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1]*Y[a-1][c-1][d-1]-X[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][b-1][d-1]*Y[a-1][b-1][d-1])*invs(Y[d-1][T[a-1][b-1][d-1]-1][T[a-1][c-1][d-1]-1]*X[a-1][c-1][d-1],n))%n
                                keepgoing=True                                
    return tuple([tmm(X),tmm(Y)])



def tmf(T,n):
### find all tri-bracket modules given a tribracket T and a system Z mod n
# input: tri-bracket T, Z mod n
# output: tribracket module based on T and n
    m = len(T)
    #creating a working list, which stores an empty co-effcient table for both X___ and Y___
    # an output list is created, and by default to be empty
    working, output = [[tmm(CTempty(m)), tmm(CTempty(m))]], []
    while len(working)!= 0:
        temporary = working[0]
        working[0:1] = []
        firstZ = findz(temporary)
        #if there is still unfilled entry left in the co-efficient tables
        if firstZ:
            for k in range(1,n):
                # obtain a list matrix version of original list working
                working2 = [lmm(temporary[0]), lmm(temporary[1])]
                # If k is invertible, then fill it, else continue the loop
                if invs(k,n):
                    #fill in k in firstZ in the co-efficient table
                    working2[firstZ[0]][firstZ[1]][firstZ[2]][firstZ[3]] = k
                    #fill in tables by checking all tri-bracket module axioms
                    working3 = tmfill(T,working2,n)
                    #if the current table passes all axiom checking
                    if working3:
                        working.append(tuple([tmm(working3[0]), tmm(working3[1])]))
        #if all co-efficients have been filled up
        else:
            #attach them to the output
            output.append(temporary)
    return output


def linkvals(X,M,n):
    #### values on links
    J=[]
    for x in glinklist:
        J.append([x,tbmodinv(glink[x],X,M,n)])
    return J


def knotvals(X,M,n):
    #### values on links
    J=[]
    for x in gknotlist:
        J.append([x,tbmodinv(gknot[x],X,M,n)])
    return J
