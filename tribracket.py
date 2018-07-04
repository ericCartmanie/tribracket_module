from sympy import *
Z = Symbol('Z')
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
    if gcd(x,n) != 1:
        return False;
    else:
        for i in range(0,n):
            if (i*x) % n == 1:
                return i

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
            for k in range(0,n):
                # obtain a list matrix version of original list working
                working2 = [lmm(temporary[0]), lmm(temporary[1])]
                # If k is invertible, then fill it, else continue the loop
                if invs(k,n):
                    #fill in k in firstZ in the co-efficient table
                    working2[firstZ[0]][firstZ[1]][firstZ[2]][firstZ[3]] = k
                    #fill in tables by checking all tri-bracket module axioms
                    working3 = triModFill(working2,T,n)
                    #if the current table passes all axiom checking
                    if working3:
                        working.append((tmm(working3[0]), tmm(working3[1])))
        #if all co-efficients have been filled up
        else:
            #attach them to the output
            output.append(temporary)
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
    return [tmm(X),tmm(Y)]


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
def listColors(G,X):
# listColors(G,X):
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


























