from collections import Counter
from copy import deepcopy

eps = 1e-9

def termorder1(x,l = 1):
    # \prec_\ell as defined in the problem
    # tiebreak by graded reverse lexicographic term order
    return (sum(x[:l]), sum(x[:-1]),[-elem for elem in x[:-1][::-1]])

def termorder2(x):
    #lexicographic term order
    return (x[:-1])

def LT(f,termorder):
    f.sort(key=termorder)
    return f[-1][:]

def LTcounter(C,termorder):
    #turn counter to list
    f = [list(monomial)+[C[monomial]] for monomial in C if abs(C[monomial]) > eps] #remove all terms with a coefficient smaller than 1e-6

    if f:
        return LT(f,termorder)  
    else:
        return 0


def polynomgcd(f,g):
    #given two monomials, find their gcd
    assert(len(f) == len(g))

    gcd = [min(x,y) for x,y in zip(f,g)]
    return gcd

def divides(f,g):
    #given two monomials, check if f is divisible by g
    assert(len(f) == len(g))
    for x,y in zip(f[:-1],g[:-1]):
        if x >= y:
            continue
        return 0
    return 1

def NF(f,G,termorder):
    #compute the normal form of polynomial f w.r.t. G a set of polynomials, and a termorder

    C = Counter()

    for term in f:
        C[tuple(term[:-1])] = term[-1]
    

    #while C is nonempty, find the leading term, and check if any of the polynomials in G divides LT.
    while LTcounter(C,termorder):

        ltf = LTcounter(C,termorder)

        for polynomial in G:
            ltg = LT([term[:] for term in polynomial],termorder)
            if divides(ltf,ltg):
                q = [0]*len(ltf)

                for i in range(len(ltf)-1):
                    q[i] = ltf[i]-ltg[i]
                q[-1] = ltf[-1]/ltg[-1]


                #mult every term in polynomial with q, and subtract that from C
                tempg = [term[:] for term in polynomial]

                for term in tempg:
                    for i in range(len(ltf)-1):
                        term[i] += q[i]
                    term[-1] *= q[-1]

                for term in tempg:
                    C[tuple(term[:-1])] -= term[-1]

                #remove all empty polynomials from C
                newC = Counter()
                for term in C:
                    if abs(C[term]) < eps:
                        continue
                    newC[term] = C[term]

                newC,C = C,newC                

                break
        else:
            break

    
    for term in C:
        if abs(C[term]) < eps:
            continue

        polynomial = [list(monomial)+[C[monomial]] for monomial in C if abs(C[monomial]) > eps]
        return polynomial
    return 0



def s(f,g,termorder):
    #f,g are polynomials

    ltf = LT(f,termorder)
    ltg = LT(g,termorder)

    gcd = polynomgcd(ltf,ltg)

    for i in range(len(ltf)-1): #divide by the gcd
        ltf[i] -= gcd[i]
        ltg[i] -= gcd[i]
    
    C = Counter()

    
    for term in f:
        temp = term[:]
        for i in range(len(ltf)-1):
            temp[i] += ltg[i]

        C[tuple(temp[:-1])] += temp[-1]
    
    for term in g:
        temp = term[:]
        for i in range(len(ltf)-1):
            temp[i] += ltf[i]
        
        C[tuple(temp[:-1])] -= temp[-1]
    
    polynomial = [list(monomial)+[C[monomial]] for monomial in C if abs(C[monomial]) > eps] #turn counter to list
    return polynomial

def criterion_check(G,termorder):
    # Check if G fulfills the Buchbergers criterion
    # return True if criterion is fulfilled,
    # otherwise a pair (i,j), the pair of polynomials that it does not fulfill for.
    for i,f in enumerate(G):
        for j,g in enumerate(G):

            if NF(s(f,g,termorder),G,termorder) == 0:
                continue
            return (i,j)
    return True

def buchbergers(F, termorder):
    G = [[term[:] for term in polynomial] for polynomial in F] # G := F


    while 1:

        check = criterion_check(G,termorder)
        if check == True:
            break
        i,j = check

        f = [term[:] for term in G[i]]
        g = [term[:] for term in G[j]]

        h = NF(s(f,g,termorder),G,termorder)

        lth = LT(h,termorder)
        for term in h:
            term[-1] /= lth[-1]
        G.append(h)

    return G



def polynomial_printer(F,termorder,vars = "txyz"):
    #default is 4 varaibles, t,x,y,z

    for f in F:
        f.sort(key=termorder)
    
    
    for f in F:
        out = []
        for term in f:
            termstring = []
            if abs(term[-1]- (-1)) < eps:
                termstring.append("-")
            elif not (abs(term[-1]-1) < eps):
                temp = term[-1]
                if abs(temp-round(temp)) < eps:
                    temp = round(temp)

                termstring.append(str(temp))

            for i in range(len(term)-1):
                if term[i]:
                    if term[i] == 1:
                        termstring.append(vars[i])
                    else:
                        termstring.append(f"{vars[i]}^{term[i]}")
            out.append("".join(termstring))
    
        print(" + ".join(out)+",")



to = termorder1


# [t,x,y,z,coefficient]
F = [[[2,0,0,0,1],[0,2,0,0,1],[0,0,2,0,1],[0,0,0,2,1]], # t^2+x^2+y^2+z^2
     [[2,0,0,0,1],[0,2,0,0,2],[0,1,1,0,-1],[0,0,0,2,-1]], #t^2 + 2x^2-xy-z^2
     [[1,0,0,0,1],[0,0,3,0,1],[0,0,0,3,-1]]] #t + y^3-z^3



#print(NF(s(F[0],F[1],termorder1),F,termorder1),"result")

G = buchbergers(F,to)
print(len(G))
polynomial_printer(G,to)

F = [[[2,0,0,0,1],[0,2,0,0,1]], #t^2 + x^2
     [[1,1,0,0,1]]] #t+x


F = [[[2,0,0,0,1]], #p1^2 
     [[1,1,0,0,1], [0,2,0,0,1]]] #p1p2+p2^2

#print(buchbergers(F,termorder1))
