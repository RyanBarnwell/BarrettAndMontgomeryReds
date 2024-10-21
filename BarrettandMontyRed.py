#barrett Reduction

import math
def barrett_red(a, b, q):
    #k calculation
    k = math.ceil(math.log(q, 2))
    
    #r = 2^k
    r = 2
    for x in range(1, k):
        r = r*2
    #mu
    mu = int((r*r)/q)
    
    z = a*b
    
    m1 = int(z/r)
    m2 = m1 * mu
    m3 = int(m2/r)
    t = z - (m3*q)
    if(t>=q):
        return t-q
    else:
        return t
print(barrett_red(1467, 2489, 7681))
print(barrett_red(500, 7, 7681))

#Montgomery Reduction
def montgomery_red(a, b, q):
    k = math.ceil(math.log(q, 2))
    #r = 2^k
    r = 2
    for x in range(1, k):
        r = r*2
    rinv = pow(r, -1, q)
    qPrime = pow(-q, -1, r)
    #convert to montgomery representation
    am = (a*r) % q
    bm = (b*r) % q
    
    t= am * bm
    u= t * qPrime % r
    cm = int((t+u*q)/r)
    
    #convert back
    c = cm * rinv % q
    return c
    
print(montgomery_red(1467, 2489, 7681))
print(montgomery_red(500, 7, 7681))

#both functions return 2888 for the first set of inputs which is known to be correct.  The second input verifies that the functions work when using the else side of
#if-else statement in the barret reduction