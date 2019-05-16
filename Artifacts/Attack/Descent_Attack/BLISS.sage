from itertools import izip
from sage.stats.distributions.discrete_gaussian_integer \
        import DiscreteGaussianDistributionIntegerSampler

BLISS=3

if BLISS==0:
    #BLISS-0 parameters
    q=7681
    n=256
    (delta1,delta2)=(0.55,0.15)
    sigma=100
    kappa=12
    C = 1.5
    B2 = 2428
    BINF = 542
    d=5
elif BLISS==1:
    #BLISS-1 parameters
    q=12289
    n=512
    (delta1,delta2)=(0.3,0)
    sigma=215
    kappa=23
    C = 1.62
    B2 = 12872
    BINF = 2100
    d=10
elif BLISS==2:
    #BLISS-2 parameters
    q=12289
    n=512
    (delta1,delta2)=(0.3,0)
    sigma=107
    kappa=23
    C = 1.62
    B2 = 11074
    BINF = 1563
    d=10
elif BLISS==3:
    #BLISS-3 parameters
    q=12289
    n=512
    (delta1,delta2)=(0.42,0.03)
    sigma=250
    kappa=30
    C = 1.62
    B2 = 10206
    BINF = 1760
    d=9
elif BLISS==4:
    #BLISS-4 parameters
    q=12289
    n=512
    (delta1,delta2)=(0.45,0.06)
    sigma=271
    kappa=39
    C = 1.62
    B2 = 9901
    BINF = 1613
    d=8
else:
    raise NotImplemented

# Define working rings.
R.<xx>=QuotientRing(ZZ[x], ZZ[x].ideal(x^n+1))
Rq.<xxx>=QuotientRing(GF(q)[x], GF(q)[x].ideal(x^n+1))

sampler = DiscreteGaussianDistributionIntegerSampler(sigma=sigma, \
                                                     algorithm='uniform+table')
# Inner product of two polynomials, seens as vectors.
def inner(a,b):
    return sum([x*y for (x,y) in izip(a.list(), b.list())])

def ffnorm(a):
    return sum(map(lambda x: x*x, a))

def s1gen():
    s1vec=[0]*n
    s2vec=[0]*n
    d1=ceil(delta1*n)
    d2=ceil(delta2*n)

    # Gen S1
    while d1>0:
        i=randint(0,n-1)
        if s1vec[i]==0:
            s1vec[i]=(-1)^randint(0,1)
            d1-=1

    while d2>0:
        i=randint(0,n-1)
        if s1vec[i]==0:
            s1vec[i]=2*(-1)^randint(0,1)
            d2-=1

    d1=ceil(delta1*n)
    d2=ceil(delta2*n)

    # Gen S2
    while d1>0:
        i=randint(0,n-1)
        if s2vec[i]==0:
            s2vec[i]=(-1)^randint(0,1)
            d1-=1
    while d2>0:
        i=randint(0,n-1)
        if s2vec[i]==0:
            s2vec[i]=2*(-1)^randint(0,1)
            d2-=1

    s1 = sum([s1vec[i]*xxx^i for i in range(n)])
    s2 = sum([s2vec[i]*xxx^i for i in range(n)])
    s11 = sum([s1vec[i]*xx^i for i in range(n)])
    s22 = sum([s2vec[i]*xx^i for i in range(n)])
    s2 = 2*s2+1
    s22 = 2*s22+1

    # Rejection sampling ?
    t =sum([(inner(s1,xxx^i*s1)+inner(s1, xxx^i*s2))*xxx^i for i in
        range(n)])
    # WATCH OUT : T IS GIVN BY LINES (BY COLs IN ORIGINAL PAPER)
    T = ([((xxx^i*t).list()) for i in range(n)])
    # Now in good direction
    T = matrix(sorted(T, key= lambda li: ffnorm(li), reverse=True)).T
    sum_vec = [sum(sorted(li, reverse=True)[:kappa]) for li in T]
    NK = sum(sorted(sum_vec, reverse=True)[:kappa])
    d1=ceil(delta1*n)
    d2=ceil(delta2*n)
    if int(NK) >= int(C^2 *5* (d1+4*d2)*kappa):
        return s1gen()
    else:
        a = s2/s1
        a1zeta = (Rq(s22)/Rq(s11)).lift().change_ring(Integers(2*q))*1/(q-2)
        a1zeta = R(a1zeta.change_ring(ZZ))
        return (s11, s22, a, a1zeta)

def sign(s1, s2, a1zeta, sidechannel=True):
    reject=True

    while reject:
        y1=R([sampler() for _ in range(n)])
        y2=R([sampler() for _ in range(n)])
        u = a1zeta*y1.lift() + y2.lift()
        u = u.lift()

        #c is a random binary polynomial of weight kappa
        dc=kappa
        cvec=[0]*n
        while dc>0:
            i=randint(0,n-1)
            if cvec[i]==0:
                cvec[i]=1
                dc-=1

        c=R(cvec)
        z1=y1+c*s1
        z2=y2+c*s2

        # Rejection sampling
        scl = (s1*c).list() + (s2*c).list()
        zl  = z1.list() + z2.list()
        NSC = sum([x^2 for x in scl])
        ISC = sum([x*y for (x,y) in izip(scl,zl)])
        M = 3
        proba = 1/RR(M*exp(-NSC/(2*sigma^2))*cosh(ISC/(sigma^2)))
        #print proba
        if random() < proba:
            reject=False

    t=0
    coshreturnprob = RR((1 + exp(-abs(ISC)/sigma^2))/2)
    while True:
        if random() < coshreturnprob: break
        t += 1

    def roundD(x):
        x0 = x % (2^d)
        if x0 >= 2^(d-1):
            x0=x0-2^d
        return (x-x0)//(2^d)

    def roundDpol(f):
        return f.change_ring(ZZ).map_coefficients(roundD)

    p = ((2*q)//2^d)
    z2dag = (roundDpol(u) - roundDpol(u - z2.lift())) % ((2*q)//2^d)
    z2dag = z2dag.map_coefficients(lambda x: x if 2*x<p else x-p)

    if sidechannel:
        return (c,z1,z2dag,z2,ISC**2,t)
    else:
        return (c,z1,z2dag)

def verify(z1,z2):
    print  matrix(((z1).list() + (z2).list())).norm()
    return matrix(((z1).list() + (z2).list())).norm()< B2 and \
           abs(max(((z1).list() + (z2).list()))) < BINF

def adjointR(f):
    #return f.lift().subs(x=-xx^(n-1))
    return -xx*R(list(reversed(f.list())))
