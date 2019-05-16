"""
This is a sage8.3 script for the polynomial approximation (Section III.A)
of the paper submission.
For approximating the cosh, please uncomment the corresponding parameters.
"""

from itertools import izip

#COSH
"""
sigma=215
ss=2*sigma^2
gp.set("ss",ss)
a=-5534960
b=5534960
f = lambda x:  (exp(x/ss)+exp(-x/ss))/2
fdiff = lambda x: (exp(x/ss)-exp(-x/ss))/2/ss
invf = lambda x: 2/(exp(x/ss)+exp(-x/ss))
invfdiff = lambda x: 2/ss*(exp(-x/ss)-exp(x/ss))/(exp(x/ss)+exp(-x/ss))^2
gp("invf(x)= 2/(exp(x/ss)+exp(-x/ss))")
gp("invfdiff(x) =  2/ss*(exp(-x/ss)-exp(x/ss))/(exp(x/ss)+exp(-x/ss))^2")
gp("f(x) = (exp(x/ss)+exp(-x/ss))/2")
gp("fdiff(x) = (exp(x/ss)-exp(-x/ss))/2/ss")
target_secu = 62
prec=512
nbits=75
"""

# EXP

sigma=215
ss=2*sigma^2
a=-64082
b=0
invf = lambda x:  exp(-x/ss)
invfdiff = lambda x:-exp(-x/ss)/ss
f = lambda x: exp(x/ss)
fdiff = lambda x: exp(x/ss)/ss
gp("invf(x)= exp(-x/ss)")
gp("invfdiff(x) = -exp(-x/ss)/ss")
gp("f(x) = exp(x/ss)")
gp("fdiff(x) = exp(x/ss)/ss")
nbits=40
prec=128
target_secu = 40


#Other parameters
gp.set("ss",ss)
Fprec=ComplexBallField(prec)
R=PolynomialRing(Fprec, 'z')
gp.set_precision(prec)
gp.set("a",a)
gp.set("b",b)
C  = 2
Cs = b-a
gp.set("Cs",Cs)

#Functions

def getapproxcoeffs(P):
    coeffs=P.coefficients(sparse=False)
    return sum([Fprec(coeffs[i]).mid()*x^i for i in range(len(coeffs))])

#The following functions compute the Sobolev scalar product
# with GP integration or integral in the ComplexBallField Fprec

def scalH2cbf(u,v):
    w = u*v + diff(u)*diff(v)
    print w
    return Fprec.integral(lambda s, _: w(s),a,b)

def scalH2pol(u,v):
    gp.set("u(x)",getapproxcoeffs(u(x)))
    gp.set("v(x)",getapproxcoeffs(v(x)))
    du= lambda x: u(x).derivative()
    dv= lambda x: v(x).derivative()
    gp.set("du(x)",getapproxcoeffs(du(x)))
    gp.set("dv(x)",getapproxcoeffs(dv(x)))
    W=gp("W(x) = 1/Cs*u(x)*v(x) + Cs*du(x)*dv(x)")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val


def scalH2polscaled(u,v):
    du= lambda x: u(x).derivative()
    dv= lambda x: v(x).derivative()
    gp.set("u(x)",getapproxcoeffs(u(x)))
    gp.set("v(x)",getapproxcoeffs(v(x)))
    gp.set("du(x)",getapproxcoeffs(du(x)))
    gp.set("dv(x)",getapproxcoeffs(dv(x)))
    W=gp("W(x) = u(x)*v(x)*invf(x)^2/Cs +(u(x)*invfdiff(x)\
        +du(x)*invf(x))*(v(x)*invfdiff(x)+dv(x)*invf(x))*Cs")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val

def scalH2polscaled1(u,v):
    du= lambda x: u(x).derivative()
    dv= lambda x: v(x).derivative()
    gp.set("u(x)",getapproxcoeffs(u(x)))
    gp.set("v(x)",getapproxcoeffs(v(x)))
    gp.set("du(x)",getapproxcoeffs(du(x)))
    gp.set("dv(x)",getapproxcoeffs(dv(x)))
    W=gp("W(x) = u(x)*v(x)*invf(x)/Cs +(u(x)*invfdiff(x)\
        +du(x)*invf(x))*dv(x)*Cs")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val

def scalH2polf(pol,f,fdiff):
    gp.set("func(x)",getapproxcoeffs(f(x)))
    gp.set("pol(x)",getapproxcoeffs(pol(x)))
    gp.set("poldiff(x)",getapproxcoeffs(pol.derivative()(x)))
    gp.set("funcdiff(x)",getapproxcoeffs(fdiff(x)))
    W=gp("W(x) = pol(x)*func(x)/Cs+(u(x)*poldiff(x)*funcdiff(x))*Cs")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val

def scalH2polgscaled(u,g,gdiff):
    gp.set("u(x)",getapproxcoeffs(u(x)))
    gp.set("g(x)",getapproxcoeffs(g(x)))
    gp.set("gdiff(x)",getapproxcoeffs(gdiff(x)))
    du= u.derivative()
    gp.set("du(x)",getapproxcoeffs(du(x)))
    W=gp("W(x) = u(x)*g(x)*invf(x)^2/Cs+(u(x)*invfdiff(x) +\
        du(x)*invf(x))*(g(x)*invfdiff(x)+ gdiff(x)*invf(x))*Cs")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val


def scalH2polfscaled(u):
    "special case of scalH2polgscaled when g=f"
    gp.set("u(x)",getapproxcoeffs(u(x)))
    W=gp("W(x) = u(x)*invf(x)/Cs")
    val = Fprec(gp("intnum(x=a,b,W(x))"))
    return val


def findegree():
    """
    Computes the maximum degree of the projection that acheives
    the target security.
    Outputs the projection
    """
    def coefffj(fj):
        return scalH2polfscaled(fj)
    orthbasis = []
    orthcoeffs= []
    z = R.gen()
    i=0
    preci = -1
    print "iteration on the degrees :\n"
    if(f(x)==f(-x)):
        while(preci <= target_secu):
            basis = z^(2*i)
            projbasis = sum([scalH2polscaled(fj,basis)/cj * fj for (cj,fj) in \
                izip(orthcoeffs,orthbasis)])
            #print projbasis
            g = basis - projbasis
            orthbasis.append(g)
            orthcoeffs.append(scalH2polscaled(g,g))
            proj = sum([coefffj(fj)/cj * fj for (cj,fj) in izip(orthcoeffs,orthbasis)])
            preci = precisionproj(proj)
            print "degree %d -> precision %d" %(2*i,preci)
            i+=1
        i=2*i
    else :
        while(preci <= target_secu):
            basis = z^(i)
            projbasis = sum([scalH2polscaled(fj,basis)/cj * fj for (cj,fj) in \
                izip(orthcoeffs,orthbasis)])
            #print projbasis
            g = basis - projbasis
            orthbasis.append(g)
            orthcoeffs.append(scalH2polscaled(g,g))
            proj = sum([coefffj(fj)/cj * fj for (cj,fj) in izip(orthcoeffs,orthbasis)])
            preci = precisionproj(proj)
            print "degree %d -> precision %d" %(i,preci)
            i+=1
    print "Found a degree !\n"
    return preci, i, proj

def polyprojf():
    """
    Computes the orthogonal projection of g on the space of polynomials
    of degree < d, wrt the scaled H^2 inner product
    """
    orthbasis, orthcoeffs = polyorth()
    def coefffj(fj):
        return scalH2polfscaled(fj)
    return sum([coefffj(fj)/cj * fj for (cj,fj) in \
                        izip(orthcoeffs,orthbasis)])

def precision():
    """
    Computes b such that the scaled H^2 distance of f to polynomials of degree
    < d is at most 2^{-b} (roughly the number of relative bits of precision in
    replacing f by its projection polyprojf(d)).
    """
    g = polyprojf()
    gp.set("g(x)",getapproxcoeffs(g(x)))
    gp.set("dg(x)",getapproxcoeffs(g.derivative()(x)))
    W=gp("W(x) = (g(x)*invf(x)-1)^2/Cs +(g(x)*invfdiff(x)+dg(x)*invf(x))^2*Cs")
    dist2 = Fprec(gp("intnum(x=a,b,W(x))"))
    dist2 = dist2 * C
    return ceil((-log(dist2,2)/2).mid()),g

def precisionproj(g):
    """
    Computes b such that the scaled H^2 distance of f to polynomials of degree
    < d is at most 2^{-b} (roughly the number of relative bits of precision in
    replacing f by its projection polyprojf(d)).
    """
    gp.set("g(x)",getapproxcoeffs(g(x)))
    gp.set("dg(x)",getapproxcoeffs(g.derivative()(x)))
    W=gp("W(x) = (g(x)*invf(x)-1)^2/Cs +(g(x)*invfdiff(x)+dg(x)*invf(x))^2*Cs")
    dist2 = Fprec(gp("intnum(x=a,b,W(x))"))
    dist2 = dist2 * C
    return ceil((-log(dist2,2)/2).mid())

def gram_mat_H2(Basis,d):
    """
    Computes the gramm matrix with respect to the H2 norm
    """
    d=len(Basis)
    gram=[[scalH2polscaled(lambda x: Basis[i],lambda x:Basis[j]) \
            for j in range(d)] for i in range(d)]
    return matrix(gram)

def rationalize(elt):
    return Rational(real(elt.mid()))

def integer_mat(G,d):
    """
    Computes the integer matrix by multiplying by the denominators
    """
    d=G.nrows()
    dens=[]
    for i in range(d):
        for j in range(d):
            dens.append(rationalize(G[i,j]).denominator())
    denom=lcm(dens)
    return matrix(ZZ,denom*(G.apply_map(rationalize))), denom

def gen_lll_reduced_lattice(m,d):
    """
    Generates Lattice with the m_i and reduce it
    """
    z = R.gen()
    if(f(x)==f(-x)):
    	B = ([m[2*i]*z^(2*i) for i in range(floor(d/2))])
    else:
    	B = ([m[i]*z^i for i in range(d)])
    G = gram_mat_H2(B,d)
    G,denom= integer_mat(G,d)
    U = G.LLL_gram()
    U = matrix(R,U)
    return (U.T)*vector(B)#(U.T)*vector(B)/denom

def babai(L,float_proj,m,d):
    """
    Takes a reduced lattice in input and finds
    the closest vector to float_proj in the lattice L
    """
    z = R.gen()
    d=len(L)
    # Compute the Gram Schmidt vectors associated to L
    reduced_basis = L
    gram_vectors = []
    for i in range(d):
        gram_vectors.append(reduced_basis[i])
        for j in range(i):
            c = scalH2polscaled(lambda x: gram_vectors[i],lambda x:  gram_vectors[j])/ \
                scalH2polscaled(lambda x: gram_vectors[j], lambda x: gram_vectors[j])
            gram_vectors[i] -= c*gram_vectors[j]

    b = float_proj
    for i in range(d):
        c= scalH2polscaled(lambda x: b,lambda x: gram_vectors[d-1-i])/ \
            scalH2polscaled(lambda x: gram_vectors[d-1-i],lambda x:   \
            gram_vectors[d-1-i])
        c=round(c.mid())
        b=b-Fprec(c)*reduced_basis[d-1-i]
    return float_proj-b

def find_best_m(float_proj,d):
    m=[0 for i in range(d)]
    for i in range(d):
        if(f(x)==f(-x)):
            if i%2 == 0:
                m[i]= 2^(ceil(log(abs(float_proj.coefficients(sparse=false)[i].mid()\
                    *2^prec),2))-prec-nbits)
            else :
                m[i]=1
        else:
            m[i]= 2^(ceil(log(abs(float_proj.coefficients(sparse=false)[i].mid()\
                *2^prec),2))-prec-nbits)
    print "log(m[i]) = ",[log(m[i],2) for i in range(d)]
    return m

def approximate():
    print "...................................."
    print "Step 1 : look for a polynomial with float coefficients"
    print "that gives more than %d bits of secu\n"%target_secu
    preci,d,float_proj = findegree()
    print "Conclusion: Approximating f(x) by a degree %d polynomial ensures %s bits of precision\n" %(d-1,preci)
    print "The corresponding polynomial is: P=", getapproxcoeffs(float_proj)
    print "\n"
    print "...................................."
    print "Step 2 : Round the float coefficients to integer ones"
    print "Generating Lattice..."
    m=find_best_m(float_proj,d)
    L = gen_lll_reduced_lattice(m,d)
    print "LLL reduction Done\n"
    print "Babai Rounding...\n"
    rounded_proj = babai(L,float_proj,m,d)
    coeffs= [round(rounded_proj.coefficients(sparse=False)[i].mid()/m[i]) \
            for i in range(len(rounded_proj.coefficients(sparse=False)))]
    print "Max size of polynomial coefficients:", ceil(float(max(map(\
        lambda x: log(abs(x),2), coeffs))))

    diffP = lambda z: (rounded_proj(z) - f(z))*invf(z)
    derdiffP = lambda z: (rounded_proj.derivative()(z) - fdiff(z))*invf(z) +\
    (rounded_proj(z) - f(z))*invfdiff(z)
    gp.set("rounded_proj(x)",getapproxcoeffs(rounded_proj(x)))
    gp.set("diffrounded_proj(x)",getapproxcoeffs(rounded_proj.derivative()(x)))
    W=gp("W(x) = ( (rounded_proj(x)-f(x))*invf(x) )^2/Cs+\
        ( \
        (diffrounded_proj(x)-fdiff(x))*invf(x)\
     + (rounded_proj(x)-f(x))*invfdiff(x) \
     )^2*Cs")
    val = (Fprec(gp("intnum(x=a,b,W(x))"))*C).mid()
    print "final security :", ceil(-log(val,2)/2)

    print "....................................\n"
    #For writing them in latex
    print "latex exact values poly:\n P="
    for i in range(len(coeffs)) :
            print "&+ %s \cdot 2^{%s} \cdot x^{%s}\\\\ "%(coeffs[i],log(m[i],2),i)
    return coeffs

#approximate()
