import sys
import numpy
load("BLISS.sage")

def dist(x,z):
    return  min(norm(x-z), norm(x+z))**2

# Generate m samples and compile them in a matrix.
def gen_A(s1,s2,az,m):
    toolbar_width = 50
    sys.stdout.write("-- Generate samples: [%s]\t" % (m))
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))

    r = []
    y = []
    som=0.;
    sigma=215.;

    while len(r)<m:
        c,z1,z2dag,z2,tx,t=sign(s1,s2,az,True)
        cadj=adjointR(c)
        w = vector(ZZ,(z1*cadj).list() +\
                      (z2dag*cadj).list())
        y.append((t))
        r.append(w)
        p=(1+exp(-sqrt(tx)/sigma^2))/2-1e-5
        cur=min(-t*log(1-p)-log(p),10);
        som+=cur;
        if len(r)%(m/50)==0:
            sys.stdout.write(u'\u2593')
            sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)+"\n")
    return len(y), vector(RDF, y), matrix(RDF, r)

# Spectral initialization for the descent.
def spectral_guess(A, y):
    npower_iter = 25
    toolbar_width = 50

    sys.stdout.write("-- Generate guess:\t\t")
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))

    T = RealDistribution('gaussian', 1)
    z0 = vector([T.get_random_element() for i in range(2*n)])
    z0 = z0/norm(z0)
    Y = (A.T*A)^(-1)
    for tt in range(1,npower_iter):
            #print(len(r)," ",som/len(r)," ",tx," ",p," ",cur," ",t)
        z0 = A.T*(diagonal_matrix(y)*(A*z0))
        z0 = Y*z0
        z0 = z0/norm(z0)
        sys.stdout.write(u'\u2593')
        sys.stdout.write(u'\u2593')
        sys.stdout.flush()
    z0 = z0/norm(z0)
    z = float(norm(s))*z0
    normest = norm(z)
    sys.stdout.write(u'\u2593')
    sys.stdout.write(u'\u2593')
    sys.stdout.write("\b" * (toolbar_width+1)+"\n")
    sys.stdout.write("-- Distance to secret: [%s]\n" % (dist(z,s)))
    return z

########################################################################
#                            Main routine                              #
########################################################################

if __name__=='__main__':
    # Read the number of samples from the standard input
    if len(sys.argv) < 2:
        sys.stdout.write("Error: please provide a number of samples\n")
        sys.exit(-1)
    nb_samples = int(sys.argv[1])

    # Generate secret key and merge the two secrets parts into one
    sys.stdout.write("-- Generate Key\n")
    s1,s2,a,az = s1gen()
    s = vector((vector(ZZ, s1).list()+vector(ZZ, s2).list()))

    # Generate the candidate secret
    m, y, A = gen_A(s1, s2, az, nb_samples)
    z = spectral_guess(A, y)

    # Save the datas as numpy binary dump and exit
    f = open("spc_guess_"+str(BLISS)+"_"+str(nb_samples),"w")
    numpy.savez(f, y=y, A=A, s=s, z=z);
    f.close()
