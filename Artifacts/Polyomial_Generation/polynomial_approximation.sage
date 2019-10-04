"""
This is a sage8.3 script for the polynomial approximation (Section III.A)
of the paper submission.
"""

from itertools import izip
from sage.all import *


class GALACTICS_approximation():

    def __init__(self,
                 a,
                 b,
                 prec,
                 target_precision,
                 nbits,
                 f,
                 fdiff=None,
                 invf=None,
                 invfdiff=None):
        """ Constructor for setting up a GALACTICS approximation framework.
        :a: min of the interval approximation
        :b: max of the interval approximation
        :prec: bits of precision for the internal computations
        :target_precision: target integer K for
                    |f-P|/|f| <= 2^-K
        :nbits: target number of bits of the polynomial P
        :f: the function of x, written as a string (for the direct use of gp())
        The following variables are optionnal but it is usually
        faster to provide them
        :fdiff: the differential of f, written as a string
        :invf: the inverse of f, written as a string
        :invfdiff: the differential of the inverse of f, written as a string
        """
        self.a = a
        self.b = b
        gp.set("a", a)
        gp.set("b", b)
        self.C = 2
        self.Cs = b - a
        gp.set("Cs", self.Cs)
        self.prec = prec
        gp.set_precision(prec)
        self.Fprec = ComplexBallField(prec)
        self.R = PolynomialRing(self.Fprec, 'z')
        self.target_precision = target_precision
        self.nbits = nbits
        gp("f(x) = " + f)
        if fdiff is None:
            gp("fdiff(x) = deriv(f)(x)")
        else:
            gp("fdiff(x) =" + fdiff)
        if invf is None:
            gp("invf(x) = 1 / f(x)")
        else:
            gp("invf(x) = " + invf)
        if invfdiff is None:
            gp("invfdiff(x) = deriv(invf)(x)")
        else:
            gp("invfdiff(x) =" + invfdiff)

    def get_coeffs(self, P):
        """ Reformats a polynomial with coefficients in FPrec
        into a polynomial with coefficients rounded to prec-bit precision.
        Used for gp(.)
        :P: a polynomial in R
        """
        coeffs = P.coefficients(sparse=False)
        return sum([self.Fprec(coeffs[i]).squash().mid() * x ^ i
                    for i in range(len(coeffs))])

    def sobolev_scalar_product(self, u, v):
        """ Computes the sobolev scalar product <u,v> between two
        functions u and v. This scalar product is scaled with the
        function f.
        :u: and :v: derivable function of x
        """
        def du(x):
            return u(x).derivative()

        def dv(x):
            return v(x).derivative()
        gp.set("u(x)", self.get_coeffs(u(x)))
        gp.set("v(x)", self.get_coeffs(v(x)))
        gp.set("du(x)", self.get_coeffs(du(x)))
        gp.set("dv(x)", self.get_coeffs(dv(x)))
        gp("W(x) = u(x) * v(x) * invf(x) ^ 2 / Cs +(u(x) * invfdiff(x)\
            + du(x) * invf(x)) * (v(x) * invfdiff(x) + dv(x) * invf(x))\
             * Cs")
        val = self.Fprec(gp("intnum(x = a, b, W(x))"))
        return val

    def sobolev_scalar_product_with_f(self, u):
        """ Computes the scaled sobolev scalar product <u,f>
        between u and the function f
        :u: a derivable function
        """
        gp.set("u(x)", self.get_coeffs(u(x)))

        def du(x):
            return u(x).derivative()
        gp.set("du(x)", self.get_coeffs(du(x)))
        gp("W(x) = u(x) * invf(x) / Cs + \
            (u(x) * invfdiff(x) + du(x) * invf(x)) \
            * (f(x) * invfdiff(x) + fdiff(x) * invf(x))\
             * Cs")
        val = self.Fprec(gp("intnum(x=a, b, W(x))"))
        return val

    def find_degree(self, force_degree=None):
        """
        Computes the maximum degree of the projection that acheives
        the target precision. It proceeds by iterating the degree of
        the projection space.
        :target_precision: variable defining the target precision
        Outputs the precision, the degree and the projection
        """
        z = self.R.gen()
        ortho_basis = []
        ortho_coeffs = []
        precision = -1
        projection = 0
        print "iteration on the degrees :\n"
        i = 0
        while(precision <= self.target_precision or (force_degree is not None and i <= force_degree)):
            can_basis = z ^ i
            projection_on_basis = sum([self.sobolev_scalar_product(b,
                                                                   can_basis) /
                                       c * b for (c, b) in
                                       izip(ortho_coeffs, ortho_basis)])
            g = can_basis - projection_on_basis
            norm_g = self.sobolev_scalar_product(g, g)
            ortho_basis.append(g)
            ortho_coeffs.append(norm_g)
            projection += self.sobolev_scalar_product_with_f(g) / norm_g * g
            precision = self.precision_proj(projection)
            print("degree %d -> precision %d" % (i, precision))
            if (gp("f(x) == f(- x)")):
                i += 2
            else:
                i += 1
        print "Found a degree !\n"
        return precision, i, projection

    def precision_proj(self, g):
        """
        Computes the relative distance between g and f. Then, finds b such that
        the relative distance is at most 2 ^ {-b}
        :g: a derivable function
        """
        gp.set("g(x)", self.get_coeffs(g(x)))
        gp.set("dg(x)", self.get_coeffs(g.derivative()(x)))
        gp("W(x) = (g(x) * invf(x) - 1) ^ 2 / Cs + (g(x) * invfdiff(x) + \
            dg(x) * invf(x)) ^ 2 * Cs")
        d = self.Fprec(gp("intnum(x = a, b, W(x))"))
        d = d * self.C
        return ceil((- log(d, 2) / 2).squash().mid())

    def compute_gram_matrix(self, Basis):
        """
        Computes the gram matrix with respect to the Sobolev norm
        Basis: a d * d array
        """
        d = len(Basis)
        gram = [[self.sobolev_scalar_product(lambda x: Basis[i],
                                             lambda x: Basis[j])
                for j in range(d)] for i in range(d)]
        return matrix(gram)

    def rationalize(self, elt):
        """
        Reformats an element in Fprec into a rational element.
        :elt: an element in Fprec
        """
        return Rational(real(elt.squash().mid()))

    def integer_mat(self, G):
        """
        Reformats a matrix in Fprec into an integer matrix
        by transforming the elements to rationals and
        multiplying by the denominators.
        :G: an matrix in Fprec
        """
        d = G.nrows()
        dens = []
        for i in range(d):
            for j in range(d):
                dens.append(self.rationalize(G[i, j]).denominator())
        denom = lcm(dens)
        return matrix(ZZ, denom * (G.apply_map(self.rationalize))), denom

    def gen_lll_reduced_lattice(self, m):
        """
        Generates Lattice with the m_i and reduce it
        :m: array of rationals
        """
        d = len(m)
        z = self.R.gen()
        if(gp("f(x) == f(- x)")):
            B = ([m[2 * i] * z ^ (2 * i) for i in range(floor(d / 2))])
        else:
            B = ([m[i] * z ^ i for i in range(d)])
        G = self.compute_gram_matrix(B)
        G, denom = self.integer_mat(G)
        U = G.LLL_gram()
        U = matrix(self.R, U)
        return (U.T) * vector(B)

    def babai_reduce(self, L, float_proj):
        """
        Finds a vector in the lattice L close to float_proj using
        Babai's alogrithm.
        :L: a LLL reduced lattice represented by a basis matrix.
        :float_proj: an element in R
        """
        d = len(L)
        reduced_basis = L
        gram_vectors = []
        for i in range(d):
            gram_vectors.append(reduced_basis[i])
            for j in range(i):
                c = (self.sobolev_scalar_product(lambda x: gram_vectors[i],
                                                 lambda x: gram_vectors[j]) /
                     self.sobolev_scalar_product(lambda x: gram_vectors[j],
                                                 lambda x: gram_vectors[j]))
                gram_vectors[i] -= c * gram_vectors[j]
        b = float_proj
        for i in range(d):
            c = (self.sobolev_scalar_product(
                lambda x: b, lambda x: gram_vectors[d - 1 - i]) /
                self.sobolev_scalar_product(
                lambda x: gram_vectors[d - 1 - i],
                lambda x: gram_vectors[d - 1 - i]))
            c = round(c.squash().mid())
            b = b - self.Fprec(c) * reduced_basis[d - 1 - i]
        return float_proj - b

    def find_best_m(self, float_proj, d):
        """
        Sets the number of bits of the first coefficient to
        nbits.
        and computes the other m coefficients according to the
        coefficients float_proj
        :float_proj: an element in R
        :d: degree of float_proj
        """
        m = [0 for i in range(d)]
        for i in range(d):
            if (gp("f(x) == f(- x)")):
                if i % 2 == 0:
                    m[i] = 2 ^ (ceil(
                        log(abs(
                            float_proj.coefficients(sparse=false)[i].squash().mid() *
                                2 ^ self.prec), 2)) - self.prec - self.nbits)
                else:
                    m[i] = 1
            else:
                m[i] = 2 ^ (ceil(
                    log(abs(float_proj.coefficients(sparse=false)[i].squash().mid() *
                            2 ^ self.prec), 2)) - self.prec - self.nbits)
        print "log(m[i]) = ", [log(m[i], 2) for i in range(d)]
        return m

    def integer_coefficients(self, poly, m):
        """
        Returns an array containing the coefficients of the polynomial
        divided by the elements m
        :poly: element in R
        :m: array of rational
        """
        return [long(poly.coefficients(sparse=False)[i].squash().mid().real() / m[i])
                for i in range(len(poly.coefficients(sparse=False)))]

    def update_m(self, coeffs, m):
        """Remove trailing zeros in the binary writting of coeffs
        by adding the corresponding power of two in m
        :coeffs: array of coefficients (long)
        :m: array of power of twos
        """
        new_coeffs = coeffs
        new_m = m
        for i in range(len(coeffs)):
            coeff = coeffs[i]
            if coeff != 0.:
                j = 0
                while (coeff >> j) % 2 == 0:
                    j += 1
                new_m[i] *= 2 ** j
                new_coeffs[i] = new_coeffs[i] >> j
        return new_coeffs, new_m

    def approximate(self, force_degree=None):
        """Main function that runs the approximation.
        It finishes by printing the latex writing of
        the projection and outputs the coefficients
        of the polynomial and the coefficients m.
        """
        print("....................................")
        print("Step 1 : look for a polynomial with float coefficients")
        print("that gives more than %i bits of secu\n" % target_precision)
        preci, d, float_proj = self.find_degree(force_degree)
        print("Conclusion: Approximating f(x) by a degree %i polynomial ensures %i \
        bits of precision\n" % (d - 1, preci))
        print("The corresponding polynomial is: P = ")
        print(self.get_coeffs(float_proj))
        print("\n")
        print("....................................")

        print("Step 2 : Round the float coefficients to integer ones")
        print("Generating Lattice...")
        m = self.find_best_m(float_proj, d)
        L = self.gen_lll_reduced_lattice(m)
        print("LLL reduction Done\n")
        print("Babai Rounding...\n")
        rounded_proj = self.babai_reduce(L, float_proj)

        # For DEBUG only: to bypass the previous steps with hard coded result for cosh
        """
        z = self.R.gen()
        P_coeffs = [[1.545137394291834906346252374027227131808897230541152890293501138075153559832365790854033152104244638692052525678175077134781439496696790486*10^-622 ,  5.07*10^-761],[-8.32790318207601907931432227244438963490847512286353681126302230067918812020611200825431107848755484723545462137461803515363519360310606816*10^-608 ,  4.15*10^-746], [2.15635486369886636569889264399565373960991257991217886010897191915387258651581025262084621605653998108363633572965325319415758772736630815*10^-593 ,  6.66*10^-731], [-3.57162770826233082866701867832544508874408828589884689266099893611640966032133296629814291636112353824858007696357057421403750083372442405*10^-579 ,  5.02*10^-717], [4.25170521774436932753684568575643333041273209162668418160502576435429361598307490789349556230772085065412476357596382648590833145302057024*10^-565 ,  4.81*10^-703], [-3.87494658557343844647320428709623292075317427703529706433134779673803538743785356346572050284466775661821261534328293589833854703937218173*10^-551 ,  5.74*10^-689], [2.81261539588233644633219735435050820460072134461518970332190829580286086022676255565002879862937310415844552135043441959724660091317691047*10^-537 ,  6.23*10^-675], [-1.66998785334931918878664172506520747050571182872952243165340349158634602369043852010021867697586474425599820475936199675317116817946173145*10^-523 ,  3.54*10^-661], [8.26766048561478380088967597657042627220060747290016083467317305349819783804985285250040963737907815364488738967931576683304405115654322452*10^-510 ,  3.24*10^-648], [-3.46160612629303387611980528229172583160647542735248373238028585889634414692875597098540264401445737592234744883809422024859782861207152236*10^-496 ,  6.13*10^-634], [1.239019258828938267962615842359633704354250720026597473645793571885785743095307135055595049034197861505590348431159096333753893590953404037*10^-482 ,  2.27*10^-621], [-3.82290248935406453243443095942625821861026379392469456579705892531956509668257163855136453013405519891488822341548934672845049383069625715*10^-469 ,  1.11*10^-607], [1.023364940015230297335316293167225812484612401003607781235336797760354179186116092421201133991122903119699976523713756846524103260979577526*10^-455 ,  8.56*10^-594], [-2.388752897587903478790080750239957998641965850430409324079409241112793531663561817496774615367115749830408446109459916315258075406534346122*10^-442 ,  8.85*10^-581], [4.88091353412388500519648554471518431303069535424231477781720239427445975414439018091556608981814733675862504999290783547013329324703511011*10^-429 ,  1.74*10^-567], [-8.75563326454453946604087105146701907414863745647933363925082260587044683884856024750401573913311906834305301421351361271339013533543389356*10^-416 ,  6.01*10^-554], [1.38187157867630752869292887657084298088202528098789802592924505530974785983754584729647568925174371420290271841866838822358845847052079636*10^-402 ,  6.24*10^-540], [-1.9215033300643080065159395091184451617971747646258444066602644407650784081106690887416801607959271420993007391865603193282945786539797215*10^-389 ,  5.58*10^-526], [2.35631505292619984678659584225434046382574863515912899364550910690672316040518212432298437188413267102960377952321753294978091328698630*10^-376 ,  1.62*10^-511], [-2.5479340578236691140941441493392661524514825501498760961967355762020227092376120276945312845457280871485567159675402170950994638422809*10^-363 ,  8.71*10^-497], [2.43160177918319972331925687675255088360548853548181855320377423889500747818973534554749608539120536421742388744057868133163033549613*10^-350 ,  3.91*10^-482], [-2.0395912264228335579932657595968977747555208392236140251137846266168791338574261930279992114873728406545969216340041258424456399284*10^-337 ,  8.47*10^-468], [1.51815017541539561445851737194997388115162876983428409902014204018362944062982990056748921660377085456005747678257251919980318428*10^-324 ,  4.13*10^-453], [-9.651447641752859109557457538498363893132362795420816622198312195015880970002274628333034075227053966840434307282540377542273388*10^-312 ,  8.99*10^-439], [5.9421688549347453446555152114769973649598327199354370651749631461002100647897914888091876489785242536777658320088485729452308*10^-299 ,  1.76*10^-424], [-2.113451714197533125645958071479087840343754032587592134515687362492156058159893513407470307657887034685672581063026415051004*10^-286 ,  4.57*10^-410], [2.3990922556075894359217731352170474668700752230004647764447497128099679263696218584565757452449377891490304704298630192018*10^-273 ,  2.19*10^-395], [1.468528250568592329232464511867626372526652843646773513512159019331408023669110816003432257973538196799786378777744923061*10^-260 ,  4.65*10^-381], [2.98021231237768923063253248113972453778233324045635138679598982618176500443735768383343350494065855003648405746904801930*10^-247 ,  5.12*10^-367], [3.7343767140819665966563325658774034338392871141187701156518567256599085059960969231276655457624193673233667744810751686*10^-234 ,  7.35*10^-353], [4.547403982211884265922226981755035640819761897136107885094979209936876099342464003180949014492100865829507012260805073*10^-221 ,  4.55*10^-339], [4.88451986349532963770604534541102935065407537639179858452816913125952345457583443271063453116971320929118556900456630*10^-208 ,  5.73*10^-325], [4.6863916693601849028953440302190595662197741953800492830563497361070682720691375808538137549907527853186210029694544*10^-195 ,  7.43*10^-311], [3.973078906761106514219778984144327287850193963861218731147459200643483350382283961635188980307120056987681323308858*10^-182 ,  6.27*10^-297], [2.95438006676864982787829690968197240906520598895554504839979238287125793035124766766720155664679063272916458648901*10^-169 ,  5.70*10^-283], [1.90897856461564741350185748332322591439927646573990570496652135158303639942860153482776437257189573098663473209903*10^-156 ,  5.48*10^-270], [1.0605432242003337707875521764745257422192745478798382720033906289567352777872263700987643579611940581667948141529*10^-143 ,  1.79*10^-256], [5.003584803218891948654584864579559413373225514283428957295243200477818725245060688447200674614995535286046565398*10^-131 ,  4.26*10^-243], [1.97577312345928981219673441794552169308611678483606915292731001829980939482263881165253447352085358698401025109*10^-118 ,  3.48*10^-229], [6.4170363707554819225034166190329752006182174459265761876227810021375854606915359364602015612970871618984008720*10^-106 ,  3.87*10^-216], [1.6783006327028066801229834062233156412600342122501470588352831173089904717847207177961599363723304159117520805*10^-93 ,  1.95*10^-203], [3.442665528828565315999441419798521959286095104129716772569886113049733927038639479451737653622484151624407064*10^-81 ,  3.94*10^-190], [5.355253700442378526237878074670842110927258821964149161817721507905104181308542192758171530148749638873359497*10^-69 ,  7.86*10^-178], [6.041820413088104630752995445573945787937541871462592032599342136152963427627330360941263705104730444956428753*10^-57 ,  8.96*10^-166], [4.647550875769300371181094311944906333100028523784170659549216099091768269887799898114033561390254555411624298*10^-45 ,  8.58*10^-154], [2.224467221428337294354444043315432673353176736722361502441262971042760448341705412047929411021619794386339208*10^-33 ,  7.52*10^-142], [5.703758070814815125615361283919845263371771960933632104045729703111107522506273426188272424042224884033203125*10^-22 ,  7.13*10^-131], [5.850004138877928257384687716357003475959892152546970578441687393933534622192382812500000000000000000000000000*10^-11 ,  6.79*10^-120],[1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 ,  1.68*10^-109]]
        rounded_proj = 0
        for i in range(len(P_coeffs)):
            rounded_proj += self.Fprec(P_coeffs[i][0]).squash().add_error(P_coeffs[i][1]) * z^(96-2*i)
        log_m =  [-99, 0, -133, 0, -170, 0, -208, 0, -247, 0, -286, 0, -326, 0, -367, 0, -408, 0, -449, 0, -491, 0, -532, 0, -574, 0, -617, 0, -659, 0, -702, 0, -745, 0, -788, 0, -831, 0, -875, 0, -919, 0, -962, 0, -1006, 0, -1050, 0, -1094, 0, -1139, 0, -1182, 0, -1228, 0, -1269, 0, -1313, 0, -1356, 0, -1399, 0, -1443, 0, -1487, 0, -1531, 0, -1575, 0, -1619, 0, -1664, 0, -1708, 0, -1753, 0, -1798, 0, -1844, 0, -1889, 0, -1935, 0, -1982, 0, -2028, 0, -2075, 0, -2123, 0, -2172, 0]
        m = [2^log_mi for log_mi in log_m ] 
        """
        coeffs = self.integer_coefficients(rounded_proj, m)
        coeffs, m = self.update_m(coeffs, m)
        z = self.R.gen()
        rounded_proj = self.R(sum([self.Fprec(coeffs[i]).squash() * m[i] * z ^ i for i in range(len(coeffs))]))
        max_size_coefficients = ceil(max(
            [log(abs(RR(coeff)), 2) for coeff in coeffs if RR(coeff) != 0.0]))
        final_security = self.precision_proj(rounded_proj)
        print(" -----------------------")
        print(" | Final security : %s  |" % final_security)
        print(" | Max size coeffs : %s |" % max_size_coefficients)
        print(" -----------------------")
        if max_size_coefficients > self.nbits:
            print("Caution: the size of the coefficients",
                  "exceeded the requirement.")
        print("....................................\n")

        # latex writing
        print("latex exact values poly:\nP &= %i \cdot 2 ^ {%i} \\\\ "
              % (coeffs[0], log(m[0], 2)))
        for i in range(1, len(coeffs)):
                if coeffs[i] != 0:
                    print("  &+ %i \cdot 2 ^ {%i} \cdot x ^ {%i}\\\\ "
                          % (coeffs[i], log(m[i], 2), i))
        print("binary latex exact values poly:\nP &= %s \cdot 2 ^ {%i} \\\\ "
              % (bin(int(coeffs[0])), log(m[0], 2)))
        for i in range(1, len(coeffs)):
                if coeffs[i] != 0:
                    print("  &+ %s \cdot 2 ^ {%i} \cdot x ^ {%i}\\\\ "
                          % (bin(int(coeffs[i])), log(m[i], 2), i))
        return coeffs, [log(m[i], 2) for i in range(1, len(coeffs))]


# -------------------------------- TESTS -------------------------------- #
"""
For approximating the cosh, please uncomment the corresponding parameters.
"""

# EXP on I3

a = -64082
b = 0
f = "exp(x / 92450)"  # 92450 = 2 * 215 ^ 2
invf = "exp(- x / 92450)"
fdiff = "exp(x / 92450)/ 92450"
invfdiff = "- exp(- x / 92450)/ 92450"
nbits = 35
prec = 128
target_precision = 40 #48


# COSH on I2
"""
a = 0  # 5534960
b = 5534960
f = "(exp(x / 92450) + exp(- x / 92450)) / 2"
invf = "2 / (exp(x / 92450) + exp(- x / 92450))"
fdiff = "(exp(x / 92450) - exp(- x / 92450)) / 2 / 92450"
invfdiff = "2 / 92450 * (exp(- x / 92450) - exp( x / 92450)) / \
    (exp(x / 92450) + exp(- x / 92450)) ^ 2"
nbits = 100
prec = 512
target_precision = 65
"""
galactics = GALACTICS_approximation(a, b, prec, target_precision,
                                    nbits, f, fdiff, invf, invfdiff)
galactics.approximate()

# vim: ft=python
