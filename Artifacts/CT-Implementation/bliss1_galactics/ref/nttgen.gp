q=12289
n=512
g=znprimroot(q)^((q-1)/(2*n))
R=2^32
bitrev(m)=subst(Polrev(binary(n+m)[2..-1]),x,2)
zetas=vector(n,k,lift(g^bitrev(k-1)*R))
zetas_inv=vector(n,k,lift((1/g)^(bitrev(k-1)+1)*R))

/*
q=8380417
n=256
g=znprimroot(q)^((q-1)/(2*n))
R=2^32
bitrev(m)=subst(Polrev(binary(n+m)[2..-1]),x,2)
zetas=vector(n,k,lift(g^bitrev(k-1)*R))
zetas_inv=vector(n,k,lift((1/g)^(bitrev(k-1)+1)*R))
*/
