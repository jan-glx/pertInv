

P(D|X,Z)=N(D|X+Z,1+Z)
P(X|Z)=N(X|1,0)
P(Z|X)=Bernulli(1/(1+exp(X)))


Z* ~ Bernulli(0.5)
X* ~ P(Z*,D)


Xt,Zt <- X*,Z* with probability alpha else Xtm1,Ztm1

P(Xt,Zt|Xtm1,Ztm1) =  alpha * Bernulli(Zt,0.5) *P(Xt,Zt|D) + (1-alpha) delta(Xt-Xtm1,Zt-Ztm1)
=P(Xtm1,Ztm1|Xt,Zt) * P(Xt,Zt) / P(Xtm1,Ztm1) = (alpha * Bernulli(Ztm1,0.5)*  P(Xtm1,Ztm1|D) + (1-alpha) delta(Xt-Xtm1,Zt-Ztm1)) * P(Xt,Zt) / P(Xtm1,Ztm1)

alpha * Bernulli(Zt,0.5) *P(Xt,Zt|D) = alpha * Bernulli(Ztm1,0.5)*  P(Xtm1,Ztm1|D) * P(Xt,Zt) / P(Xtm1,Ztm1) + (1-alpha) delta(Xt-Xtm1,Zt-Ztm1) * (P(Xt,Zt) / P(Xtm1,Ztm1)-1)

P

Q(X,Z|Xt,Zt) = 1/2 P(X|Z,D)

rho(X*,Z*|Xt,Zt)=min 1, P(X*,Z*|D)/P(Xt,Zt|D)*Q(Xt,Zt)/Q(X*,Z*)

rho(X*,Z*|Xt,Zt)=min 1, P(X*,Z*|D)/P(Xt,Zt|D)*Q(Xt,Zt)/Q(X*,Z*)
