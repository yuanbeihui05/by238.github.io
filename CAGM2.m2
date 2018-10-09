kk=ZZ/32749
ringP3=kk[x_0..x_3]
ringP1=kk[s,t]
cubicMap = map(ringP1,ringP3,{s^3,s^2*t,s*t^2,t^3})
idealCubic= kernel cubicMap
f = vars ringP3
OmegaP3 = kernel f
g = generators OmegaP3
OmegaP3=image g
presentation OmegaP3
OmegaP3res = kernel (f**(ringP3^1/idealCubic))
delta1 = jacobian idealCubic
delta2 = delta1 // (gens OmegaP3res)
delta = map(OmegaP3res, module idealCubic, delta2)

ringP4 = kk[x_0..x_4]
idealX = ideal(x_1+x_3,x_2+x_4)
idealL1 = ideal(x_1,x_2)
idealL2 = ideal(x_3,x_4)
idealY = intersect(idealL1,idealL2)
degree(idealX+idealY)
dim(idealX+idealY)
hilbertSeries (idealY+idealX)

S = QQ[x,y,z];
I = ideal(x^5+y^3+z^3,x^3+y^5+z^3,x^3+y^3+z^5);
I:saturate(I)

PP3 = QQ[t,x,y,z,w]
L = ideal(x,y)
M = ideal(x-t*z,y+t^2*w)
X = intersect(L,M)
Xzero = substitute(saturate(X,t),{t=>0})
degree (ideal(x^2,y)/ideal(x,y^2,z))
degree ideal(x,y^2,z)
