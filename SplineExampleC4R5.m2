--Dual Graph is a 4-cycle--
needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"
r=3
V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}}
V' ={{0,0},{2,0},{1,1},{-1,2},{0,-3}}
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}}
S = QQ[x,y,z]
A = QQ[X_0..X_4]
L = {X_0*X_1*X_2,X_0*X_2*X_3,X_0*X_3*X_4,X_0*X_4*X_1}
D = simplicialComplex L
I = ideal(D)
M = splineModule(V,F,r,BaseRing=>S)
M' = splineModule(V',F,r,BaseRing=>S)

--Barycentric Coordinate map--

h_0=X_0+X_1+X_2+X_3+X_4
h_1=X_1+X_2-X_3-2*X_4
h_2=X_2+2*X_3-3*X_4

G_0=map(A,S,{X_1+X_2,X_2,X_0+X_1+X_2})
G_1=map(A,S,{X_2-X_3,X_2+2*X_3,X_0+X_2+X_3})
G_2=map(A,S,{-X_3-2*X_4,2*X_3-3*X_4,X_0+X_3+X_4})
G_3=map(A,S,{X_1-2*X_4,-3*X_4,X_0+X_1+X_4})
H=map(A,S,{h_1,h_2,h_0})

ModuleMatrix=super M_{0,1,2,3}
barycoord = (mm) -> matrix(for i from 0 to 3 list {G_i(mm_(i,0)),G_i(mm_(i,1)),G_i(mm_(i,2)),G_i(mm_(i,3))})
BM=barycoord(ModuleMatrix)

----------- r=1 ----------------
szg1 = -7*X_2^2-28*X_2*X_3+98*X_3^2
szg2 = 12*X_1*X_2^2+9*X_2^3+18*X_2^2*X_3
szg3 = 9*X_1*X_2^2-27*X_2^2*X_3

----------- r=2 ----------------
szg1 = -108*X_1*X_2^3-81*X_2^4-324*X_2^3*X_3
szg2 = -126*X_1*X_2^3-77*X_2^4-238*X_2^3*X_3+420*X_2^2*X_3^2-196*X_2*X_3^3-2744*X_3^4+4116*X_3^3*X_4
szg3 = 1204*X_1*X_2^3+1155*X_2^4+5628*X_2^3*X_3+6048*X_2^2*X_3^2-9408*X_2*X_3^3-16464*X_3^4+43904*X_3^3*X_4

----------- r=3 ----------------
szg1 = -1050*X_1*X_2^4-749*X_2^5-4340*X_2^4*X_3-4760*X_2^3*X_3^2+15680*X_2^2*X_3^3+6860*X_2*X_3^4-96040*X_3^5+144060*X_3^4*X_4
szg2 = 113190*X_1*X_2^4+79919*X_2^5+459620*X_2^4*X_3+480200*X_2^3*X_3^2-1756160*X_2^2*X_3^3+528220*X_2*X_3^4+7126168*X_3^5-13109460*X_3^4*X_4
szg3 = -1512*X_1^2*X_2^4-1848*X_1*X_2^5-721*X_2^6-3108*X_2^5*X_3-1428*X_2^4*X_3^2-2464*X_2^3*X_3^3-15288*X_2^2*X_3^4+32928*X_2*X_3^5+76832*X_3^6-230496*X_3^5*X_4+172872*X_3^4*X_4^2


C1=QQ[HP_0,HP_1,HP_2,SZG1,SZG2,SZG3]
F= map (A,C1,{h_0,h_1,h_2,szg1,szg2,szg3})
J=preimage(F,I)
gens J
mingens J

gens preimage(F,ideal(X_2,X_3,X_4))
gens preimage(F,ideal(X_2,X_3,X_1))
gens preimage(F,ideal(X_1,X_2,X_4))
gens preimage(F,ideal(X_1,X_3,X_4))

h1 = super M_{1}
h2 = super M_{2}
h3 = super M_{3}

mult = (h1,h2) -> matrix(for i from 0 to 3 list {h1_(i,0) * h2_(i,0)})

m11=mult(h1,h1) // (gens M)
m12=mult(h1,h2) // (gens M)
m13=mult(h1,h3) // (gens M)
m22=mult(h2,h2) // (gens M)
m23=mult(h2,h3) // (gens M)
m33=mult(h3,h3) // (gens M)

MMinBary=matrix{{1, SZG1,SZG2,SZG3}}

baryrel = (mm) -> matrix(for i from 0 to 3 list {G(mm_(i,0))})

G=map(C1,S,{HP_1,HP_2,HP_0})

M11=baryrel(m11)
M12=baryrel(m12)
M13=baryrel(m13)
M22=baryrel(m22)
M23=baryrel(m23)
M33=baryrel(m33)

J11=(MMinBary_(0,1))^2-(MMinBary*M11)_(0,0)
J12=(MMinBary_(0,1))*(MMinBary_(0,2))-(MMinBary*M12)_(0,0)
J13=(MMinBary_(0,1))*(MMinBary_(0,3))-(MMinBary*M13)_(0,0)
J22=(MMinBary_(0,2))^2-(MMinBary*M22)_(0,0)
J23=(MMinBary_(0,2))*(MMinBary_(0,3))-(MMinBary*M23)_(0,0)
J33=(MMinBary_(0,3))^2-(MMinBary*M33)_(0,0)
J=ideal(J11,J12,J13,J22,J23,J33)

preimage(F,I):J
J:preimage(F,I)
