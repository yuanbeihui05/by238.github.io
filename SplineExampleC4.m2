--Dual Graph is a 4-cycle--
needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"
V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}}
V' ={{0,0},{2,0},{1,1},{-1,2},{0,-3}}
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}}
S = QQ[x,y,z]
A = QQ[X_0..X_4]
L = {X_0*X_1*X_2,X_0*X_2*X_3,X_0*X_3*X_4,X_0*X_4*X_1}
D = simplicialComplex L
I = ideal(D)
M = splineModule(V,F,1,BaseRing=>S)
M' = splineModule(V',F,1,BaseRing=>S)


h1 = super S1_{1}
h2 = super S1_{2}
h3 = super S1_{3}

mult = (h1,h2) -> matrix(for i from 0 to 3 list {h1_(i,0) * h2_(i,0)})

m11=mult(h1,h1) // (gens M)
m12=mult(h1,h2) // (gens M)
m13=mult(h1,h3) // (gens M)
m22=mult(h2,h2) // (gens M)
m23=mult(h2,h3) // (gens M)
m33=mult(h3,h3) // (gens M)

--Courant Coordinate map--

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
barycoord(ModuleMatrix)
szg1 = -7*X_2^2-28*X_2*X_3+98*X_3^2
szg2 = 12*X_1*X_2^2+9*X_2^3+18*X_2^2*X_3
szg3 = 9*X_1*X_2^2-27*X_2^2*X_3

C1=QQ[HP_0,HP_1,HP_2,SZG1,SZG2,SZG3]
F1= map (A,C1,{h_0,h_1,h_2,szg1,szg2,szg3})

MMinBary=matrix{{1, szg1,szg2,szg3}}

baryrel = (mm) -> matrix(for i from 0 to 3 list {H(mm_(i,0))})
baryrel(m11)

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

J1=preimage(F1,I)
primaryDecomposition(J1)

J2=preimage(F1,ideal(X_1,X_4))
gens gb J2
primaryDecomposition(J2)
J3=preimage(F1,ideal(X_1,X_2))
J4=preimage(F1,ideal(X_2,X_3))
J5=preimage(F1,ideal(X_3,X_4))
J1:intersect(J2,J3,J4,J5)
