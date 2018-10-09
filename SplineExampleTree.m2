--Dual Graph is a tree--
needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"
V = {{0,0},{-5,-1},{-1,1},{1,1},{0,7},{5,-1},{9,1}}
F = {{0,1,2},{0,2,3},{0,3,5},{3,5,6},{2,3,4}}
S = QQ[x,y,z]
L = {X_0*X_1*X_2,X_0*X_2*X_3,X_0*X_3*X_5,X_3*X_5*X_6,X_2*X_3*X_4}
D = simplicialComplex L
I = ideal(D)
M = splineModule(V,F,1,BaseRing=>S)
--Barycentric Coordinate map--
A=QQ[X_0..X_6]

h_0=X_0+X_1+X_2+X_3+X_4+X_5+X_6
h_1=(-5)*X_1+(-1)*X_2+X_3+5*X_5+9*X_6
h_2=(-1)*X_1+X_2+X_3+7*X_4+(-1)*X_5+X_6

G_0=map(A,S,{(-5)*X_1+(-1)*X_2,(-1)*X_1+X_2,X_0+X_1+X_2})
G_1=map(A,S,{(-1)*X_2+X_3,X_2+X_3,X_0+X_2+X_3})
G_2=map(A,S,{X_3+5*X_5,X_3+(-1)*X_5,X_0+X_3+X_5})
G_3=map(A,S,{X_3+5*X_5+9*X_6,X_3+(-1)*X_5+X_6,X_3+X_5+X_6})
G_4=map(A,S,{(-1)*X_2+X_3,X_2+X_3+7*X_4,X_2+X_3+X_4})

H=map(A,S,{h_1,h_2,h_0})

barycoord = (mm) -> matrix(for i from 0 to 3 list {G_i(mm_(i,0)),G_i(mm_(i,1)),G_i(mm_(i,2)),G_i(mm_(i,3)),G_i(mm_(i,4))})
barycoord(super M_{0,1,2,3,4})


