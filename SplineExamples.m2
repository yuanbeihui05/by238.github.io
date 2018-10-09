--Examples using the AlgebraicSplines package--
installPackage("AlgebraicSplines")
check "AlgebraicSplines"

--Morgan-Scot
V = {{-1,-1},{1,-1},{0,1},{10,10},{-10,10},{0,-10}};

V'= {{-1,-1},{1,-1},{0,1},{10,10},{-10,10},{1,-10}};

F = {{0,1,2},{2,3,4},{0,4,5},{1,3,5},{1,2,3},{0,2,4},{0,1,5}};
F'= {{2,3,4},{0,4,5},{1,3,5},{1,2,3},{0,2,4},{0,1,5}};
S = QQ[x,y,z]
M = splineModule(V,F,1,BaseRing=>S);

M' = splineModule(V',F,1,BaseRing=>S);

M''=splineModule(V',F',1,BaseRing=>S);
prune HH splineComplex(V,F,1)
prune HH splineComplex(V',F,1)
prune splineModule(V',F',0)
splineDimensionTable(0,5,M)
splineDimensionTable(0,5,M')

--Morgan-Scot like triangulation
V = {{1,0},{0,1},{-1,0},{0,-1},{5,5},{-5,5},{-5,-5},{5,-5}};

--tweak upper right vertex
V' = {{1,0},{0,1},{-1,0},{0,-1},{5,6},{-5,5},{-5,-5},{5,-5}}
F = {{0,1,2},{0,2,3},{0,1,4},{1,2,5},{2,3,6},{0,3,7},{0,4,7},{1,4,5},{2,5,6},{3,6,7}};

prune HH splineComplex(V,F,1)
prune HH splineComplex(V',F,1)

--Dual Graph is a cycle--
needsPackage "AlgebraicSplines"
V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}}
V' ={{0,0},{2,0},{1,1},{-1,2},{0,-3}}
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}}
S = QQ[x,y,z]
M = splineModule(V,F,1,BaseRing=>S)
M' = splineModule(V',F,1,BaseRing=>S)
h1 = super M_{1}
f12 = super M_{2}
g12 = super M_{3}
mult = (h1,h2) -> matrix(for i from 0 to 3 list {h1_(i,0) * h2_(i,0)})
mult(h1,h1)
mult(h1,h1) // (gens M)
mult(h1,f12) // (gens M)

L012=matrix{{(y)^2*(x-y)^2},{0},{0},{0}}
L023=matrix{{0},{(x-y)^2*(2*x+y)^2},{0},{0}}
L034=matrix{{0},{0},{(2*x+y)^2*(3*x+2*y)^2},{0}}
L041=matrix{{0},{0},{0},{(3*x+2*y)^2*(y)^2}}
L023'=matrix{{0},{(x-y)*(2*x+y)},{0},{0}}

L012 // (gens M)
L023 // (gens M)
L034 // (gens M)
L041 // (gens M)

L023' // (gens M)

S2 = splineModule(V,F,2,BaseRing=>S)
S3 = splineModule(V,F,3,BaseRing=>S)
S4 = splineModule(V,F,4,BaseRing=>S)

ModuleMatrix=super M_{0,1,2,3}
e14=ModuleMatrix_(0,3)
e14=9*x*y^2-9*y^3
A=QQ[X_0..X_4]

f14= (x,y)-> 9*x*y^2-9*y^3
f24= (x,y)-> 4*x^3-3*x*y^2-y^3
f13= (x,y)-> 12*x*y^2-3*y^3
f23= (x,y)-> 4*x^2*y+4*x*y^2+y^3
f12= (x,y)-> -7*y^2
f22= (x,y)-> 14*x^2-28*x*y+7*y^2
f32= (x,y)-> 18*x^2-24*x*y+8*y^2


f= (x,y)-> e14
f(x,y)

f14(X_1+X_2,X_2)
f24(X_1+X_2,X_2)
f13(X_1+X_2,X_2)
f23(X_1+X_2,X_2)
f12(X_1+X_2,X_2)
f22(X_1+X_2,X_2)
f32(X_1+X_2,X_2)

f14(X_2-X_3,X_2+2*X_3)
f24(X_2-X_3,X_2+2*X_3)
f13(X_2-X_3,X_2+2*X_3)
f23(X_2-X_3,X_2+2*X_3)
f12(X_2-X_3,X_2+2*X_3)
f22(X_2-X_3,X_2+2*X_3)
f32(X_2-X_3,X_2+2*X_3)

f32(-X_3-2*X_4,2*X_3-3*X_4)

f(X_1+X_2,X_2)

--Schlegel diagram of cube--
V = {{1,1},{-1,1},{-1,-1},{1,-1},{2,2},{-2,2},{-2,-2},{2,-2}}
--move two top vertices up simultaneously
V' = {{1,1},{-1,1},{-1,-1},{1,-1},{2,3},{-2,3},{-2,-2},{2,-2}}
--just move one of top vertices up
V''= {{1,1},{-1,1},{-1,-1},{1,-1},{2,3},{-2,2},{-2,-2},{2,-2}}
F = {{0,1,2,3},{0,1,4,5},{0,3,4,7},{2,3,6,7},{1,2,5,6}}
prune splineModule(V,F,0)
prune splineModule(V',F,0)
prune splineModule(V'',F,0)
prune HH splineComplex(V,F,0)
prune HH splineComplex(V',F,0)
prune HH splineComplex(V'',F,0)

--Hal and Dalbec
--free for r<=4, not free for r=5, free for r=6
V={{0,0,0},{3,0,0},{0,3,0},{0,0,3},{-2,-2,1},{1,0,-3}};

F={{0,1,2,3},{0,1,2,5},{0,2,3,4},{0,2,4,5},{0,1,3,4},{0,1,4,5}};

scan(10,r->print{r,pdim splineModule(V,F,r)})

--Centrally triangulated octahedron--
V={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};

F={{0,1,2,3},{0,1,2,6},{0,2,3,4},{0,2,4,6},{0,1,3,5},{0,3,4,5},{0,4,5,6},{0,1,5,6}};

prune splineModule(V,F,1)
--not coned
prune splineModule(V,F,1,Homogenize=>false)
--user created ring
S=QQ[x,y,z]
M=splineModule(V,F,1,Homogenize=>false,BaseRing=>S)
--get dimension of spline spaces from degrees 0 to 5
splineDimensionTable(0,5,M)
--compare hilbert function and polynomial
hilbertComparisonTable(0,5,M)
postulationNumber M
--tweak top vertex
V'={{0,0,0},{1,0,0},{0,1,0},{1,1,1},{-1,0,0},{0,-1,0},{0,0,-1}};

M'=prune splineModule(V',F,1,Homogenize=>false,BaseRing=>S)
splineDimensionTable(0,5,M')
hilbertComparisonTable(0,5,M')
postulationNumber M'
--spline Complexes--
C=cellularComplex(F,InputType=>"Simplicial",Homogenize=>false,BaseRing=>S)
SC=splineComplex(V,F,1,Homogenize=>false,BaseRing=>S)
IC=idealsComplex(V,F,1,Homogenize=>false,BaseRing=>S)
SC==(coker inducedMap(C,IC))
prune HH C
prune HH SC
SC'=splineComplex(V',F,1,Homogenize=>false,BaseRing=>S)
prune HH SC'

--centrally sub-divided cube--
V={{0,0,0},{1,1,1},{-1,1,1},{-1,1,-1},{1,1,-1},{1,-1,1},{-1,-1,1},{-1,-1,-1},{1,-1,-1}};

F={{0,1,2,3,4},{0,1,4,5,8},{0,5,6,7,8},{0,2,3,6,7},{0,1,2,5,6},{0,3,4,7,8}};

S=QQ[x,y,z]
SC=splineComplex(V,F,1,Homogenize=>false,BaseRing=>S)
prune HH SC
associatedPrimes(annihilator HH_2 SC)

--Cube and Octahedron--
V = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {-2, -2, -2}, {-2, 2, -2}, {2, 2, -2}, {2, -2, -2}, {-2, -2, 2}, {-2, 2, 2}, {2, 2, 2}, {2, -2, 2}};

F = {{0, 1, 2, 3, 4, 5}, {0, 8, 9, 12, 13}, {1, 6, 7, 10, 11}, {2, 7, 8, 11, 12}, {3, 6, 9, 10, 13}, {4, 10, 11, 12, 13}, {5, 6, 7, 8, 9}, {0, 2, 8, 12}, {0, 3, 9, 13}, {0, 4, 12, 13}, {0, 5, 8, 9}, {1, 2, 7, 11}, {1, 3, 6, 10}, {1, 4, 10, 11}, {1, 5, 6, 7}, {2, 4, 11, 12}, {3, 4, 10, 13}, {3, 5, 6, 9}, {2, 5, 7, 8}, {0, 2, 4, 12}, {0, 2, 5, 8}, {0, 3, 4, 13}, {0, 3, 5, 9}, {1, 2, 4, 11}, {1, 2, 5, 7}, {1, 3, 4, 10}, {1, 3, 5, 6}};

S = QQ[w,x,y,z]
prune splineModule(V,F,0,BaseRing=>S)
SC = splineComplex(V,F,0,BaseRing=>S);

prune HH_2 SC

--generalizedSplines function--
--derivations on graphic arrangements
S=QQ[(symbol x)_0..(symbol x)_3]
E={{0,1},{1,2},{0,2},{0,3},{2,3}}
k=2;

L=apply(E,e->((x_(first e)-x_(last e))^k))
prune generalizedSplines(E,L)

--integer splines on four cycle
E={{a,b},{b,c},{c,d},{a,d}}
L=apply(E,i->random(3,15))
generalizedSplines(E,L)
--labels taken mod 9
generalizedSplines(E,L,RingType=>9)

--splines on cell determined by arcs
S=QQ[x,y,z]
E={{a,b},{b,c},{a,c}};

r=1;

L={y*z-x^2,x*z+y^2,y*z^2-x^3};

Lr=apply(L,f->f^r);

SM = generalizedSplines(E,Lr)
scan(5,r->(
    Lr = apply(L,f->f^(r+1));
    
    SM = generalizedSplines(E,Lr);
    
    print hilbertPolynomial(SM,Projective=>false);
    
    print ("PostulationNumber="|toString(postulationNumber SM))
    ))
splineDimensionTable(0,30,SM)
hilbertComparisonTable(10,30,SM)
postulationNumber SM
