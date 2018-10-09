needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"

mult = (h1,h2) -> (
    if numRows h1 != numRows h2 then error "expected columns vector with same number of rows";
    matrix(for i from 0 to numRows h1 - 1 list {h1_(i,0) * h2_(i,0)})
    )

end--

--Dual Graph is a 4-cycle--
restart
load "./4cycle.m2"

-- input simplicial complex in RR^2
-- Input:4cycles
V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}} -- INPUT
V = {{0,0},{1,0},{0,1},{-1,0},{0,-1}} -- |W|=4
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}} -- INPUT

--Input: 6cycles
V={{0,0},{1,0},{1,1},{-1,1},{-1,0},{-1,-1},{1,-1}} --|W|=6
F={{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,6},{0,6,1}}

--Input: 8cycles
V={{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}} --|W|=8
F={{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,6},{0,6,7},{0,7,8},{0,8,1}}

R = QQ[x,y,z]
A = QQ[X_0..X_(#V-1)]
L = for f in F list product(f, i -> X_i)
D = simplicialComplex L
Idelta = ideal(D)

r = 24 -- INPUT
M = splineModule(V,F,r,BaseRing=>R)
hilbertSeries(M)
G = gens M
nsplines = numColumns G - 1 -- number of nontrivial splines.

trim ideal for i from 1 to 3 list G_(0,i)
oo:ideal(y^(r+1))
radical oo

S1 = R[s_1 .. s_nsplines, 
    Join=>false, 
    Degrees=> drop(degrees source G, 1)
    ]

ss := matrix{{1_S1}} | vars S1
I1 = ideal flatten  for i from 1 to nsplines list 
  for j from i to nsplines list (
    mij := mult(G_{i},G_{j}) // G;
    s_i * s_j - (ss * substitute(mij,S1))_(0,0)
    )
assert isHomogeneous I1
netList decompose I1

--Barycentric Coordinate map--

-- hxyz is the x, y, z functions in terms of Courant functions
hxyz = {sum for i from 0 to #V-1 list V_i_0 * X_i,
    sum for i from 0 to #V-1 list V_i_1 * X_i,
    sum for i from 0 to #V-1 list X_i
    }

-- B consists of the local coordinates for each face
B = for f in F list (
    map(A,R,{sum for i in f list V_i_0 * X_i,
            sum for i in f list V_i_1 * X_i,
            sum for i in f list X_i
            })
    )

-- input is the matrix of splines, output is the some matrix in terms of
-- local coordinates.
barycoord = (mm) -> matrix(
    for i from 0 to numRows mm - 1 list
      for j from 0 to numColumns mm - 1 list B_i(mm_(i,j))
    )
BC = barycoord G

-- generators of S^r in terms of Courant functions.
phi = map(A, S1, join(
        for i from 1 to numColumns BC-1 list (
            sum unique flatten((flatten entries BC_{i})/terms)
            ),
        hxyz))

J = preimage(phi, Idelta)
assert isHomogeneous J
CJ = decompose J
assert(intersect CJ == J)
assert(J == I1)

-- the minimal components of J and preimages of the components of S^0 
-- should be the same.
compsA = for P in decompose Idelta list preimage(phi,P)
compsB = CJ
compsC = for P in decompose Idelta list P
compsA = compsA/(I -> sort flatten entries gens gb I)
compsB = compsB/(I -> sort flatten entries gens gb I)
assert(set compsA === set compsB)

-- now intersect components 2 by 2:
netList compsA
hashTable for ij in subsets(toList(0..#compsA-1), 2) list ij => trim(compsA_(ij_0) + compsA_(ij_1))
hashTable for ij in subsets(toList(0..#compsC-1), 2) list ij => trim(preimage(phi,compsC_(ij_0) + compsC_(ij_1)))

hashTable for ij in subsets(toList(0..#compsA-1), 3) list ij => trim(compsA_(ij_0) + compsA_(ij_1)+compsA_(ij_2))
hashTable for ij in subsets(toList(0..#compsC-1), 3) list ij => trim(preimage(phi,compsC_(ij_0) + compsC_(ij_1)+compsC_(ij_2)))

---------------------
C1=QQ[HP_0,HP_1,HP_2,SZG1,SZG2,SZG3]
F= map (A,C1,{h_0,h_1,h_2,szg1,szg2,szg3})

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

J:preimage(F,I)
