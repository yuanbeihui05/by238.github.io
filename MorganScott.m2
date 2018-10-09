needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"

mult = (h1,h2) -> (
    if numRows h1 != numRows h2 then error "expected columns vector with same number of rows";
    matrix(for i from 0 to numRows h1 - 1 list {h1_(i,0) * h2_(i,0)})
    )

end--

--Morgan-Scott--
restart
load "4cycle.m2"

-- input simplicial complex in RR^2

V= {{-1,-1},{1,-1},{0,1},{10,10},{-10,10},{1,-10}};-- INPUT
F = {{0,1,2},{2,3,4},{0,4,5},{1,3,5},{1,2,3},{0,2,4},{0,1,5}};-- INPUT Morgan-Scott

F'={{2,3,4},{0,4,5},{1,3,5},{1,2,3},{0,2,4},{0,1,5}};-- the dual graph is a cycle but Delta is not a star
F'={{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,3,5},{1,2,3}};-- the dual graph is a cycle but Delta is not a star


R = QQ[x,y,z]
A = R[X_0..X_(#V-1)]
L = for f in F list product(f, i -> X_i)
D = simplicialComplex L
Idelta = ideal(D)

loadPackage "LexIdeals"
isCM Idelta

r = 1 -- INPUT
--Ic=ideal(x^(r+1),y^(r+1),z^(r+1),(x+y+z)^(r+1),(x+2*y+3*z)^(r+1),(x+4*y+11*z)^(r+1))
--netList for r from 0 to 10 list {r, hilbertSeries ideal(x^(r+1),y^(r+1),z^(r+1),(x+y+z)^(r+1),(x+2*y+3*z)^(r+1),(x+4*y+11*z)^(r+1)),hilbertSeries splineModule(V,F,r,BaseRing=>R)}
--netList for r from 0 to 10 list {r, hilbertSeries ideal(x^(r+1),y^(r+1),z^(r+1),(y+2*z)^(r+1),(x+2*y+3*z)^(r+1),(x+4*y+11*z)^(r+1)),hilbertSeries splineModule(V,F,r,BaseRing=>R)}

M = splineModule(V,F,r,BaseRing=>R)
M'= splineModule(V,F',r,BaseRing=>R)
G = gens M
nsplines = numColumns G - 1 -- number of nontrivial splines.

--netList for r from 0 to 10 list {r,hilbertSeries(splineModule(V,F,r,BaseRing=>R)), hilbertSeries(splineModule(V,F',r,BaseRing=>R))}

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
Ixyz=ideal {x-sum for i from 0 to #V-1 list V_i_0 * X_i,
        y-sum for i from 0 to #V-1 list V_i_1 * X_i,
	    z-sum for i from 0 to #V-1 list X_i
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
phi = map(A, S1, 
        for i from 1 to numColumns BC-1 list (
            sum unique flatten((flatten entries BC_{i})/terms)
            ))

J = preimage(phi, Idelta+Ixyz)
assert isHomogeneous J
isCM J
CJ = decompose J
assert(intersect CJ == J)
assert(J == I1)

-- the minimal components of J and preimages of the components of S^0 
-- should be the same.
compsA = for P in decompose Idelta list preimage(phi,P)
compsB = CJ
compsA = compsA/(I -> sort flatten entries gens gb I)
compsB = compsB/(I -> sort flatten entries gens gb I)
assert(set compsA === set compsB)

-- now intersect components 2 by 2:
netList compsA
hashTable for ij in subsets(toList(0..#compsA-1), 2) list ij => trim(compsA_(ij_0) + compsA_(ij_1))
hashTable for ij in subsets(toList(0..#compsA-1), 3) list ij => trim(compsA_(ij_0) + compsA_(ij_1))