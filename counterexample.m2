S = QQ[x,y,z]
linform = (e, V) -> (
    f := det matrix{{x,y,z},append(V_(e_0), 1), append(V_(e_1), 1)};
    (trim ideal f)_0
    )
containsVertex = method()
containsVertex(ZZ, List) := (v, E) -> positions(E, e -> member(v, e))

V = {{0,-1}, {-2,0}, {2,0}, {0,-2}, {0,2}}
E = {{0, 1}, {0, 2}, {1,4}, {1,3}, {0,4}, {2,4}, {2,3}}
linforms = for e in E list linform(e, V)

H0matrix = method()
H0matrix ZZ := (r) -> (
    e0 := containsVertex(0, E); 
    e1 := containsVertex(1, E);
    e2 := containsVertex(2, E);
    Z0 := syz matrix{(linforms_e0)/(f -> f^(r+1))};
    Z1 := syz matrix{(linforms_e1)/(f -> f^(r+1))};
    Z2 := syz matrix{(linforms_e2)/(f -> f^(r+1))};
    M := (submatrix(Z1, {0}, ) || matrix{{0,0}}) 
         | submatrix(Z0, {0,1}, ) 
         | (matrix{{0,0}} || submatrix(Z2, {0}, ));
    map(S^{2: -r-1},,M)
    )
end--

restart
load "counterexample.m2"
containsVertex(0, E)
containsVertex(2, E)

M = coker H0matrix 3
assert isHomogeneous M
for i from 0 to 20 list hilbertFunction(i, M)
M = H0matrix 11
res matrix basis(25, M);
gens oo
rank oo
det oo
M = coker H0matrix 7
assert isHomogeneous M
for i from 0 to 20 list hilbertFunction(i, M)
