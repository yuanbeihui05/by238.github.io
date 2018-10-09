needsPackage "AlgebraicSplines"

--4-cycle--
V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}} -- INPUT
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}} -- INPUT

for r from 1 to 10 do (J_r=idealsComplex(V,F,r));
netList for r from 1 to 10 list{r,(res J_r_0).dd_1}
