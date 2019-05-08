--Spline Examples--
needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"

--no totally interior edge--
V = {{0,0},{1,0},{-1,1},{0,-1}}
F = {{0,1,2},{0,2,3},{0,1,3}}

V = {{0,0},{1,0},{1,1},{-1,2},{-2,-3}}
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4}}

--1 totally interior edge--
V = {{0,0},{133,0},{277,119},{-313,151},{-413,-292},{216,-141}} -- INPUT 1 total interior edge
V = {{0,0},{1,0},{2,1},{-2,1},{-2,-1},{2,-1}} -- Input 1 total interior nongenric
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4},{1,4,5},{1,2,5}} -- INPUT type(4,4)
F = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,2,5}} -- INPUT type(5,3)

V = {{0,0},{17,3},{52,103},{-105,32},{22,-117}} -- INPUT type(3,4)
F = {{0,1,2},{1,2,3},{0,1,3},{0,3,4},{0,2,4}} -- INPUT type(3,4)

--2 adjacent totally interior edges --
V' = {{random 10,random 10},{-100,0},{100,0},{-200,100},{200,100},{-200,-100},{200,-100}}
V = {{0,0},{-1,0},{1,0},{-2,1},{2,1},{-2,-1},{2,-1}}
V'' = {{0,0},{-100+random 10,random 10},{100+random 10,random 10},{-200+random 10,100+random 10},{200+random 10, 100+random 10},{-200+random 10,-100+random 10},{200-random 10,-100-random 10}}
V''={{0,0},{-100,0},{100,3},{-200,100},{200,100},{-200,-100},{200,-100}}
F = {{0,1,3},{0,1,5},{1,3,5},{0,3,4},{0,5,6},{0,2,4},{0,2,6},{2,4,6}};

V = {{0,-1},{-2,0},{2,0},{-4,2},{-4,-2},{0,-2},{4,-2},{4,2},{0,2}}
V'= {{0,0},{-200+random 10,random 10},{200+random 10,random 10},{-400+random 10,200+random 10},{-400+random 10,-200+random 10},{random 10,-200+random 10},{400+random 10,-200+random 10},{400+random 10,200+random 10},{0,200}}
V''= {{0,-1},{-2,0},{3,0},{-4,2},{-4,-2},{0,-2},{6,-2},{6,2},{0,2}}
F = {{0,1,5},{0,2,5},{0,1,8},{0,2,8},{1,3,8},{1,3,4},{1,4,5},{2,5,6},{2,6,7},{2,7,8}}

--Morgan-Scott--
V = {{-1,0},{2,0},{0,2},{-8,11},{1,-10},{10,10}} -- INPUT Morgan-Scott(generic)
V = {{-1,0},{1,0},{0,2},{-10,10},{0,-10},{10,10}} -- INPUT Morgan-Scott
F = {{0,1,2},{0,2,3},{0,1,4},{1,2,5},{0,3,4},{1,4,5},{2,3,5}}-- INPUT Morgan-Scott

--Computing spline modules with the package "AlgebraicSplines"--
S = QQ[x,y,z,MonomialOrder=>{Position=>Down}]
r=5
M = splineModule(V,F,r,BaseRing=>S);
splineDimensionTable(0,5,M)

--functions--
HG1D=(V,F,r)->(prune HH splineComplex(V,F,r,BaseRing=>S))_1;-- H1(R/J)
HG0D=(V,F,r)->(prune HH idealsComplex(V,F,BaseRing=>S))_0; -- H0(J)
HF=(V,F,d,r)-> { HG0=(prune HH idealsComplex(V,F,r))_0; --Hilbert Function of H1(R/J)
    return hilbertFunction(d,HG0)}
MaxMinorPhi=(V,F,r)->{HJ0=(prune HH splineComplex(V,F,r))_1; --The Fitting ideal Ig(phi), where H1(R/J)=coker(phi).
    phi=(res HJ0).dd_1;
    IgPhi=minors(numgens target phi,phi);
    return IgPhi}
    
--print the result--    
-- Betti table of H1(R/J)--
netList for r from 1 to 10 list {r, betti res HG1D(V,F,r), betti res HG1D(V',F,r)}-- betti res HG1D(V''',F',r)}
netList for r from 6 to 10 list {r,betti res HG1D(V,F,r),
betti res HG1D(V',F,r), betti res HG1D(V'',F,r),betti res HG1D(V''',F,r)}
-- Betti table of Ig(phi)--
netList for r from 6 to 10 list{r, betti res MaxMinorPhi(V,F,r),
    betti res MaxMinorPhi(V',F,r), betti res MaxMinorPhi(V'',F,r), betti res MaxMinorPhi(V''',F,r)}
--Hilbert function of H1(R/J)--
LTitle ={"r\\d-r"}|for d from 1 to 10 list d;
netList ({LTitle}| for r from 1 to 10 list ({r}| 
	for d from (r+1) to (r+10) list HF(V,F,d,r)))



-- playing 15 Nov 2018 --
LTitle ={"r\\d-r"}|for d from 1 to 10 list d;
netList ({LTitle}| for r from 6 to 9 list ({r}| 
	for dr from (r+1) to (r+10) list HF(V,F,dr,r)))

netList ({LTitle}| for r from 6 to 9 list ({r}| 
	for dr from (r+1) to (r+10) list HF(V',F,dr,r)))

netList ({LTitle}| for r from 6 to 9 list ({r}| 
	for dr from (r+1) to (r+10) list HF(V'',F,dr,r)))

netList ({LTitle}| for r from 6 to 9 list ({r}| 
	for dr from (r+1) to (r+10) list HF(V''',F,dr,r)))

R = QQ[x,y]
r = 3
I = ideal(x^(r+1), y^(r+1), (x+y)^(r+1))
Z = syz gens I
minors(2,oo)
J = ideal Z^{0}
ideal gens gb J
R = frac(QQ[a,b,c])[x,y]
I = ideal((x)^(r+1), (x-b*y)^(r+1), (x-c*y)^(r+1))
Z = syz gens I
J = ideal Z^{0}
ideal gens gb J
I = ideal(x^(r+1), y^(r+1), (x-y)^(r+1), (x+y)^(r+1))
Z = syz gens I

-- playing Nov 25, 2018 --
ResH1= res HG1D(V,F,2)
betti ResH1
Matr=ResH1.dd_1
gens gb image Matr

---H. Schenck's code ---
HSC =(n)->(apply(n, i->(C=splineComplex(V,F,i);
     hilbertSeries(HH_1 C, Reduce=>true))))
netList for r from 11 to 15 list{r, regularity HH_1 splineComplex(V,F,r)}
DHSC =(n)->(apply(n,i->(C=splineComplex(V,F,i+1);
	    {i+1,degree numerator hilbertSeries(HH_1 C, Reduce=>true)})))
HS =(r)-> hilbertSeries (HH_1 splineComplex(V,F,r), Reduce =>true)
netList for r from 21 to 30 list{r, degree numerator HS(r)}
