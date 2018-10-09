needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"

--INPUT--
V = {{0,0},{2,0},{3,1},{-2,1},{-3,-2},{3,-1}} -- INPUT 1 total interior edge
V = {{0,0},{133,0},{277,119},{-313,151},{-413,-292},{216,-141}} -- INPUT 1 total interior edge
V = {{-1,0},{2,0},{0,2},{-8,11},{1,-10},{10,10}} -- INPUT Morgan-Scott
F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4},{1,4,5},{1,2,5}} -- INPUT 1 total interior edge 4+4
F = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,2,5}} -- INPUT 1 total interior edge 5+3
F = {{0,1,2},{0,2,3},{0,1,4},{1,2,5},{0,3,4},{1,4,5},{2,3,5}}-- INPUT Morgan-Scott
r=2 -- INPUT smoothness

--Generating the Simplicial Complex with Facets F, and the associated chain complex--
SC = ZZ[x_0..x_(#V-1)]
MF=matrix(F)
SComplex= simplicialComplex for i from 0 to #F-1 list (x_(MF_(i,0))*x_(MF_(i,1))*x_(MF_(i,2)))
CSC=chainComplex SComplex
--Positions of interior vertices in the matrix--
LMatrixOfOnes={}
for i from 0 to numColumns (CSC.dd_2)-1 do (LMatrixOfOnes=LMatrixOfOnes|{{1}})
LabelIntEdges=substitute(CSC.dd_2*matrix(LMatrixOfOnes),ZZ/2)
ColumnOfIntEdges={}-- initial List of the Position of Interior edges in the matrix
ColumnOfBdryEdges={}
for i from 0 to numRows LabelIntEdges-1 do (if LabelIntEdges_(i,0)==0 then ColumnOfIntEdges=(ColumnOfIntEdges|{i}) else ColumnOfBdryEdges=(ColumnOfBdryEdges|{i}))--ColumnVar of the submatrix

--Positions of interior vertices in the matrix--
MBdryEdges=submatrix(CSC.dd_1,,ColumnOfBdryEdges)
RowOfBdryVertices={}
RowOfIntVertices={}
for i from 0 to numRows MBdryEdges-1 do(
    Flag=false,
    for j from 0 to numColumns MBdryEdges-1 do(if MBdryEdges_(i,j)==0 then Flag=Flag else Flag=true),
    if Flag then RowOfBdryVertices=(RowOfBdryVertices|{i}) else RowOfIntVertices=(RowOfIntVertices|{i})
    )

--for i from 0 to #F-1 do (SCF_i=simplicialComplex {x_(MF_(i,0))*x_(MF_(i,1))*x_(MF_(i,2))})
--MDEdges=faces(1,SComplex)
--BdryEdges=terms(substitute(substitute((MDEdges*CSC.dd_2*matrix(LMatrixOfOnes))_(0,0),ZZ/2[x_0..x_(#V-1)]),SC))-- a list of boundary edges of Delta, in terms of x_i
--LDEdges=for i from 0 to (numColumns MDEdges)-1 list(MDEdges_(0,i))-- a list of all edges of Delta, in terms of x_i
--IntEdges=LDEdges - set BdryEdges -- a list of all interior edges of Delta, in terms of x_i
--MBVertices=faces(0,simplicialComplex BdryEdges)
--LBVertices=for i from 0 to (numColumns MBVertices)-1 list(MBVertices_(0,i))--a list of all boundary vertices of Delta, in terms of x_i


--Generating Fat points ideals--
R=QQ[x,y,z]
M=matrix(V)*matrix{{x},{y}}
for i from 0 to #V-1 do(
    for j from 0 to #V-1 do(
	I_{i,j}=ideal(M_(i,0)+z, M_(j,0)+z),
	MP_{i,j}=cokernel(gens(I_{i,j}^(r+1)))
	))--ideal of points in PP^2--
for i from 0 to #RowOfIntVertices-1 do(Q_i=promote(ideal(1),R),
    for j from 0 to numColumns(CSC.dd_1)-1 do(
	if not((CSC.dd_1)_(RowOfIntVertices_i,j)==0) then Q_i=intersect(Q_i,I_(indices((faces(1,SComplex))_(0,j)))^(r+1))
	))--intersection of all fat point ideals around v-- 
--Generating modules
MQ=cokernel(gens(Q_0))
for i from 1 to #RowOfIntVertices-1 do(MV_i=cokernel(gens (Q_i)),MQ=MQ++MV_i
    )
MI=MP_(indices((faces(1,SComplex))_(0,(ColumnOfIntEdges_0))))
for i from 1 to #ColumnOfIntEdges-1 do(
j=ColumnOfIntEdges_i,
MI=MI++MP_(indices((faces(1,SComplex))_(0,j)))
)
phi=map(MI,MQ,(i,j) -> promote((submatrix(CSC.dd_1,RowOfIntVertices,ColumnOfIntEdges))_(j,i),R))
gens ker phi
netList for i from r to 4*r+1 list{i,hilbertFunction(i,ker phi)}

--Some other computations--
for r from 1 to 30 do (N1_r=intersect(I1_2^(r+1),I1_3^(r+1),I1_4^(r+1),I1_5^(r+1)))
for r from 1 to 30 do (N2_r=intersect(I2_2^(r+1),I2_5^(r+1)))
netList for r from 1 to 10 list{r,hilbertFunction(2*r+1,intersect(N1_r+N2_r,I1_1^(r+1))),hilbertFunction(2*r+1,intersect(N1_r,I1_1^(r+1))+intersect(N2_r,I1_1^(r+1)))}
netList for r from 1 to 10 list{r,hilbertFunction(3*r+1,intersect(N1_r+N2_r,I1_1^(r+1))),hilbertFunction(3*r+1,intersect(N1_r,I1_1^(r+1))+intersect(N2_r,I1_1^(r+1)))}
r=8
netList for i from 1 to 10 list{i,hilbertFunction(i,intersect(N1_r+N2_r,I1_1^(r+1))),hilbertFunction(i,intersect(N1_r,I1_1^(r+1))+intersect(N2_r,I1_1^(r+1)))}
netList for r from 16 to 20 list{r,intersect(N1_r+N2_r,I1_1^(r+1))==(intersect(N1_r,I1_1^(r+1))+intersect(N2_r,I1_1^(r+1)))}
intersect(N1_r+N2_r,I1_1^(r+1))==(intersect(N1_r,I1_1^(r+1))+intersect(N2_r,I1_1^(r+1)))

for r from 1 to 10 do (J_r=idealsComplex(V,F,r));
for r from 1 to 10 do (K_r=idealsComplex(V',F,r));
for r from 1 to 10 do (J'_r=idealsComplex(V,F',r));
for r from 1 to 10 do (K'_r=idealsComplex(V',F',r));
for r from 1 to 10 do (G_r=res(HH_0 J_r));
