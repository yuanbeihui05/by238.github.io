needsPackage "AlgebraicSplines"
loadPackage "SimplicialComplexes"

--INPUT--
V = {{0,0},{2,0},{3,1},{-2,1},{-3,-2},{3,-1}} -- INPUT 1 total interior edge
V = {{0,0},{133,0},{277,119},{-313,151},{-413,-292},{216,-141}} -- INPUT 1 total interior edge
V = {{0,0},{17,3},{52,103},{-105,32},{22,-117}} -- INPUT 3+4
V = {{0,0},{0,1},{3,2},{-3,2},{0,-1}} -- INPUT 3+4(not generic)

V = {{-1,0},{2,0},{0,2},{-8,11},{1,-10},{10,10}} -- INPUT Morgan-Scott(generic)
V = {{-1,0},{1,0},{0,2},{-10,10},{0,-10},{10,10}} -- INPUT Morgan-Scott
V = {{-11,0},{23,0},{0,21},{-83,117},{13,-102},{107,103},{6,13}} -- INPUT Morgan-Scott+(generic)
V = {{-1,0},{2,0},{0,-2},{-8,11},{1,-10},{10,10}} -- INPUT Morgan-Scott(3+4+5)


F = {{0,1,2},{0,2,3},{0,3,4},{0,1,4},{1,4,5},{1,2,5}} -- INPUT 1 total interior edge 4+4
F = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,2,5}} -- INPUT 1 total interior edge 5+3
F = {{0,1,2},{1,2,3},{0,1,3},{0,3,4},{0,2,4}} -- INPUT 3+4
F = {{0,1,2},{0,2,3},{0,1,4},{1,2,5},{0,3,4},{1,4,5},{2,3,5}}-- INPUT Morgan-Scott
F = {{0,2,3},{0,1,4},{1,2,5},{0,3,4},{1,4,5},{2,3,5},{0,2,6},{1,2,6},{0,1,6}}-- INPUT Morgan-Scott+
F = {{0,1,2},{0,2,4},{1,2,4},{0,1,3},{0,3,4},{1,4,5},{1,3,5}}-- INPUT 3+4+5
r=2 -- INPUT smoothness

-- code MS+BY
SC = ZZ[x_0..x_(#V-1)]
SComplex = simplicialComplex for f in F list product(f/(i -> SC_i))
CSC = chainComplex SComplex
LMatrixOfOnes = matrix{numColumns CSC.dd_2 : {1_SC}}
bdedges = positions(flatten entries(CSC.dd_2 * LMatrixOfOnes), i -> i % 2 == 1)
iedges = positions(flatten entries(CSC.dd_2 * LMatrixOfOnes), i -> i % 2 == 0)
ivertices = positions(entries submatrix(CSC.dd_1,bdedges), x -> all(x, a -> a == 0))
RMatrixOfOnes=matrix{{#ivertices : 1_SC}}
piedges = positions(flatten entries(RMatrixOfOnes*submatrix(CSC.dd_1,ivertices,)),i -> i % 2 == 1)
tiedges = iedges_(positions(flatten entries(RMatrixOfOnes*submatrix(CSC.dd_1,ivertices,iedges)),i -> i % 2 == 0))

R = QQ[x,y,z]
VM = matrix V |  matrix{numRows matrix V : {1_R}}
linforms = flatten entries(VM * transpose vars R)
verticesOfEdges = (flatten entries faces(1,SComplex))/(m -> (support m)/index)
iedgeIdeals = for e in iedges list trim ideal (for a in verticesOfEdges#e list linforms_a)

SplineInfo = new Type of HashTable
new SplineInfo from hashTable{
        "bdedges" => bdedges,
        "iedges" => iedges,
	"piedges" => piedges,
	"tiedges" => tiedges,
	"ivertices" => ivertices
	    }

-- problem in this code: verticesOfEdges is all edges, we only want interior edges.--no problem if we just keep using all edges. B.Y.
phi = (k) -> (
      edgeideals := for i in iedgeIdeals list i^k;
      vertexideals := for v in ivertices list trim intersect(
      edges := positions(verticesOfEdges, e -> member(v,e));
      ids := iedgeIdeals_edges;
      for i in ids list i^k
	    );
   --return (edgeideals, vertexideals);
   tar := directSum for I in iedgeIdeals list comodule I^k;
   sou := directSum for v from 0 to (#ivertices-1) list comodule vertexideals_v;
   --return(tar,sou);
   mapphi:= map(tar,sou,(i,j)->promote((submatrix(CSC.dd_1,ivertices,iedges))_(j,i),R));
   return(mapphi);
      )
iPIedges = (k)->(
    ids:=iedgeIdeals_piedges;
    result=trim intersect(for i in ids list i^k);
    return(result);
    )

iIedges = (k)->(
    ids:=iedgeIdeals_iedges;
    result=trim intersect(for i in ids list i^k);
    return(result);
    )
     
phi 2
KerOfPhi=for k from 11 to 15 list ker(phi k); 
CokerOfPhi= for k from 11 to 15 list coker(phi k);
TarOfPhi=for k from 11 to 15 list target(phi k);
SouOfPhi=for k from 11 to 15 list source(phi k);


netList for d from 30 to 40 list{
    d,  for j from 0 to 4 list{
	hilbertFunction(d,KerOfPhi_j),hilbertFunction(d,SouOfPhi_j),hilbertFunction(d,TarOfPhi_j),hilbertFunction(d,CokerOfPhi_j)}}

netList for k from 1 to 20 list{
    k,for j from 0 to 4 list {
	hilbertFunction(k,LK_j),hilbertFunction(k,iPIedges(j+1)),hilbertFunction(k,iIedges(j+1))}}

netList for j from 0 to 4 list{j+1,numerator hilbertSeries(KerOfPhi_j)}
netList for j from 0 to 4 list{j+1,numerator hilbertSeries(SouOfPhi_j)}
netList for j from 0 to 4 list{j+1,numerator hilbertSeries(TarOfPhi_j)}
netList for j from 0 to 4 list{j+1,numerator hilbertSeries(CokerOfPhi_j)}
netList for j from 0 to 9 list{j+1,factor(numerator hilbertSeries(CokerOfPhi_j))}
netList for j from 0 to 4 list{j+1,betti (res KerOfPhi_j), betti (res SouOfPhi_j), betti(res TarOfPhi_j), betti(res CokerOfPhi_j)}


