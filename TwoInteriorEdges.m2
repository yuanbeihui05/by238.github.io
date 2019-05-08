newPackage(
        "TwoInteriorEdges",
        Version => "0.1", 
        Date => "",
        Authors => {{Name => "", 
                  Email => "", 
                  HomePage => ""}},
        Headline => "specialized code to investigate Beihui's counterexamples",
        DebuggingMode => true
        )

export {
    "computeCounterexampleModp",
    "containsVertex",
    "degreePart",
    "H0matrix",
    "standardGraph",
    "linearForm"
    }

degreePart = method()
degreePart(ZZ, Matrix) := Matrix => (deg,f) -> (
    lift(matrix basis(deg,f), coefficientRing ring f)
    )

standardGraph = () -> (
    -- return a (V, E) pair
    V := {
        {-2,0}, {0,-1}, {2,0},  -- interior vertices, in line order
        {0,-2}, {0,2}, -- connect to 0th vertex
        {0,2}, -- connect to vertex #1
        {0,2}, {0,-2} -- connect to vertex #2
        };
    E := {{0,1}, {1,2}, -- interior edges
        {0,3}, {0,4}, -- connect to vertex #0
        {1,5}, -- connect to vertex #1
        {2,6}, {2,7} -- connect to vertex #2
        };
    (V,E)
    )

containsVertex = method()
-- given an index 'v' into the list of vertices V, find the edges in E incident to v.
containsVertex(ZZ, List) := (v, E) -> positions(E, e -> member(v, e))

linearForm = method(Options => {Ring => QQ[getSymbol "x", getSymbol "y", getSymbol "z"]})
linearForm(List, List) := opts -> (e, V) -> (
    -- e is a list of 2 indices into V, generally an element of E
    -- returns a linear form in the ring opts.Ring
    R := opts#Ring;
    x := R_0;
    y := R_1;
    z := R_2;
    f := det matrix{{x,y,z},append(V_(e_0), 1), append(V_(e_1), 1)};
    (trim ideal f)_0
    )

H0matrix = method(Options => options linearForm)
H0matrix(ZZ, List, List) := opts -> (r,V,E) -> (
    -- assumptions: 
    --  (1) the first 3 vertices of V are the interior ones, in a line
    --  (2) each of these vertices is connected to precisely 3 others.
    S := opts#Ring;
    e0 := containsVertex(0, E); 
    e1 := containsVertex(1, E);
    e2 := containsVertex(2, E);
    if #e0 =!= 3 or #e1 =!= 3 or #e2 =!= 3 then error "expected 3 edges from the first three vertices";
    linforms := for e in E list linearForm(e, V, opts);
    Z0 := syz matrix{(linforms_e0)/(f -> f^(r+1))};
    Z1 := syz matrix{(linforms_e1)/(f -> f^(r+1))};
    Z2 := syz matrix{(linforms_e2)/(f -> f^(r+1))};
    M := (submatrix(Z0, {0}, ) || matrix{{0,0}}) 
         | submatrix(Z1, {0,1}, ) 
         | (matrix{{0,0}} || submatrix(Z2, {0}, ));
    map(S^{2: -r-1},,M)
    )

computeCounterexampleModp = method()
computeCounterexampleModp ZZ := (r) -> (
  deg1 := floor((9*r+2)/4);
  << "For r = " << r << ", the top degree of HH_0(J) should = " << deg1 << endl;
  (V,E) := standardGraph();
  f := H0matrix(r, V, E);
  Rp := ZZ/32003[gens ring f];
  fp := sub(f, Rp);
  d1 := degreePart(deg1,fp);
  << "  In degree " << deg1 << ": (#rows,#cols,rank)=" << (numrows d1, numcols d1, rank d1) << endl;
  d2 := degreePart(deg1+1,fp);
  << "  In degree " << deg1+1 << ": (#rows,#cols,rank)=" << (numrows d2, numcols d2, rank d2) << endl;
  (fp, d1, d2)
  )        
beginDocumentation()

doc ///
Key
  TwoInteriorEdges
Headline
  private code for investigating counterexamples to Schenck-Stiller
Description
  Text
    Here is an example use.  Note that the (V,E) used here has some real restrictions:
    the first 3 vertices in V are the interior ones, in a line, and the rest are the vertices
    which provide exactly 3 edges from each of these vertices.  If not, @TO "H0matrix"@ will give
    an error.  This example produces a nice version of the presentation matrix of the
    $H_0 J$ module.
  Example
    (V,E) = standardGraph()
    f = H0matrix(1, V, E)
///

doc ///
Key
  H0matrix
  (H0matrix, ZZ, List, List)
Headline
  Create a presentation matrix for H0 J
Usage
  f = H0matrix(r, V, E)
Inputs
  r:ZZ
    the smoothness of the splines, a non-negative integer.
  V:List
    a list of vertices.  The first 3 are the interior ones, in line.
  E:List
    a list of edges.  Each of the 3 interior vertices must have exactly 3 edges incident.
Outputs
  f:Matrix
    the presentation matrix of $HH_0(J)$.
Description
  Text
    Here is the standard counterexample of Beihui.

    The conjecture states that the last degree of $M = H_0(J)$
    is in degree $2r$.  Thus for $r = 2$, we obtain a counterexample
    since this is nonzero in degree 5.
    
    The conjecture is that the last nonzero degree of M is floor((9r+2)/4)
  Example
    (V,E) = standardGraph()
    f = H0matrix(2, V, E)
    M = coker f
    for i from 0 to 10 list hilbertFunction(i, M)
    assert(hilbertFunction(5,M) != 0)
Caveat
  This is not general code
SeeAlso
  "AlgebraicSplines::AlgebraicSplines"
///

TEST ///
  (V,E) = standardGraph()
  assert(V == {{-2, 0}, {0, -1}, {2, 0}, {0, -2}, {0, 2}, {0, 2}, {0, 2}, {0, -2}})
  assert(E == {{0, 1}, {1, 2}, {0, 3}, {0, 4}, {1, 5}, {2, 6}, {2, 7}})
  linforms = for e in E list linearForm(e,V)
  assert(linforms == {x+2*y+2*z, x-2*y-2*z, x+y+2*z, x-y+2*z, x, x+y-2*z, x-y-2*z})
///

TEST ///
  (V,E) = standardGraph()
  f = H0matrix(2, V, E)
  M = coker f
  assert(hilbertFunction(5,M) != 0)
  assert(hilbertFunction(6,M) == 0)
///

TEST ///
  (V,E) = standardGraph()
  f = H0matrix(3, V, E)
  M = coker f
  assert(hilbertFunction(7,M) != 0)
  assert(hilbertFunction(8,M) == 0)
///

TEST ///
  (V,E) = standardGraph()
  M = coker H0matrix(4, V, E)
  assert(hilbertFunction(9,M) != 0)
  assert(hilbertFunction(10,M) == 0)
///

TEST ///
  (V,E) = standardGraph()
  M = coker H0matrix(5, V, E)
  assert(hilbertFunction(11,M) != 0)
  assert(hilbertFunction(12,M) == 0)
///

TEST ///
  (V,E) = standardGraph()
  M = coker H0matrix(6, V, E)
  assert(hilbertFunction(14,M) != 0)
  assert(hilbertFunction(15,M) == 0)
///

end--

restart
needsPackage "TwoInteriorEdges"
(V,E) = standardGraph()
f = H0matrix(2, V, E)
f = H0matrix(7, V, E)

restart
uninstallPackage "TwoInteriorEdges"
restart
check "TwoInteriorEdges"
installPackage "TwoInteriorEdges"
viewHelp TwoInteriorEdges



restart
needsPackage "TwoInteriorEdges"
(V,E) = standardGraph()
f = H0matrix(3, V, E)
degreePart(6,f)
rank degreePart(7,f)
rank degreePart(8,f)
hilbertFunction(6, coker f) == 12-6

-- r=7
  restart
  needsPackage "TwoInteriorEdges"
  (V,E) = standardGraph()
  f = H0matrix(7, V, E)
  degreePart(12,f)
  hilbertFunction(12, coker f) == 30-6
  degreePart(13,f)
  hilbertFunction(13, coker f) == 42-18
  degreePart(16,f) -- 90x90
  rank degreePart(16,f) == 88 -- !!
  d17 = degreePart(17,f);
  assert((numrows d17, numcols d17, rank d17) == (110, 126, 110))

-- r=11
  restart
  needsPackage "TwoInteriorEdges"
  (V,E) = standardGraph()
  f = H0matrix(11, V, E)
  Rp = ZZ/32003[gens ring f]
  fp = sub(f, Rp)
  isHomogeneous fp
  d1 = degreePart(25,fp);
  (numrows d1, numcols d1, rank d1)
  d2 = degreePart(26,fp);
  (numrows d2, numcols d2, rank d2)

-- r=15
  restart
  needsPackage "TwoInteriorEdges"
  (fp, d1, d2) = computeCounterexampleModp 3;
  (fp, d1, d2) = computeCounterexampleModp 7;
  (fp, d1, d2) = computeCounterexampleModp 11;  
  (fp, d1, d2) = computeCounterexampleModp 15;
  (fp, d1, d2) = computeCounterexampleModp 19;  
  (fp, d1, d2) = computeCounterexampleModp 23;
  
