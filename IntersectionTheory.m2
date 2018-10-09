-- Intersection theory in Macaulay2.

-- There is a package, Schubert2 (loosely based on an older Maple
-- package, called Schubert, written by Stein-Arild Stromme and
-- Sheldon Katz).

restart
needsPackage "Schubert2"

-- 1. Smooth projective varieties
-- Schubert2 handles smooth projective varieties and vector bundles
-- on them.  However, instead of equations, a variety is given
-- by its intersection ring (which might only be up to numerical equivalence,
-- or only containing some divisors on the variety, etc.).

P3 = projectiveSpace 3
intersectionRing P3

-- this looks a bit complicated, we will see what this refers to soon.
-- basically, h corresponds to the divisor of a rational section of the bundle S
-- and the H_(2,i) are the chern classes of the quotient bundle Q.
-- This is not a minimal presentation at all:

AP3 = intersectionRing P3
minimalPresentation AP3

-- since our varieties are projective, they come with a map to A(Spec k) = ZZ.
-- this map is called "integral":
integral(h^3) -- negative of the class of a point
integral(H_(2,1)^3) -- class of a point
integral(H_(2,3)) -- class of a point
integral(h) -- integral returns zero if the class is not zero-dimensional.

-- 2. Vector bundles (on a smooth projective variety).  Similarly, a
-- vector bundle (or a coherent sheaf) on a variety X is not given by
-- equations, or as a module.  Instead, it is given solely through its
-- chern classes.

(S,Q) = bundles P3
S
Q
-- one obtains the total chern class of a bundle using "chern"
chern S
chern_1 S
-- this is really 1+h:
1+h == chern S
chern Q
-- The Whitney sum formula says that the following product is 1,
-- since we have an exact sequence 0 --> S --> V --> Q --> 0,
-- where V is a the trivial vector bundle of rank 4, on P^3.
(chern S) * (chern Q)

-- 3. Splitting principle.  This is used to compute chern classes of related bundles.
chern dual Q
chern exteriorPower(2,Q)
chern ((dual S) ** Q)
chern symmetricPower(3,Q)
chern symmetricPower(3, dual S)

-- 4. Example: Let's compute the number of lines on a cubic surface.
-- We first define the Grassmannian:
G = flagBundle{2,2}
 intersectionRing G
 -- later, we can understand this presentation.  For now, we just use
 -- this ring, whatever it is.
 IG = intersectionRing G
 (S,Q) = bundles G
 symmetricPower(3,dual S)
 chern symmetricPower(3,dual S)
 integral chern_4 symmetricPower(3,dual S)
 -- 27 lines!!
 
-- 5. How to create our own abstract varieties, bundles?
X = base(4, Bundle => (E,2,c), Bundle => (F,3,d))
  A = intersectionRing X  
  chern E
  chern F
  chern dual F
  chern(E**F)
  chern(exteriorPower(2,E))
  chern(exteriorPower(2,F))
  chern(exteriorPower(3,F))
  chern(Hom(E,F))
  chern ((dual E) ** F)
  oo == ooo
  segre E
  segre F

-- Riemann-Roch on a curve
X = abstractVariety(1,QQ[K, Degrees => {1}][D,Join=>false])
X.TangentBundle = abstractSheaf(X,Rank=>1,ChernClass=>1-K)
intersectionRing X
OO(D)
chern oo
chern symmetricPower(3, OO(D) ++ OO(D))
chi OO(D)

-- Riemann-Roch on a surface
X = abstractVariety(2,QQ[K,c_2, Degrees => {1..2}][D,Join=>false])
X.TangentBundle = abstractSheaf(X,Rank=>2,ChernClass=>1-K+c_2)
intersectionRing X
OO(D)
chern oo
chern symmetricPower(3, OO(D) ++ OO(D))
chi OO(D)

-- 6. Projective bundles, flag bundles
-- Given a vector bundle E on a variety X, P(E) is an important variety.
-- Let's consider that.

-- The number of space conics intersecting 8 given lines
Gd = flagBundle({1,3}, VariableNames => {,"d"})
  intersectionRing Gd
(Sd,Qd) = bundles Gd
f = projectiveBundle(dual symmetricPower_2 Qd, VariableNames => {,{e}})
integral (2*d1 + e)^8

f1 = projectiveBundle(dual symmetricPower_2 Qd, VariableNames => {,"e"})
flagBundle({5,1}, dual symmetricPower_2 Qd, VariableNames => {,{symbol e}})
flagBundle({4,2}, dual symmetricPower_2 Qd, VariableNames => {,{e1,e2}})
  intersectionRing oo


