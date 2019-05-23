-- load lines 104 - lines 145 (inclusive, get def of L)
--  of file resol.m2
r = 7
M=PHM(r,L)
LisHomogeneous M
reduceHilbert hilbertSeries coker M
N = compress(M^{0} ++ M^{1})
isHomogeneous N

target M === target N
reduceHilbert hilbertSeries coker N

K = (image N)/(image M)
K1 = prune K
reduceHilbert hilbertSeries K

pts = for p in L list trim minors(2, vars R || diff(vars R, p))
degree intersect (pts/(i -> i^2))
degree intersect (pts/(i -> i^3))
degree intersect (pts/(i -> i^4))
