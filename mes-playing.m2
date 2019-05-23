S = QQ[x,y,z]
I = ideal(x^3, y^3, (x+y)^3)
Z = syz gens I
J = trim ideal Z^{2}
for i from 0 to 10 list basis(i, S^1/J)
kk = frac(QQ[a,b])
S = kk[x,y,z]
I = ideal(x^3, y^3, (a*x+b*y)^3)
syz gens I


for r from 2 to 6 list degree HG1D(V,F,r)
H7 = HG1D(V,F,7);
degree H7
append(o12, 122)

{6,28,80}
  {22, 52}
{16,52,122}
  {52-16, 122-52}
basis(H7)


C = idealsComplex(V,F,2);
prune HH_0 C
