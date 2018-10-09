a=5, b=3, r=3
Bbbk=QQ
R=Bbbk[x,y,z,A,B, Degrees=>{1,1,1,a,b}]
polyn1=product flatten(entries (random(Bbbk^a,Bbbk^2)*matrix{{x},{z}}))
polyn2=product flatten(entries (random(Bbbk^b,Bbbk^2)*matrix{{y},{z}}))

J=ideal(A-polyn1,B-polyn2)
isHomogeneous J
I=(ideal(A,y))^r+(ideal(B,x))^r
isHomogeneous (I+J)
netList{{betti res J, betti res (ideal(A,y))^r,betti res ((ideal(A,y))^r+J), betti res I, betti res(I+J),betti res tor1}}
tor1=Tor_1(cokernel (gens I),cokernel (gens J))
betti res oo

f1=1-2*x^3-2*x^4-x^6+6*x^7+2*x^8-2*x^9+4*x^10-6*x^11-4*x^12+6*x^13-6*x^14+2*x^15+5*x^16-6*x^17+4*x^18-2*x^20+2*x^21-x^22
f2=(1-x^3-x^6+x^7-x^9+x^10-x^12+x^13)^2*(1-2*x^4+x^8)
factor((f1-f2)/((1-x^4)^2*(1-x)^3))
factor(f2)
factor(f1)
factor(f1-f2)
netList for r from 1 to 10 list{r,intersect((ideal(A,y))^r+(ideal(B,x))^r,(ideal(x,y))^r+J)==intersect((ideal(A,y))^r,(ideal(x,y))^r+J)+intersect((ideal(B,x))^r,(ideal(x,y))^r+J),
}
{trim intersect((ideal(A,y))^r+(ideal(B,x))^r,(ideal(x,y))^r+J),
trim intersect((ideal(A,y))^r,(ideal(x,y))^r+J)+intersect((ideal(B,x))^r,(ideal(x,y))^r+J)}



