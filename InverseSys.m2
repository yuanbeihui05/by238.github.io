R=QQ[x,y,z,Degrees=>{{1},{1},{1}}]

k=2--INPUT
g=product(flatten entries flatten (random(QQ^2,QQ^2)*matrix{{y},{z}}))
f=product(flatten entries flatten (random(QQ^3,QQ^2)*matrix{{x},{z}}))
I=for k from 1 to 20 list ((ideal(x,g))^k+(ideal(y,f))^k+(ideal(x,y))^k);
Sou= for k from 1 to 20 list {(ideal(x,y*g))^k,(ideal(y,x*f))^k};
Tar=for k from 1 to 20 list {(ideal(x,g))^k,(ideal(y,f))^k,(ideal(x,y))^k};
J= for k from 1 to 20 list {(ideal(x,g))^k+(ideal(x,y))^k,(ideal(y,f))^k+(ideal(x,y))^k};


hilbertFunction(5,J_0_1)+hilbertFunction(5,J_0_0)
hilbertFunction(5,I_5)
gens gb I_4
netList for j from 0 to 9 list{j+1, leadTerm I_j}
netList for j from 0 to 9 list{j+1, numerator hilbertSeries(I_j)}
netList for d from 1 to 20 list {d, 
    for j from 0 to 9 list{hilbertFunction(d,J_j_0)+hilbertFunction(d,J_j_1)+hilbertFunction(d,J_j_2),hilbertFunction(d,I_j)}}
netList for j from 0 to 19 list{j+1,hilbertFunction(2*(j+1),Sou_j_0)+hilbertFunction(2*(j+1),Sou_j_1),hilbertFunction(2*(j+1),Tar_j_0)+hilbertFunction(2*(j+1),Tar_j_1)+hilbertFunction(2*(j+1),Tar_j_2),hilbertFunction(2*(j+1),I_j)}

---U+V+W----
netList for d from 11 to 30 list{d,
    for j from 10 to 14 list {hilbertFunction(d,R)-hilbertFunction(d,Tar_j_2),-hilbertFunction(d,J_j_0)+hilbertFunction(d,Tar_j_2),-hilbertFunction(d,J_j_1)+hilbertFunction(d,Tar_j_2),hilbertFunction(d,R)-hilbertFunction(d,I_j)}}
--M=matrix{{1,0,0},{0,1,0},{0,0,1},{1,1,1},{1,2,3},{1,4,11}}--INPUT
n=numRows M -1
for i from 0 to n do l_i=M_(i,0)*x+M_(i,1)*y+M_(i,2)*z
for i from 0 to n do P_i=ideal(M_(i,2)*x-M_(i,0)*z,M_(i,2)*y-M_(i,1)*z,M_(i,1)*x-M_(i,0)*y)

--Ideal of fat points--
for r from 0 to 6 do J_r=intersect(for i from 0 to n list (P_i^(r+1)))
netList for r from 0 to 6 list {r,hilbertSeries J_r}

--For 2-Adj case --
I1=(k)->(ideal(x,y))^k;
I2=(k)->(ideal(x,z))^k;
J1=(k)->(ideal(x^2-z^2,y))^k;
J2=(k)->(ideal(x,y-z))^k;
J3=(k)->(ideal(z,x^2-y^2))^k;
IAd=(k)->(intersect (I1(k)+J1(k),J2(k))+I2(k)+J3(k));
H0=(d,k)->(hilbertFunction(d,IAd(k))+hilbertFunction(d,J1(k)+I1(k)+J2(k)));
H1=(d,k)->(hilbertFunction(d,J1(k))+hilbertFunction(d,I1(k))+hilbertFunction(d,J2(k))+hilbertFunction(d,I2(k))+hilbertFunction(d,J3(k)));
H2=(d,k)->(hilbertFunction(d,intersect(J1(k),I1(k)))+hilbertFunction(d,intersect(I1(k),I2(k),J2(k)))+hilbertFunction(d,intersect(I2(k),J3(k))));
netList for i from 0 to 7 list for k from 4*i+4 to 5*i+5 list (4*i+3, k, H2(4*i+k+3,k)+H0(4*i+k+3,k)-H1(4*i+k+3,k),H2(4*i+k+3,k),H1(4*i+k+3,k),H0(4*i+k+3,k))
netList for k from 1 to 7 list {k, gens gb IAd(k)}
