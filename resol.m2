
S=QQ[x,y,z,MonomialOrder=>{Position=>Down}]
--the 3-line case--
S=ZZ/32003[x,y]
r=3
M=matrix{{x^r,(x-y)^r,y^r}}
C=res cokernel M
betti C
--playing--
S=ZZ/32003[x,y,z]
r=2
M1=matrix{{x,y,z,0,0,0},{0,0,0,x,y,z}}
C1=res cokernel M1
C1.dd


M2=matrix{{x,y,z,0,0},{0,0,x,y,z}}
C2=res cokernel M2
C2.dd
betti C2
hilbertSeries cokernel M2
M2'=matrix{{x,y,z,0,0,0},{0,0,0,x^2,y,z}}
regularity coker M2'

M3=matrix{{x,y,z,0},{0,x,y,z}}
C3=res cokernel M3
C3.dd
betti C3
F3=res minors(2,M3)
F3.dd
betti F3
hilbertSeries cokernel M3
regularity coker M3
M3'=matrix{{x,y,z,0,0,0,0},{0,x,y,z,y^2,y*x,x^2}}
regularity coker M3'


M4=matrix{{x,y,z},{y,z,x}}
C4=res cokernel M4
C4.dd
hilbertSeries cokernel M4

regularity ideal(x^3,y^3,z^3,x^2*y^2-z^4,x*z^4-y^5)

--the 9r/4 conjecture-- 
S=ZZ/32003[x,y,z]
BettiOfKerSyz1=(r)->{
    l10=x;--x+2*y+2*z;
    l20=y;-- -x+2*y+2*z;
    l13=y-3*z-2*x;--x+y+2*z;
    l14=y-3*z; --x-y+2*z;
    l26=x-3*z;-- -x-y+2*z;
    l27=x-2*y-3*z;-- -x+y+2*z;
    l05=x-y;--x;
    M1=matrix{{l13^(r+1),l14^(r+1),l10^(r+1)}};
    M2=matrix{{l20^(r+1),l26^(r+1),l27^(r+1)}};
    M0=matrix{{l10^(r+1),l05^(r+1),l20^(r+1)}};
    Syz1v0=(res image M0).dd_1;
    Syz1v1=(res image M1).dd_1;
    Syz1v2=(res image M2).dd_1;
    Syz1=matrix{{Syz1v1_(2,0),Syz1v1_(2,1),Syz1v0_(0,0),Syz1v0_(0,1),0,0},
    {0,0,Syz1v0_(2,0),Syz1v0_(2,1),Syz1v2_(0,0),Syz1v2_(0,1)}};
    b=betti(res image Syz1);
    return b
}
i=1
netList for j from 0 to 3 list for i from 0 to 3 list BettiOfKerSyz1(4*i+j)
netList for i from 0 to 3 list BettiOfKerSyz1(4*i)
gens gb image Syz1

--the 9r/4 conjecture in chosen coordinates--
Syz0=(r)->{
    w=-x+y+z;
    u=x-z; 
    M0=matrix{{x^(r+1),(z)^(r+1),(x-z)^(r+1)}};
    M1=matrix{{w^(r+1),(w-u)^(r+1),u^(r+1)}};
    M2=matrix{{x^(r+1),(x-y)^(r+1),y^(r+1)}};
    Syz0v0=(res image M0).dd_1;
    Syz0v1=(res image M1).dd_1;
    Syz0v2=(res image M2).dd_1;
    Syz0=matrix{{Syz0v1_(2,0),Syz0v1_(2,1),Syz0v0_(2,0),Syz0v0_(2,1),0,0},
	{0,0,Syz0v0_(0,0),Syz0v0_(0,1),Syz0v2_(0,0),Syz0v2_(0,1)}};
    return Syz0}
-- Hal's coordinates --
Syz0 = (r)->{
    u=x+y;
    v=x-y;
    L={v-3*z-2*u,v-3*z,u,u-v,v,u-3*z,u-2*v-3*z};
    M0=matrix{{L#2^(r+1),L#3^(r+1),L#4^(r+1)}};
    M1=matrix{{L#0^(r+1),L#1^(r+1),L#2^(r+1)}};
    M2=matrix{{L#4^(r+1),L#5^(r+1),L#6^(r+1)}};
    Syz0v0=(res image M0).dd_1;
    Syz0v1=(res image M1).dd_1;
    Syz0v2=(res image M2).dd_1;
    Syz0=matrix{{Syz0v1_(2,0),Syz0v1_(2,1),Syz0v0_(0,0),Syz0v0_(0,1),0,0},
	{0,0,Syz0v0_(2,0),Syz0v0_(2,1),Syz0v2_(0,0),Syz0v2_(0,1)}};
    return Syz0}

Syz0(3) -- it does not generate a skew one and a symmetric one automatically B.Y.

gens (image Syz0(3))

--Hal's code for 4k+4=r+1, to generate a skew one and a symmetric one --
R=QQ[z,x,y,MonomialOrder=>Lex]
NewMM=(k)->{
    A0=matrix apply(2*k+1,i->(apply(2*k+3,j->(((-1)^(2*k+3+i-j))*binomial(4*k+4,2*k+3+i-j)))));
    A1=submatrix(A0,,{0..k,k+2..2*k+2});--Adrop in Hal's notes
    C=transpose gens kernel A1;
    firstHalf=submatrix(C,,{0..k});
    secondHalf=submatrix(C,,{k+1..2*k+1});
    K1=firstHalf|matrix{{0}}|secondHalf;
    Mat1=A0||K1;
    K2=transpose gens kernel Mat1;
    newMat=Mat1||K2;
    return newMat}

NiceMiddle = (k)->(L=NewMM k;
    R2=QQ[a,b];
    Rdeg=super basis(2*k+2,R2);
    gamma=submatrix((L**R2)*transpose Rdeg,{2*k+1,2*k+2},);
    S1=submatrix(syz matrix{{a^(4*k+4),b^(4*k+4),(a-b)^(4*k+4)*gamma_(0,0)}},{0,1},{0..0});
    S2=submatrix(syz matrix{{a^(4*k+4),b^(4*k+4),(a-b)^(4*k+4)*gamma_(1,0)}},{0,1},{0..0});
    m23=map(R,R2,matrix{{x,y}});
    m23(S1|S2))
NiceSyz = (k,L0,L2)->(NL=NewMM k;
    R2=QQ[a,b];
    Rdeg=super basis(2*k+2,R2);
    gamma=submatrix((NL**R2)*transpose Rdeg,{2*k+1,2*k+2},);
    S1=submatrix(syz matrix{{a^(4*k+4),b^(4*k+4),(a-b)^(4*k+4)*gamma_(0,0)}},{0,1},{0..0});
    S2=submatrix(syz matrix{{a^(4*k+4),b^(4*k+4),(a-b)^(4*k+4)*gamma_(1,0)}},{0,1},{0..0});
    m23=map(R,R2,matrix{{L0,L2}});
    m23(S1|S2))

PHM=(i,L)->(D=matrix{{0_R,0_R}};
    E=submatrix(NiceSyz(floor((i-3)/4),L#0,L#2),{1..1},)||D;
    B=NiceSyz(floor((i-3)/4),L#2,L#4);
    C=D||submatrix(NiceSyz(floor((i-3)/4),L#4,L#6),{0..0},);
    M=E|B|C;
    M1=apply(2,l->(apply(6,j->M_(l,j))));
    MM1=map(R^{2:(-i-1)},,matrix M1);
    MM1
    )
(u,v)=(x+y,x-y);
(u,v)=(x,y);
L={v-3*z,v-3*z-2*u,2*u,2*u-2*v,2*v,u-3*z-2*v,u-3*z}
r=7
MHJ=PHM(r,L)
Dsum=(r,L)->(submatrix (PHM(r,L),{0..0},{0..3})++ submatrix (PHM(r,L),{1..1},{2..5}));
FactorMat=(M)->(r1=rank source M;
    r2=rank target M;
    apply(r1,i->(apply(r2,j->(factor M_(j,i))))))
FactorMat(PHM(3,L))
id_(R^2)
r=3

Syz1v1=(r)->(syz submatrix(PHM(r,L),{0..0},{0..3})++id_(R^2))
Syz1v2=(r)->(id_(R^2)++(syz submatrix(PHM(r,L),{1..1},{2..5})))
prune coker (Syz1v1(r)|Syz1v2(r))
IHJ=(r)->(
    S1v1=(syz submatrix (PHM(r,L),{0..0},{0..3}))++id_(R^2);
    S1v2=id_(R^2)++(syz submatrix (PHM(r,L),{1..1},{2..5}));
    coker submatrix(PHM(r,L)*S1v1,{1..1},))

IHJ1=(r)->(submatrix (PHM(r,L),{1..1},)*submatrix(Syz1v1(r),,{0..(rank source Syz1v1(r))-3}))
IHJ2=(r)->(submatrix (PHM(r,L),{1..1},)*submatrix(Syz1v1(r),,{(rank source Syz1v1(r))-2..(rank source Syz1v1(r))-1}))
netList for i from 3 to 3 list for k from 0 to 3 list (betti res IHJ(4*k+i))
netList for i from 3 to 3 list for k from 0 to 3 list (betti res prune coker(Syz1v1(4*k+i)|Syz1v2(4*k+i)))
netList for i from 3 to 3 list for k from 0 to 3 list (reduceHilbert hilbertSeries prune coker (Syz1v1(4*k+i)|Syz1v2(4*k+i)))
netList for i from 3 to 3 list for k from 4 to 4 list (betti res (gin IHJ(4*k+i)))

syz MHJ
syz (Syz1v1|Syz1v2)
Syz2=syz (submatrix(Syz1v1,{2..3},{0..(rank target Syz1v1)-1})|submatrix(Syz1v2,{2..3},{2..(rank target Syz1v2)+1}))
submatrix (Syz1v1,,{0..(rank target Syz1v1)-1})*submatrix(Syz2,{0..(rank target Syz1v1)-1},)
submatrix (Syz1v2,,{2..(rank target Syz1v2)+1})*submatrix(Syz2,{(rank target Syz1v1)..(rank target Syz1v1)+(rank target Syz1v2)-1},)
PHM(3,L)*Syz1v1
IHJ=coker submatrix (PHM(r,L)*Syz1v1,{1..1},)
betti res IHJ
hilbertSeries(IHJ)
netList for k from 0 to 6 list for i from 2*(4*k+3) to 9*k+8 list hilbertFunction(i,IHJ(4*k+3))
netList for k from 0 to 4 list for i from 2*(4*k+3) to 9*k+8 list hilbertFunction(i,cokernel PHM(4*k+3,L))
netList for k from 0 to 4 list for i from (4*k+3) to 9*k+8 list (hilbertFunction(i,cokernel PHM(4*k+3,L))-hilbertFunction(i,cokernel Dsum(4*k+3,L)))
--the 1-total-interiore-edge case--
S=ZZ/32003[x,y,z]
BettiSyz1=(r)->{
    l01=y;
    l02=x+2*y;
    l03=x-2*y;
    l14=x+y-z;
    l15=-x+y+z;
    M0=matrix{{l02^(r+1),l03^(r+1),l01^(r+1)}};
    M1=matrix{{l01^(r+1),l14^(r+1),l15^(r+1)}};
    Syz1v0=(res image M0).dd_1;
    Syz1v1=(res image M1).dd_1;
    Syz1=matrix{{Syz1v0_(2,0),Syz1v0_(2,1),Syz1v1_(0,0),Syz1v1_(0,1)}};
    b=betti(res image Syz1);
    return b}

netList for j from 0 to 3 list for i from 0 to 3 list BettiSyz1(4*i+j)
    


