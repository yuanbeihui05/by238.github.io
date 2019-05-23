--Hal's code for 4k+4=r+1, to generate a skew one and a symmetric one --
R=QQ[z,x,y,MonomialOrder=>Lex]
--NewMM, input an integer, output a square matrix with last two rows corresponded to coeffients of B--
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