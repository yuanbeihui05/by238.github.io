R=QQ[x,y,z]
m_1=ideal(2*x-3*z,y)
m_2=ideal(3*x-4*z,y)
m_3=ideal(x,7*y-5*z)
m_4=ideal(x,13*y-17*z)
m_0=ideal(x,y)
for r from 1 to 10 do (I_r=intersect(m_1^(r+1),m_2^(r+1),m_3^(r+1),m_4^(r+1)));
trim (I_4+m_0^5)
trim (I_1+m_0^2)
