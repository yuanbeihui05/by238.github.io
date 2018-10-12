clearAll
R = ZZ/32003[x,y,z]
QuotientIdeal=(a,b,r)->(
    lin1 = for i from 0 to a-2 list (x-random(ZZ/32003)*z);
    lin2 = for i from 0 to b-2 list (y-random(ZZ/32003)*z);
    I1 = ideal for f in lin1 list f^(r+1);
    I2 = ideal for f in lin2 list f^(r+1);
    Q1 = I1: ideal(z^(r+1));
    Q2 = I2: ideal(z^(r+1));
    Q = trim (Q1 + Q2);
    return Q)
QI=r->QuotientIdeal(3,4,r)
LQ= for r from 8 to 18 list {{r, 
betti res QI(r), 
betti res ideal leadTerm QI(r)},
{r, ,ideal leadTerm QI(r)}};
netList ({{"r","betti res Q", "betti res InQ"}}|flatten LQ)













