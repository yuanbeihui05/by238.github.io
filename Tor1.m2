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

-- The chart of QI(r)
LQ= for r from 1 to 30 list {{r, 
betti res QI(r), 
betti res ideal leadTerm QI(r)},
{r,leadTerm gens QI(r) ,ideal leadTerm QI(r)}};
netList ({{"r","betti res Q", "betti res InQ"}}|flatten LQ)


--The chart of V_{r,k}--
r0=20 -- the max of r
k0=20 -- the max of k
DimJ0= for r from 1 to r0 list{{r}| 
    for k from 1 to k0 list hilbertFunction(k-1, comodule QI(r))};
netList flatten {{{"r/k"}| for k from 1 to k0 list k}| flatten DimJ0}
ChartOfV=(r0,k0)->(
    DimJ0= for r from 1 to r0 list{{r}|
	for k from 1 to k0 list hilbertFunction(k-1,comodule QI(r))};
    chart= {{"r/k"}|for k from 1 to k0 list k}|flatten DimJ0;
    return chart)
netList ChartOfV(15,15)









