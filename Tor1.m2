
+ M2 --no-readline --print-width 139
Macaulay2, version 1.10
with packages: ConwayPolynomials, Elimination, IntegralClosure, InverseSystems, LLLBases, PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : r=15

o1 = 15

i2 : (m,n)=(3,4)

o2 = (3, 4)

o2 : Sequence

i3 : S=ZZ/32003[x,y,z]

o3 = S

o3 : PolynomialRing

i4 : lin1=for i from 0 to m-1 list (x-i*z)

o4 = {x, x - z, x - 2z}

o4 : List

i5 : lin2=for i from 0 to n-1 list (y-i*z)

o5 = {y, y - z, y - 2z, y - 3z}

o5 : List

i6 : I1=ideal for f in lin1 list f^(r+1)

             16   16      15        14 2       13 3        12 4        11 5        10 6         9 7         8 8         7 9        6 10  
o6 = ideal (x  , x   - 16x  z + 120x  z  - 560x  z  + 1820x  z  - 4368x  z  + 8008x  z  - 11440x z  + 12870x z  - 11440x z  + 8008x z   -
     --------------------------------------------------------------------------------------------------------------------------------------
          5 11        4 12       3 13       2 14        15    16   16      15        14 2        13 3        12 4         11 5       10 6  
     4368x z   + 1820x z   - 560x z   + 120x z   - 16x*z   + z  , x   - 32x  z + 480x  z  - 4480x  z  - 2883x  z  - 11764x  z  + 464x  z  +
     --------------------------------------------------------------------------------------------------------------------------------------
          9 7        8 8       7 9        6 10         5 11        4 12         3 13         2 14           15        16
     7818x z  - 1589x z  - 731x z  + 7424x z   + 15176x z   - 1979x z   - 11091x z   + 13897x z   - 12240x*z   + 1530z  )

o6 : Ideal of S

i7 : I2=ideal for f in lin2 list f^(r+1)

             16   16      15        14 2       13 3        12 4        11 5        10 6         9 7         8 8         7 9        6 10  
o7 = ideal (y  , y   - 16y  z + 120y  z  - 560y  z  + 1820y  z  - 4368y  z  + 8008y  z  - 11440y z  + 12870y z  - 11440y z  + 8008y z   -
     --------------------------------------------------------------------------------------------------------------------------------------
          5 11        4 12       3 13       2 14        15    16   16      15        14 2        13 3        12 4         11 5       10 6  
     4368y z   + 1820y z   - 560y z   + 120y z   - 16y*z   + z  , y   - 32y  z + 480y  z  - 4480y  z  - 2883y  z  - 11764y  z  + 464y  z  +
     --------------------------------------------------------------------------------------------------------------------------------------
          9 7        8 8       7 9        6 10         5 11        4 12         3 13         2 14           15        16   16      15   
     7818y z  - 1589y z  - 731y z  + 7424y z   + 15176y z   - 1979y z   - 11091y z   + 13897y z   - 12240y*z   + 1530z  , y   - 48y  z +
     --------------------------------------------------------------------------------------------------------------------------------------
          14 2         13 3         12 4        11 5         10 6        9 7         8 8       7 9         6 10        5 11        4 12  
     1080y  z  - 15120y  z  - 12595y  z  - 5325y  z  + 13286y  z  + 7066y z  - 15847y z  - 412y z  - 11936y z   - 9562y z   - 4049y z   -
     --------------------------------------------------------------------------------------------------------------------------------------
          3 13         2 14          15        16
     1186y z   + 14478y z   + 7010y*z   + 2686z  )

o7 : Ideal of S

i8 : K1 = ideal (z^(r+1))+I1

             16   16   16      15        14 2       13 3        12 4        11 5        10 6         9 7         8 8         7 9  
o8 = ideal (z  , x  , x   - 16x  z + 120x  z  - 560x  z  + 1820x  z  - 4368x  z  + 8008x  z  - 11440x z  + 12870x z  - 11440x z  +
     --------------------------------------------------------------------------------------------------------------------------------------
          6 10        5 11        4 12       3 13       2 14        15    16   16      15        14 2        13 3        12 4         11 5
     8008x z   - 4368x z   + 1820x z   - 560x z   + 120x z   - 16x*z   + z  , x   - 32x  z + 480x  z  - 4480x  z  - 2883x  z  - 11764x  z 
     --------------------------------------------------------------------------------------------------------------------------------------
           10 6        9 7        8 8       7 9        6 10         5 11        4 12         3 13         2 14           15        16
     + 464x  z  + 7818x z  - 1589x z  - 731x z  + 7424x z   + 15176x z   - 1979x z   - 11091x z   + 13897x z   - 12240x*z   + 1530z  )

o8 : Ideal of S

i9 : K2 = ideal (z^(r+1))+I2

             16   16   16      15        14 2       13 3        12 4        11 5        10 6         9 7         8 8         7 9  
o9 = ideal (z  , y  , y   - 16y  z + 120y  z  - 560y  z  + 1820y  z  - 4368y  z  + 8008y  z  - 11440y z  + 12870y z  - 11440y z  +
     --------------------------------------------------------------------------------------------------------------------------------------
          6 10        5 11        4 12       3 13       2 14        15    16   16      15        14 2        13 3        12 4         11 5
     8008y z   - 4368y z   + 1820y z   - 560y z   + 120y z   - 16y*z   + z  , y   - 32y  z + 480y  z  - 4480y  z  - 2883y  z  - 11764y  z 
     --------------------------------------------------------------------------------------------------------------------------------------
           10 6        9 7        8 8       7 9        6 10         5 11        4 12         3 13         2 14           15        16   16
     + 464y  z  + 7818y z  - 1589y z  - 731y z  + 7424y z   + 15176y z   - 1979y z   - 11091y z   + 13897y z   - 12240y*z   + 1530z  , y  
     --------------------------------------------------------------------------------------------------------------------------------------
          15         14 2         13 3         12 4        11 5         10 6        9 7         8 8       7 9         6 10        5 11  
     - 48y  z + 1080y  z  - 15120y  z  - 12595y  z  - 5325y  z  + 13286y  z  + 7066y z  - 15847y z  - 412y z  - 11936y z   - 9562y z   -
     --------------------------------------------------------------------------------------------------------------------------------------
          4 12        3 13         2 14          15        16
     4049y z   - 1186y z   + 14478y z   + 7010y*z   + 2686z  )

o9 : Ideal of S

i10 : J1=I1:K1_0

              4      3 2         2 3           4        5   5        3 2        2 3          4        5   2 4       5         6
o10 = ideal (x z - 4x z  - 10100x z  - 11795x*z  + 1461z , x  + 6893x z  - 7213x z  - 2592x*z  - 3265z , x z  - 2x*z  - 13576z )

o10 : Ideal of S

i11 : J2=I2:K1_0

                3         4   2 2       4   3         4   4       4
o11 = ideal (y*z  + 16000z , y z  + 840z , y z - 4214z , y  + 227z )

o11 : Ideal of S

i12 : Q=ideal leadTerm (J1+J2)

                3   2 2   3    4   5   4    5   2 4   3 3
o12 = ideal (y*z , y z , y z, y , z , x z, x , x z , x z )

o12 : Ideal of S

i13 : syz gens Q

o13 = {4} | y  0  0  z2 0  0  x2z 0  x3 0  0    0   0   |
      {4} | -z y  0  0  0  0  0   0  0  0  x4   0   0   |
      {4} | 0  -z y  0  0  0  0   0  0  0  0    x4  0   |
      {4} | 0  0  -z 0  0  0  0   0  0  0  0    0   x5  |
      {5} | 0  0  0  -y 0  x2 0   0  0  0  0    0   0   |
      {5} | 0  0  0  0  x  0  0   0  0  z2 -y2z -y3 0   |
      {5} | 0  0  0  0  -z 0  0   0  0  0  0    0   -y4 |
      {6} | 0  0  0  0  0  -z -y  x  0  0  0    0   0   |
      {6} | 0  0  0  0  0  0  0   -z -y -x 0    0   0   |

              9       13
o13 : Matrix S  <--- S

i14 : res Q

       1      9      13      5
o14 = S  <-- S  <-- S   <-- S  <-- 0
                                    
      0      1      2       3      4

o14 : ChainComplex

i15 : 
      res Q

       1      9      13      5
o15 = S  <-- S  <-- S   <-- S  <-- 0
                                    
      0      1      2       3      4

o15 : ChainComplex

i16 : res Q

       1      9      13      5
o16 = S  <-- S  <-- S   <-- S  <-- 0
                                    
      0      1      2       3      4

o16 : ChainComplex

i17 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o17 = ------------------------
                     3
              (1 - T)

o17 : Expression of class Divide

i18 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o18 = ------------------------
                     3
              (1 - T)

o18 : Expression of class Divide

i19 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o19 = ------------------------
                     3
              (1 - T)

o19 : Expression of class Divide

i20 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o20 = ------------------------
                     3
              (1 - T)

o20 : Expression of class Divide

i21 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o21 = ------------------------
                     3
              (1 - T)

o21 : Expression of class Divide

i22 : res Q

       1      9      13      5
o22 = S  <-- S  <-- S   <-- S  <-- 0
                                    
      0      1      2       3      4

o22 : ChainComplex

i23 : hilbertSeries Q

            4     7    9    10
      1 - 4T  + 5T  - T  - T
o23 = ------------------------
                     3
              (1 - T)

o23 : Expression of class Divide

i24 : fact

o24 = fact

o24 : Symbol

i25 : factor(1-4x^4+5x^7-x^9-t^10)
stdio:26:11:(3): warning: character 'x' immediately following number
stdio:26:16:(3): warning: character 'x' immediately following number
stdio:26:10:(3):[1]: error: no method for adjacent objects:
--            4 (of class ZZ)
--             4
--    SPACE   x  (of class S)

i26 : factor(1-4*x^4+5*x^7-x^9-x^10)

             3             2                   2                  2
o26 = (x - 1) (x + 13298)(x  + 11481x - 7879)(x  + 7227x - 1189)(x  + x + 1)(-1)

o26 : Expression of class Product

i27 : (1-4*x^4+5*t^7-x^9-x^10)/(1-x)^3
stdio:28:13:(3): error: no method for binary operator ^ applied to objects:
--            t (of class Symbol)
--      ^     7 (of class ZZ)

i28 : (1-4*x^4+5*x^7-x^9-x^10)/(1-x)^3

       7     6     5      4      3     2
o28 = x  + 4x  + 9x  + 11x  + 10x  + 6x  + 3x + 1

o28 : frac(S)

i29 : 
      
      betti

o29 = betti

o29 : MethodFunctionWithOptions

i30 : betti res Q

             0 1  2 3
o30 = total: 1 9 13 5
          0: 1 .  . .
          1: . .  . .
          2: . .  . .
          3: . 4  3 .
          4: . 3  2 .
          5: . 2  5 2
          6: . .  2 2
          7: . .  1 1

o30 : BettiTally

i31 : betti res Q

             0 1  2 3
o31 = total: 1 9 13 5
          0: 1 .  . .
          1: . .  . .
          2: . .  . .
          3: . 4  3 .
          4: . 3  2 .
          5: . 2  5 2
          6: . .  2 2
          7: . .  1 1

o31 : BettiTally

i32 : hilbertFunction (Q,7)
stdio:35:1:(3): error: no method found for applying hilbertFunction to:
     argument 1 :            3   2 2   3    4   5   4    5   . (of class Ideal)
                   ideal (y*z , y z , y z, y , z , x z, x , x.
     argument 2 :  7 (of class ZZ)

i33 : hilbertFunction(comodule Q,7)
stdio:36:1:(3): error: no method found for applying hilbertFunction to:
     argument 1 :  cokernel | yz3 y2z2 y3z y4 z5 x4z x5 x2z4. (of class Module)
     argument 2 :  7 (of class ZZ)

i34 : hilbertFunction(7,Q)

o34 = 1

i35 : hilbertFunction(6,!)
stdio:38:19:(3): error: syntax error at '!'

i35 : hilbertFunction(6,Q)

o35 = 4

i36 : 