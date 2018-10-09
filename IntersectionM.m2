loadPackage "LocalRings"
R = ZZ/32003[x,y];
C = ideal"y-x2"; -- parabola
D = ideal"y-1";  -- line
E = ideal"y";    -- line

use R;
P = ideal"y-1,x-1";
RP = localRing(R, P);
assert(length (RP^1/promote(C+D, RP)) == 1)
assert(length (RP^1/promote(C+E, RP)) == 0)

use R;
P = ideal"x,y";  -- origin
RP = localRing(R, P);
assert(length(RP^1/promote(C+D, RP)) == 0)
assert(length(RP^1/promote(C+E, RP)) == 2)

restart
needsPackage "LocalRings"
R = ZZ/32003[x,y];
C = ideal"y-x3";
D = ideal"y-x2";
E = ideal"y";

use R;
P = ideal"x,y";
RP = localRing(R, P);
assert(length(RP^1/promote(C+D, RP)) == 2)
assert(length(RP^1/promote(C+E, RP)) == 3)

use R;
P = ideal"x-1,y-1";
RP = localRing(R, P);
assert(length(RP^1/promote(C+D, RP)) == 1)
assert(length(RP^1/promote(C+E, RP)) == 0)
	 