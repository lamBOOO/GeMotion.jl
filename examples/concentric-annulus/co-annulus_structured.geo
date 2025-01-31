// Command line Parameters
If(!Exists(p))
  p = 5;
EndIf

// Settings
res = 100;

// Parameters
R1 = 2/3;
R2 = 5/3;

Point(1) = {0, 0, 0, res};

Point(1000) = { R1, 0, 0, res};
Point(1001) = {-R1, 0, 0, res};

Point(1100) = { R2, 0, 0, res};
Point(1101) = {-R2, 0, 0, res};

Physical Point("inner") = {1000,1001};
Physical Point("outer") = {1100,1101};
Physical Point("all") = {1000,1001,1100,1101};

Circle(2000) = {1000,1,1001};
Circle(2001) = {1001,1,1000};

Circle(2100) = {1100,1,1101};
Circle(2101) = {1101,1,1100};

Line(2200) = {1000,1100};
Line(2201) = {1001,1101};

Line Loop(3000) = {2000,2001}; Physical Curve("inner",3000) = {2001, 2000};
Line Loop(3100) = {2100,2101}; Physical Curve("outer",3100) = {2101, 2100};
Line Loop(3200) = {2201,2101,-2200,-2001};
Line Loop(3201) = {2201,-2100,-2200,2000};

Physical Curve("all",3200) = {2001, 2000, 2101, 2100};

Plane Surface(4000) = {3200};
Plane Surface(4001) = {3201};
Physical Surface("mesh",4000) = {4000,4001};

Transfinite Curve(2000) = p * 10 + 1 Using Progression 1;
Transfinite Curve(2001) = p * 10 + 1 Using Progression 1;
Transfinite Curve(2100) = p * 10 + 1 Using Progression 1;
Transfinite Curve(2101) = p * 10 + 1 Using Progression 1;
Transfinite Curve(2200) = p * 10 + 1 Using Bump 0.1;
Transfinite Curve(2201) = p * 10 + 1 Using Bump 0.1;
Transfinite Surface(4000);
Transfinite Surface(4001);
Recombine Surface(4000);
Recombine Surface(4001);
