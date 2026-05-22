// Command line Parameters
If(!Exists(p))
  p = 4;
EndIf

If(!Exists(eccentricity))
  eccentricity = 0.8;
EndIf

If(!Exists(phi_deg))
  phi_deg = 0;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

// Parameters
R1 = 2/3;
R2 = 5/3;
shift = eccentricity * (R2 - R1);
phi = phi_deg * Pi / 180;
cx = shift * Cos(phi);
cy = shift * Sin(phi);

// Add midpoints of circles
Point(1) = {0, 0, 0, res};
Point(2) = {cx, cy, 0, res};

Point(1000) = {cx, cy + R1, 0, res};
Point(1001) = {cx, cy - R1, 0, res};

Point(1100) = {0,  R2, 0, res};
Point(1101) = {0, -R2, 0, res};

Physical Point("inner") = {1000,1001};
Physical Point("outer") = {1100,1101};
Physical Point("heatfunction_zero_points") = {1000,1001,1100,1101};
Physical Point("all") = {1000,1001,1100,1101};

Circle(2000) = {1000,2,1001};
Circle(2001) = {1001,2,1000};

Circle(2100) = {1100,1,1101};
Circle(2101) = {1101,1,1100};

Line Loop(3000) = {2100,2101};
Physical Curve("outer",3100) = {2101, 2100};
Line Loop(3001) = {2001,2000};
Physical Curve("inner",3000) = {2001, 2000};

Plane Surface(4000) = {3000,3001};
Physical Surface("mesh",4000) = {4000};

Physical Curve("all",3200) = {2001, 2000, 2101, 2100};
