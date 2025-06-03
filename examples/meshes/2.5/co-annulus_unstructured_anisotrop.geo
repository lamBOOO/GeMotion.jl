// Command line Parameters
If(!Exists(p))
  p = 5;
EndIf

// Settings
res = 100;
Mesh.CharacteristicLengthMax = 1.0 * 2^(-p);
Mesh.MshFileVersion = 2.0;

// Parameters
R1 = 2./3.;  R2 = 5./3.;
Point(1)    = {0,0,0,1};
Point(1000) = { R1, 0, 0, 1};
Point(1001) = {-R1, 0, 0, 1};
Point(1100) = { R2, 0, 0, 1};
Point(1101) = {-R2, 0, 0, 1};

Physical Point("inner") = {1000,1001};
Physical Point("outer") = {1100,1101};
Physical Point("all") = {1000,1001,1100,1101};

Circle(2000) = {1000,1,1001};   // inner
Circle(2001) = {1001,1,1000};
Circle(2100) = {1100,1,1101};   // outer
Circle(2101) = {1101,1,1100};

Line Loop(30) = {2100,2101}; Physical Curve("outer",3100) = {2100,2101};    // outer loop  (positive)
Line Loop(31) = {2001,2000}; Physical Curve("inner",3000) = {2001,2000};   // inner loop  (negative → hole)
Plane Surface(40) = {30,31};    // << ONE surface with a hole
Physical Surface("mesh",4000) = {40};

Physical Curve("all",3200) = {2000, 2001, 2100, 2101};

// Characteristic length away from walls
fac = Sqrt(2)^(-(p+4-1)) ;
clen = 0.2585;
// res = fac * clen ;
// Mesh.CharacteristicLengthMax = res ;

// ------------------------------------------------------------------------
// SIZE FIELD
Field[1] = Distance ;
Field[1].CurvesList = {2000,2001} ;
Field[1].NumPointsPerCurve = 2000 ;

Field[2] = Threshold ;
Field[2].InField  = 1 ;
Field[2].SizeMin    = fac * 1/3 * clen ;   // at inner wall R = R1
Field[2].SizeMax    = fac * 1 * clen ;   // in the far field
Field[2].DistMin  = 0.03  ;
Field[2].DistMax  = 0.03 ;
Field[2].StopAtDistMax  = 0 ;
Field[2].Sigmoid  = 1 ;

Background Field = 2 ;
