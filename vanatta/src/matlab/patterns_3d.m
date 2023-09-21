c_sound = 1500;
f = 18.5e3;
lambda = c_sound/f;

c_light = 3e8;
f_equiv = c_light / lambda;

el = dipoleCylindrical;

el.Length = lambda/2;
el.Radius = 1e-3; % 1 mm

ra = linearArray;
ra.Element = el;
ra.NumElements = 4;
ra.ElementSpacing = lambda;
% 
pattern(ra, f_equiv);