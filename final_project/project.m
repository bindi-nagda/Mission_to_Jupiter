% Name: Bindi Nagda Space Flight HW 7

% Look-up values 
re = 149.6*10^6;        % Radius of Circular Earth Orbit
rm = 227.9*10^6;       % Radius of Circular Mars Orbit
%rm = 5.79*10^7;         % Radius of Mercury
rj = 778.6*10^6;        % Radius of Circular Jupiter Orbit
rjup = 71490;           % Radius of Jupiter
ue = 398600;
uj = 126686*10^3;
um = 42828;
us = 132712*10^6;
rmars = 3396;          % Radius of Mars
%rmer = 2439.7;
tha = 0;                % True anomaly at Earth Departure
thb = 180;               % True anomaly at Mercury Arrival
r0 = 86029+6378;          % Radius of Earth parking orbit

% First Heliocentric Transfer Orbit
et = (re-rm)/((rm*cosd(thb))-(re*cosd(tha)));  % Eccentricity of first transfer orbit
ht = sqrt(us*re*(1+(et*cosd(tha))));     % Specific angular momentum of first transfer orbit
at = (ht*ht/(us*(1-et^2)));              % Semi-major axis of first transfer orbit
vt1 = us/ht*(1+(et*cosd(thb)));          % Helicentric Tangential Velocity at thb
vr1 = us/ht*(et*sind(thb));              % Heliocentric Radial Velocity at thb

vt0 = us/ht*(1+(et*cos(tha)));           % Helicentric Tangential Velocity at tha
vr0 = us/ht*(et*sin(tha));               % Heliocentric Radial Velocity at tha
vd = sqrt(vt0^2+vr0^2);                  % Heliocentric velocity of spacecraft at Earth departure
ve = sqrt(us/re);                        % Velocity of Earth
dVD = vd - ve;                           % Hyperbolic excess velocity at departure from Earth SOI
C3 = dVD*dVD;
v0 = sqrt(dVD^2 + (2*ue/r0));            % Velocity of s/c on departure hyperbola
v_e = sqrt(ue/(r0));                     % Velocity of s/c on parking orbit around Earth
DV1 = abs(v0 - v_e);                     % Injection burn to escape from Earth SOI
% kit is using v0 as v1 and v_e as v0.
%plot(alt,DV1);

% Mars fly-by
% Patch condition for fly-by
vmars = sqrt(us/rm);                     % Velocity of Mars
v_inf_v = vt1 - vmars;                   % Transversal hyperbolic excess velocity in Mars frame of reference
v_inf_s = -vr1;                          % Radial hyperbolic excess velocity in Mars frame of reference
%v_arrival = sqrt(vt1^2+vr1^2);
v_inf_m = sqrt(v_inf_v^2+v_inf_s^2);     % Hyperbolic excess velocity at Mars SOI
rp = rmars+100;                          % Periareum altitude for fly-by
eh = 1+((rp*v_inf_m^2)/(um));            % Eccentricity of Mars arrival hyperbola

% Determining sign of turn angle, delta, for trailing side fly-by
if (vr1 > 0)
  delta = 2*asin(1/eh);
else 
  delta  = -2*asin(1/eh);
end

phi1 = atan(v_inf_s/v_inf_v);

% Trailing Side Fly-by
phi2 = phi1 + delta;

% Heliocentric tangential and radial velocity components after Mars fly-by
vt2 = vmars+v_inf_m*cos(phi2);
vr2 = -v_inf_m*sin(phi2);
v_inf_m2 = sqrt((vt2-vmars)^2+vr2^2);
v2 = v_inf_m2+vmars;

% Second Heliocentric Transfer Orbit
h2 = rm*vt2;              % Specific angular momentum of second transfer orbit
th2 = atan(((h2/(vr2*rm))-(us/(vr2*h2)))^-1); % True anomaly on second ...
                                              % transfer orbit after fly-by
e2 = vr2*h2/(us*sin(th2));              % Eccentricity of second transfer orbit
a2 = (h2^2)/(us*(1-e2^2));              % Semi-major axis of second transfer orbit
rp2 = a2*(1-e2);                        % Perihel radius
ra2 = a2*(1+e2);                        % Apohel radius

% True anomaly of s/c at Jupiter arrival on second transfer orbit 
th3 = acos((h2^2/(us*rj)-1)/e2);  

% Heliocentric tangential and radial velocity components at Jupiter Arrival
vt3 = us/h2*(1+(e2*cos(th3)));
vr3 = us/h2*(e2*sin(th3));
va = sqrt(vt3^2+vr3^2);
ga = atan(vr3/vt3);

vjup = sqrt(us/rj);                            % Velocity of Jupiter

dVA = sqrt(va^2 + vjup^2-(2*va*vjup*cos(ga))); % Hyperbolic Excess Velocity 
                                               % at Jupiter SOI

vj_inf_v = vt3 - vjup;            % Heliocentric transversal hyperbolic excess velocity
vj_inf_s = -vr3;                  % Heliocentric radial hyperbolic excess velocity
v_inf_j = sqrt(vj_inf_v^2+vj_inf_s^2); % Excess Hyperbolic Velocity at Jupiter SOI
 
v_j = sqrt(uj/(5*rjup));          % Orbital velocity of 5*jupiter radii
v3 = sqrt(dVA^2+(2*uj)/(5*rjup)); % Velocity on Jupiter arrival hyperbola
DV2 = abs(v_j - v3);              % Injection burn to get to target orbit.. 
                                  % ...from second tranfer orbit
                                    
DVT = DV1+DV2;                    % Total delta-V required

% B-plane targetting  
B = (5*rjup)*(sqrt(1+(2*uj)/(5*rjup*v_inf_j^2)));  % B = Aiming Radius
BT = B*cosd(55);           
BR = B*sind(55);
 
% Plot of first transfer orbit from Earth to Mars 
th = 0:0.01:85*pi/180;
R1 = (ht^2/us)./(1+et*cos(th));
x1 = cos(th).*R1;
y1 = sin(th).*R1;
plot(x1,y1,'g','linewidth',3)
ylim ([-10*10^8 10*10^8]);
grid on
hold on
 
% Calculate offset angle between apse lines of the transfer orbits
t = (thb*pi/180)-th2;

% Plot of second transfer orbit to Jupiter after Mars fly-by
theta = (th2+t):0.01:th3+t;                        % Offset angle included
R2 = (h2^2/us)./(1+e2*cos(theta-t));
x2 = cos(theta).*R2;
y2 = sin(theta).*R2;
plot(x2,y2,'r','linewidth',1)
xlabel('Distance (km)');
ylabel('Distance (km)');
title ('Helicentric trajectory from Earth to Jupiter with Mars fly-by')
 
% Plot circular orbits for Earth, Mercury and Jupiter
thm = 0:0.01:2*pi;
xe = re*cos(thm);
ye = re*sin(thm);
plot(xe,ye,'b')

xm = rm*cos(thm);
ym = rm*sin(thm);
plot(xm,ym,'b')

xj = rj*cos(thm);
yj = rj*sin(thm);
plot(xj,yj,'b')

% Draw the planets
% Earth
ex = re*cos(tha);
ey = re*sin(tha);
scatter(ex,ey,'b','linewidth',6);
%text(ex+50000000, ey, 'Earth', 'Fontsize', 12);

% Mars
mx = rm*cos(85*pi/180);
my = rm*sin(85*pi/180); 
scatter(mx,my,'r','linewidth',2);
%text(mx+50000000, my+50000000,'Mars Fly-by', 'Fontsize', 12);

% Jupiter
jx = rj*cos((th3+t));
jy = rj*sin((th3+t)); 
scatter(jx,jy,'o','linewidth',10);
%text(jx+20000000, jy-50000000, 'Jupiter', 'Fontsize', 12);

% Sun
% scatter(0,0,'y','linewidth',29)
