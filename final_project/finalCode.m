%FINAL CODE FOR EARTH-MARS-JUPITER 

clear all;
clc

%Radius of planets and moons
rEarth=6378; %Radius of earth
rMars=3396;  %Radius of Mars
rJupiter=71490; %Radius of Jupiter
rGanymede = 2634.49613; %Radius of Ganymede
rCallisto = 2410; %Radius of Callisto

%Semi-major axis of planets and moons
rSEarth = 149600000; %Radius of Earth's orbit wrt Sun
rSMars = 227900000; %Radius of Mar's orbit wrt Sun 
rSJupiter = 778600000; %Radius of Jupiter's orbit wrt Sun
rSCallisto = 1882709; %Radius of Callisto's orbit wrt Jupiter
rSGanymede = 1070412; %Radius of Ganymede's orbit wrt Jupiter

%Graviational parameters of sun, planets, and moons
G = 6.67408*10^-20; % Universal gravitational constant 
muSun=132712000000; %Standard gravitational parameter of Sun
muMars=42828; %Standard gravitational parameter of Mars
muJupiter=126686000; %Standard gravitational parameter of Jupiter
muEarth=398600; %Standard gravitational parameter of Earth
muCallisto = G *10759000*10^16; %Standard gravitational parameter of Callisto
muGanymede = G *14819000*10^16; %Standard gravitational parameter of Ganymede

thetaM=(85.5517*pi)/180;% Optimized value of true anomaly for Mars fly-by

%Orbital elements of the transfer orbit
eT1= (rSEarth - rSMars) / ((rSMars * cos(thetaM))-(rSEarth * cos(0)));%Eccentricity of transfer orbit between Earth SOI to Mars SOI
hT1=sqrt(muSun * rSEarth * ( 1 + (eT1 * cos(0))));%Angular momentum for transfer orbit between Earth SOI and Mars SOI
aT1 = (hT1^2)/((muSun)*(1 - (eT1^2)));

vP1= (muSun / hT1) * (1 + eT1 * cos(thetaM) );
vR1= (muSun / hT1) * (eT1 * sin(thetaM) );
v_Mars_beforeflyby = (((vP1)^2)+((vR1)^2))^0.5;
vM = sqrt(muSun / rSMars);% Speed of Mars in the heliocentric plane

v_d = sqrt(muSun*((2/rSEarth)-(1/aT1)));
vE = sqrt(muSun/rSEarth);
v_inf_dep = v_d - vE;

%Entereing Mars' Sphere of influence 
vInfinityV= vP1 - vM;
vInfinityS= -1 * vR1;

vInfinity1 = sqrt( (vInfinityV)^2 + (vInfinityS)^2 );

%Mars flyby
rp= 110 + rMars; %Perilume radius around Mars (Optimized)

eF= 1 + ((rp * (vInfinity1)^2 ) / muMars); %Flyby orbit's eccentricity

delta1 = 2 * asin(1/eF); %Turning angle

phi1= (atan(vInfinityS / vInfinityV)) + (2*pi);

%Trailing side calculations
% This is the practical choice since we want to have a higher heliocentric
% speed after the flyby.

phi2= phi1 + delta1;

vP2= vM + (vInfinity1 * cos(phi2) );
vR2= - vInfinity1 * sin(phi2) ;

v_Mars_afterflyby = (((vP2)^2)+((vR2)^2))^0.5;
hT2= rSMars * vP2;% Angular momentum for transfer orbit between Mars SOI and Jupiter SOI

theta2 = atan(( ( hT2 / (vR2 * rSMars)) - (muSun / (vR2 * hT2)))^-1);
eT2= (vR2 * hT2) / (muSun * sin(theta2));% Eccentricity for transfer orbit between Mars SOI and Jupiter SOI
aT2 = (hT2^2) / ( muSun * ( 1 - (eT2)^2));% SEMI_MAJOR_AXIS after Mars flyby 

theta3 =  (acos( (((hT2)^2 / (rSJupiter * muSun)) - 1) / eT2));

% True anomaly for Jupiter and transfer orbit intersection

vP3 = (muSun / hT2) * (1 + eT2 * cos(theta3) );
vR3= (muSun / hT2) * (eT2 * sin(theta3) );

vJ = sqrt(muSun / rSJupiter);% Velocity of Jupiter in the heliocentric plane

%Entereing Jupiter's Sphere of influence 
vInfinityV2= vP3 - vJ;% HYPERBOLIC EXCESS VELOCITY TRANSVERSE entering Jupiter's SOI
vInfinityS2= -1 * vR3;% HYPERBOLIC EXCESS VELOCITY RADIAL entering Jupiter's SOI
vInfinity2 = sqrt( (vInfinityV2)^2 + (vInfinityS2)^2 );

rpJ = 1149700; %Perilume radius around Jupiter (Optimized)

%Jupiter Arrival
eJupiter = 1 + ((rpJ .* (vInfinity2).^2)./ muJupiter);
a_jupiter = (-muJupiter) ./ (vInfinity2).^2;
h_Jupiter = ( muJupiter .* a_jupiter .* ( 1 - eJupiter.^2)).^(0.5);

theta_intersection= acos ( (1./eJupiter) .* ( ((h_Jupiter.^2)./(muJupiter .* rSCallisto)) - 1 ) );%Intersection of incoming hyperbola with Callisto's SOI 
SOI_Jupiter = 48200000; %SOI radius of Jupiter
theta_intersection1= acos ( (1./eJupiter) .* ( ((h_Jupiter.^2)./(muJupiter .* SOI_Jupiter)) - 1 ) );%Intersection of heliocentric ellipse with Jupiter's SOI

%Callisto SOI orbit - Jupiter incoming hyperbola
v_hyp_aR = (muJupiter./h_Jupiter).*eJupiter.*sin(theta_intersection);
v_hyp_aT = (muJupiter./h_Jupiter).*(1 + (eJupiter.*cos(theta_intersection)));
v_SpaceCraft_Jupiter = sqrt ( (v_hyp_aR).^2 + (v_hyp_aT).^2 );
vCallisto = sqrt(muJupiter./rSCallisto);

%Velocities inside Callisto's SOI
v_inf_CallistoV = v_hyp_aT - vCallisto;
v_inf_CallistoS = - v_hyp_aR;
v_inf_Callisto = ((v_inf_CallistoV).^2 + (v_inf_CallistoS).^2 ).^0.5;

rp_C= 660 + rCallisto;% perilume radius around Callsito (Optimized)
eF_a= 1 + ( (rp_C .* (v_inf_Callisto).^2 ) ./  (muCallisto));% Flyby orbit's eccentricity
delta2 = 2 .* asin(1./eF_a);% Aiming angle
phi3= (atan(v_inf_CallistoS ./ v_inf_CallistoV)) + (2*pi);

%Leading side calculations
% This is the practical choice since we want to have a slower juliocentric
% speed after the flyby.

phi4= phi3 - delta2;

%Leaving Callisto's SOI interface
vP_Callisto= vCallisto + (v_inf_Callisto .* cos(phi4));
vR_Callisto= - v_inf_Callisto .* sin(phi4) ;

hT_Callisto= rSCallisto * vP_Callisto;% Angular momentum for transfer orbit between Amlathea SOI and Callisto SOI
theta_Callisto_flyby = (atan(( ( hT_Callisto / (vR_Callisto * rSCallisto)) - (muJupiter / (vR_Callisto * hT_Callisto)))^-1));
eT_Callisto_Ganymede= (vR_Callisto * hT_Callisto) / (muJupiter * sin(theta_Callisto_flyby));% Eccentricity for transfer orbit between Callisto SOI and Ganymede SOI
aT_Callisto_Ganymede = (hT_Callisto^2) / (muJupiter * ( 1 - (eT_Callisto_Ganymede)^2));% SEMI_MAJOR_AXIS after Amalthea flyby
v_SpaceCraft_Callisto =  ( (vP_Callisto).^2 + (vR_Callisto).^2 )^0.5;

% True anomaly for Ganymede and transfer orbit intersection
theta_ganymede =  (acos( (((hT_Callisto)^2 / (rSGanymede * muJupiter)) - 1) / eT_Callisto_Ganymede));

%Velocities entering Ganymede's SOI interface
vP_ganymedeSOI = (muJupiter / hT_Callisto) * (1 + eT_Callisto_Ganymede * cos(theta_ganymede) );
vR_ganymedeSOI = (muJupiter / hT_Callisto) * (eT_Callisto_Ganymede * sin(theta_ganymede) );

v_ganymede = sqrt(muJupiter / rSGanymede);% Velocity of Jupiter in the heliocentric plane

%Entereing Ganymede's Sphere of influence 
vInfinityVganymede= vP_ganymedeSOI - v_ganymede;
vInfinitySganymede= -1 * vR_ganymedeSOI;
v_inf_ganymede = sqrt( (vInfinityVganymede)^2 + (vInfinitySganymede)^2 );


r_p = 200 + rGanymede; %Target radius around Ganymede (Optimized)
e_hyperbolic_into_ganymede = 1 + ( (r_p .* (v_inf_ganymede).^2) ./ muGanymede); %Eccentricity of the hyperbola inside Ganymede SOI
v_p_around_ganymede = sqrt ( (v_inf_ganymede).^2 + ((2 .* (muGanymede))./r_p) ); %Velocity at the perilume of the hyperbola 
v_final_orbit = sqrt(muGanymede./r_p); %Velocity of the 200km orbit around Ganymede
lastburnV = v_p_around_ganymede - v_final_orbit; %Burn to get into the 200km around Ganymede

%Propellant Mass calculations
vstar = 300 * 9.81; %Exit velocity

%Mass calculation using rocket equation
a = exp(lastburnV.*1000./vstar);
spacecraft_drymass = 500; %Given
mp = ((spacecraft_drymass .* a) - spacecraft_drymass) ./ (1.1 - (0.1 .*a)); %Mass of the propellant
m_0 = spacecraft_drymass + mp + (0.1.*mp); %Total spacecraft launch mass
m_f = m_0 - mp; %Mass after burn



