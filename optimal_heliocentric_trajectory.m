%FINAL CODE FOR EARTH-MARS-JUPITER 

clear all;

rEarth=6378; %Radius of earth
rMars=3396;  %Radius of Mars
rJupiter=71490;%Radius of Jupiter
rSEarth=149600000;%Radius of Earth's orbit wrt Sun
rSMars=227900000;%Radius of Mar's orbit wrt Sun 
rSJupiter=778600000;%Radius of Jupiter's orbit wrt Sun
muSun=132712000000;%Standard gravitational parameter of Sun
muMars=42828;%Standard gravitational parameter of Mars
muJupiter=126686000;%Standard gravitational parameter of Jupiter
muEarth=398600;%Standard gravitational parameter of Earth
%thetaM=(85.5598*pi)/180;% Given true anomaly at Mars SOI in radians 

for rp_alt = 200:600:200;    % To vary orbital altitude once we reach Ganymede. Has a big effect on mass.
    j = 1;
    for rpJ = 477100:100:477500;

        for malt = 100:10:350

            thetaM = degtorad(85.45:0.0001:85.5599);

            %TRANSFER ORBIT ELEMENTS 

            eT1= (rSEarth - rSMars)./ ((rSMars.* cos(thetaM))-(rSEarth.* cos(0)));
            %Eccentricity of transfer orbit between Earth SOI to Mars SOI

            hT1=sqrt(muSun.* rSEarth.* ( 1 + (eT1.* cos(0) ) ) );
            %Angular momentum for transfer orbit between Earth SOI and Mars SOI

            aT1 = (hT1.^2)./((muSun).*(1 - (eT1.^2)));

            vP1= (muSun./ hT1).* (1 + eT1.*cos(thetaM) );
            vR1= (muSun./ hT1).* (eT1.* sin(thetaM) );

            vM = sqrt(muSun./ rSMars);% Speed of Mars in the heliocentric plane

            v_d = sqrt(muSun.*((2./rSEarth)-(1./aT1)));
            vE = sqrt(muSun./rSEarth);
            v_inf_dep = v_d - vE;

            %Entereing Mar's Sphere of influence 
            vInfinityV= vP1 - vM;
            vInfinityS= -1.* vR1;

            vInfinity1 = sqrt( (vInfinityV).^2 + (vInfinityS).^2 );

            %Mars flyby

            rp= malt + rMars;% perilume radius

            eF= 1 + ( (rp.* (vInfinity1).^2 )./  muMars );% Flyby orbit's eccentricity

            delta1 = 2 .* asin(1./eF);% Aiming angle

            phi1= atan(vInfinityS./ vInfinityV) + (2.*pi);

            %Trailing side calculations
            % This is the practical choice since we want to have a higher heliocentric
            % speed after the flyby.

            phi2= phi1 + delta1;

            vP2= vM + (vInfinity1.* cos(phi2) );
            vR2= - vInfinity1.* sin(phi2) ;

            hT2= rSMars .* vP2;% Angular momentum for transfer orbit between Mars SOI and Jupiter SOI

            theta2 = atan(( ( hT2./ (vR2.* rSMars)) - (muSun./ (vR2.* hT2))).^-1);

            eT2= (vR2 .* hT2) ./ (muSun .* sin(theta2));% Eccentricity for transfer orbit between Mars SOI and Jupiter SOI

            aT2 = (hT2.^2) ./ ( muSun .* ( 1 - (eT2).^2));% SEMI_MAJOR_AXIS after Mars flyby

            rP2 = aT2.* (1 - eT2);% PERIHEL_RADIUS after Mars flyby

            rA2 = aT2.* (1 + eT2);% APOGEE_RADIUS after Mars flyby

            %B-Plane targeting 

            theta3 =  (acos( (((hT2).^2 / (rSJupiter .* muSun)) - 1)./ eT2));
            % True anomaly for Jupiter and transfer orbit intersection

            vP3 = (muSun ./ hT2) .* (1 + eT2 .* cos(theta3) );
            vR3= (muSun ./ hT2) .* (eT2.*sin(theta3) );

            vJ = sqrt(muSun ./ rSJupiter);% Velocity of Jupiter in the heliocentric plane

            %Entereing Jupiter's Sphere of influence 
            vInfinityV2= vP3 - vJ;% HYPERBOLIC EXCESS VELOCITY TRANSVERSE entering Jupiter's SOI
            vInfinityS2= -1 * vR3;% HYPERBOLIC EXCESS VELOCITY RADIAL entering Jupiter's SOI

            vInfinity2 = sqrt( (vInfinityV2).^2 + (vInfinityS2).^2 );

            %The B-plane 


            b_1 = sqrt((rpJ + (muJupiter./(vInfinity2).^2)).^2 - ((muJupiter./(vInfinity2).^2).^2) );

            inclination = (0.*pi)./180;

            b_T= b_1 .* sin(inclination);
            b_R= b_1 .* cos(inclination);

            %Jupiter Arrival
            eJupiter = 1 + ((rpJ .* (vInfinity2).^2)./ muJupiter);
            a_jupiter = (-muJupiter) ./ (vInfinity2).^2;
            h_Jupiter = ( muJupiter .* a_jupiter .* ( 1 - eJupiter.^2)).^(0.5);

            theta_intersection= acos ( (1./eJupiter) .* ( ((h_Jupiter.^2)./(muJupiter .* 1070412)) - 1 ) );

            %Ganymede orbit - Jupiter incoming hyperbola
            gamma = atan((eJupiter.*sin(theta_intersection))./(1 + (eJupiter.*cos(theta_intersection))));
            v_hyp_ganymedeR = (muJupiter./h_Jupiter).*eJupiter.*sin(theta_intersection);
            v_hyp_ganymedeT = (muJupiter./h_Jupiter).*(1 + (eJupiter.*cos(theta_intersection)));
            v_SpaceCraft_Jupiter = sqrt ( (v_hyp_ganymedeR).^2 + (v_hyp_ganymedeT).^2 );%velocity of spacecraft at Jupiter's perilume


            %Ganymede Arrival
            rG = 1070400;  % Ganymede's orbital radius
            vGanymede = sqrt(muJupiter./ rG);

            v_inf_ganymede = v_SpaceCraft_Jupiter - vGanymede;

            G = 6.67408*10^-20;% Universal gravitational constant   
            mewGanymede = G .* 14819000*10.^16;

            r_p = rp_alt+ 2634.49613;

            e_hyperbolic_into_ganymede = 1 + ( (r_p .* (v_inf_ganymede).^2) ./ mewGanymede);

            v_p_around_ganymede = sqrt ( (v_inf_ganymede).^2 + ((2 .* (mewGanymede))./r_p) );

            v = sqrt(mewGanymede./r_p);

            lastburnV = v_p_around_ganymede - v;

            vstar = 300 * 9.81;

            a = exp(lastburnV.*1000./vstar);

            mp = ((500 .* a) - 500) ./ (1.1 - (0.1 .*a));

            m_0 = 500 + mp + (0.1.*mp);


                for i = 1:1:1100
                    if isreal(theta3(i))
                         mars_ang(malt/10-9,i,j) = thetaM(i);
                         mass(malt/10-9,i,j) = m_0(i);
                    end
                end


    end

        j = j+1;
    end

        mass(mass==0)=NaN;
      
        % Size of X = 1100, Size of Y = 5 (by inspection of mass matrix)
        [X,Y] = meshgrid(85.45:0.0001:85.5599,100:10:180);  
        
        Z = mass(:,:,3);
        C = Z;
        surf(X,Y,Z,C)
        zlim([2700,4900]);
        xlim([85.45, 85.56]);
        zlabel ('Mass (kg)')
        xlabel('ThetaM (deg)')
        ylabel('Mars Fly-by Altitude (km)')
        title(['Variation of mass with ThetaM and Altitude, for different Ganymede orbital ';'altitude values at fixed RpJ = 477300 km                                   ']);
        colormap hsv
        colorbar
        hold on
        scatter3(85.5598,100,mass(1,1099,3))  
end   

figure;

% Plot of heliocentric trajectory using values that yield best mass (by
% inspection of surface plot)

% BEST POSSIBLE VALUES:
% fly by altitude: malt = 110 km
% rp_alt = 200 km
% thetaM = 85.5517 deg
% Corresponds to element 1018 for all other 1x1018 arrays at malt = 110

% Plot of first transfer orbit from Earth to Mars 
th = 0:0.01:85.3*pi/180;
R1 = (hT1(1018)^2/muSun)./(1+eT1(1018)*cos(th));
x1 = cos(th).*R1;
y1 = sin(th).*R1;
plot(x1,y1,'g','linewidth',3)
ylim ([-10*10^8 10*10^8]);
grid on
hold on

% The following are the 1018th element of their respective arrays that yield
% lowest mass (for Mars fly-by altitude = 100 and ThetaM= 85.5177))
thetaM = 85.551699999999997;
theta2 = 1.275897618738569;
theta3 = 3.139098908060187;
eT2 = 0.651843179630087;
hT2 = 5.997924774895925e+09;

% Calculate offset angle between apse lines of the transfer orbits
t = (thetaM*pi/180)-theta2;

% Plot of second transfer orbit to Jupiter after Mars fly-by
theta = (theta2+t):0.01:(theta3+t);                        % Offset angle included
R2 = (hT2^2/muSun)./(1+eT2*cos(theta-t));
x2 = cos(theta).*R2;
y2 = sin(theta).*R2;
plot(x2,y2,'r','linewidth',1)
xlabel('Distance (km)');
ylabel('Distance (km)');
title ('Helicentric trajectory from Earth to Ganymede with Mars fly-by')
 
% Plot circular orbits for Earth, Mars and Jupiter
thm = 0:0.01:2*pi;
xe = rSEarth*cos(thm);
ye = rSEarth*sin(thm);
plot(xe,ye,'b')

xm = rSMars*cos(thm);
ym = rSMars*sin(thm);
plot(xm,ym,'b')

xj = rSJupiter*cos(thm);
yj = rSJupiter*sin(thm);
plot(xj,yj,'b')

% Draw the planets
% Earth
ex = rSEarth*cos(0);
ey = rSEarth*sin(0);
scatter(ex,ey,'b','filled','linewidth',4);
text(ex+50000000, ey, 'Earth', 'Fontsize', 12);

% Mars
mx = rSMars*cos(85*pi/180);
my = rSMars*sin(85*pi/180); 
scatter(mx,my,'r','filled','linewidth',2);
text(mx+50000000, my+50000000,'Mars Fly-by', 'Fontsize', 12);

% Jupiter
jx = rSJupiter*cos((theta3+t));
jy = rSJupiter*sin((theta3+t)); 
scatter(jx,jy,'o','filled','linewidth',6);
text(jx+20000000, jy-50000000, 'Jupiter', 'Fontsize', 12);

%Sun
scatter(0,0,'y','linewidth',29)
        
        
 