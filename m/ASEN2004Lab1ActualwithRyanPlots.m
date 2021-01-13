% ASEN 2004 Aero Lab part 1
% 01/26/ 2020
% Ankrit Uprety
% Cole Macpherson
% Scott Mansfield
% Ryan Block
% Zach Lesan
%% Housekeeping
clear;
clc;

%% Read in Excel File
dataAirfoil = xlsread("Tempest UAS Airfoil and CFD Data for Lab 1.xlsx");
alphaAirfoil = dataAirfoil(:,1);
ClAirfoil = dataAirfoil(:,2);
CdAirfoil = dataAirfoil(:,3);

dataBody = xlsread("ASEN2004Lab1CFDData.xlsx");
alphaBody = dataBody(:,1);
CLBody = dataBody(:,2);
CDBody = dataBody(:,3);

%% Define Variables
b = 3.22; %m
c = 0.23; %m
S = b*c; %m^2
AR = 16.5;%aspect ratio
tipToTail = 1.56; %m  tip to tail
wEmptyWithPropulsion = 4.5; %kg
GTOW = 6.4; %kg  gross takeoff weight
maxLoadFactor = 14;
wingLoading = 10.2; %kg/m^2
stallSpeed = 12; %m/s
cruiseSpeed = [17 25]; %m/s
endurance = 45*60; %sec
range = 40; %km
climbRate = [2.5 11.5]; %m/s
dFuse = 0.16;%m diamter of fuselage
wingArea = 0.63;

Swet = 1.464;
e = 0.9;
trc = 0.087; %thickness to chord ratio
Cfe = (0.0055+0.003)/2;
s = 1-2*((dFuse/b)^2);
pAlt = 1.0324; %kg/m^3

%% Calculations
% calculate a0
a00 = zeros(20,1);
for i = 1:length(alphaAirfoil)-1
    a00(i) = (ClAirfoil(i+1)-ClAirfoil(i))/(alphaAirfoil(i+1)-alphaAirfoil(i));
end
a00(1:end) = mean(a00(1:12));
a0 = a00(4); 
% calculate a
a = a0/(1+((57.3*a0)/(pi*e*AR)));
% find Cl where there is zero lift
lastNegClAirfoil = 0;
lastNegClAirfoilIndex = 0;
for i = 1:length(ClAirfoil)
    if ClAirfoil(i) < 0
        lastNegClAirfoil = ClAirfoil(i);
        lastNegClAirfoilIndex = i;
    end
end


slopeAtLastNegClAirfoil = a; 
% calculate alpha at zero lift
alphaZeroLift = (slopeAtLastNegClAirfoil*alphaAirfoil(lastNegClAirfoilIndex)-lastNegClAirfoil)/slopeAtLastNegClAirfoil;
% find the minimum cd of the airfoil
CdAirfoilMinIndex = find(min(CdAirfoil) == CdAirfoil);
% Calculate angle of attack at minimum drag for wing
alphaWingMinD = alphaAirfoil(CdAirfoilMinIndex);
%calculate the cl at minimum drag
ClMinD = a*(alphaWingMinD-alphaZeroLift);
% calculate the skin friction drag of the whole aircraft
CDmin = Cfe*((Swet)/(wingArea));
% calculate CD0 (derivations in the appendix)
CD01 = ((pi*AR*CDmin)+((1/(0.99*s))*(ClMinD.^2)))./((pi*AR)-(0.38.*pi.*AR.*(ClMinD).^2));
% calculate e0 (derivations in the appendix)
e0 = 1./((1/(0.99*s))+(0.38.*CD01.*pi.*AR));
% Calculate the coefficient of lift
CL = a*(alphaAirfoil(1:end)-alphaZeroLift);
% calculate k1 and k2
k1 = 1/(pi*AR*e0);
k2 = (-2*ClMinD*k1);
% Calculate whole drag
CD = (CD01 + (k1*(CL.^2))+(k2*CL));
% calculate percent difference in the two drags
percent = 100 - ((CDBody./CD(1:18))*100);

% calculate velocity at max range
CLR = (CD01*pi*e0*AR)^(1/2);
Vr = ((2*GTOW*9.81)/(CLR*S*pAlt))^(1/2);


% Calculate the lift over drag of both data
lod = CLBody./CDBody;
lode = CL./CD;
% calculate percent different in lod
pd = abs(100 - (abs((lod./lode(1:18)))*100)); % Calculate percent difference in L/D graphs
%find the cl value where lod is maximum
maxCl = CLBody(find(max(lod) == lod));
% calculate velocity at max range for cfd
Vact = ((2*GTOW*9.81)/(maxCl*S*pAlt))^(1/2);

Rpercent =100 - ((Vr/Vact)*100); % quantify difference

% calculate velocity at max endurance
CLE = (3*CD01*pi*e0*AR)^(1/2);
Ve = ((2*GTOW*9.81)/(CLE*S*pAlt))^(1/2);
lodact = (CLBody.^(3/2))./CDBody;
ClThree = CLBody(find(max(lodact) == lodact));
Vacte = ((2*GTOW*9.81)/(ClThree*S*pAlt))^(1/2);
ePercent = 100 - ((Ve/Vacte)*100);
% calculate angle of attack at max endurance
enduranceangle = alphaBody(find(max(lodact) == lodact)); % find angle of attack for given data

for i = 1: length(CL)
if abs(CL(i)-CLE)<=0.12
maxendangle = alphaBody(i);
i = 18;
end
end
%%
CDWing = CdAirfoil(1:21) + (CL.^2)./(pi*e*AR);

%% Plots
figure (1)
plot(alphaAirfoil, ClAirfoil);
xlabel("Angle of Attack (deg)");
ylabel("Cl");
title("Angle of Attack vs Cl");
xline(0);
yline(0);


figure (2)
plot(ClAirfoil, CdAirfoil);
xlabel("Cl");
ylabel("Cd");
title("Cl vs Cd");

figure (3)
hold on
plot(CL(1:18), CD(1:18));
plot(CLBody, CDBody);
xlabel("CL");
ylabel("CD");
title("CL vs CD");
legend("Calculated CD","Given CD");
hold off

figure (4)
hold on
plot(alphaAirfoil(1:21), CL, 'LineWidth', 2);
plot(alphaAirfoil, ClAirfoil, 'LineWidth', 2);
plot(alphaBody(1:18), CLBody(1:18), 'LineWidth', 2);
xline(0,'--','LineWidth',1.5);
xlabel("Angle of Attack",'FontSize',14);
ylabel("CL"','FontSize',14);
legend("Calculated CL","2-D Wing CL","Given CL CFD","AoA = 0",'FontSize',14);
title("CL vs Angle of Attack",'FontSize',16);
hold off

figure (5)
hold on
plot(CL, CDWing, 'LineWidth', 2);
plot(CL, CD, 'LineWidth', 2);
plot(CLBody,CDBody,'LineWidth', 2);
xlabel("CL",'FontSize',14);
ylabel("CD"','FontSize',14);
legend("Estimated 3-D Finite Wing","Whole Aircraft","CFD Data",'FontSize',14);
title("CD vs CL",'FontSize',16);
hold off

lode = CL./(CD);
lod = (CLBody(1:18)./CDBody(1:18));

figure (6)
hold on
plot(alphaBody(1:18), lode(1:18), 'LineWidth', 2);
plot(alphaBody(1:18), lod, 'LineWidth', 2);
xlabel("Angle of Attack",'FontSize',14);
ylabel("L/D"','FontSize',14);
legend("Calculated L/D","CFD",'FontSize',14);
title("L/D vs Angle of Attack",'FontSize',16);
hold off


figure(7)
hold on
plot(alphaAirfoil(1:21), CD, 'LineWidth', 2);
plot(alphaAirfoil, CdAirfoil, 'LineWidth', 2);
plot(alphaBody(1:18), CDBody(1:18), 'LineWidth', 2);
xline(0,'--','LineWidth',1.5);
xlabel("Angle of Attack",'FontSize',14);
ylabel("CL"','FontSize',14);
legend("CD","2-D Wing CD","Given CD CFD","AoA = 0",'FontSize',14);
title("CD vs Angle of Attack",'FontSize',16);
hold off


