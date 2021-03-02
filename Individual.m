clc
clear 

CL = [0,0.1,0.22,0.35, 0.47, 0.6, 0.65, 0.7, 0.78, 0.8];
CD = [0.02,0.023,0.03,0.04,0.05,0.06, 0.08, 0.11, 0.13,0.15];
AoA = [1,2,3,4,5,6,7,8,9,10];
AR =14;
e =  1.78*(1-0.045*AR^0.68)-0.64;
a0 = 0.1;
a = a0 * (a0/(1+((57.3*a0)/(pi*e*AR))));
CL3D = a*AoA;
CDwing = CD + ((CL3D.^2)./(pi*e*AR));

AreaF = 0.637;
AreaW = 2/AR;
AreaT = 0.104;


[~,idxMin] = min(CDwing);
CLminD = CL3D(idxMin);

% Hard coded Values
Cfe = 0.003; % equivalent skin friction coefficient
Swet = AreaF + AreaW + AreaT; % [m^2] wetted area (Tempest)
Sref = 1/AR; % [m^2] projected area of the wing (Tempest)

% minimum drag (eq. 10)
CDmin = Cfe*Swet/Sref;

% arbitrary constant (eq. 5)
k1 = 1/(pi*e*AR);

% whole aircraft drag (eq. 6)
CDA = CDmin + k1*(CL3D - CLminD).^2;

CLCD = CL3D./CDA;

CLCDMax = 0.5 * sqrt((pi *e*AR)/CDmin)


%Weight
Density = 0.295;
WeightF = 0.0207 * Density;
WeightW = AreaW *0.015 * Density;
WeightT = AreaT * 0.01 * Density;
WeightAdd = 0.4

WeightTot = WeightT + WeightW + WeightF + WeightAdd;

Range = 7 * (CLCDMax);

time = 8.4; %sec

Vsink = 7 /time;
Vinf = WeightTot * Vsink / CDmin

figure(1)
hold on
plot( AoA, CL)
plot(AoA, CD )
