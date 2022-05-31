%% EE406 Experiment 3 Calculations

clc;
clear all;
close;

set(0,'defaultTextInterpreter','latex')
set(0,'DefaultLineLineWidth',3)
set(0,'defaultAxesFontSize',15)

%% Symbolic definitions 

syms  Kg Kt Eff_m Eff_g r_mp Lm  Mc2 Mw Jm Rm Beq Km Vm V s

%% Values of the Parameters

Kg =  3.71; % Planetray Gear Ratio (1/1)
Kt = 0.00767 ; % Motor Torque Constant (Nm/A)
Eff_m = 1; % Motor Efficiency (1/1)
Eff_g = 1; % Planetary Gearbox Efficiency (1/1)
r_mp = 6.35E-3 ; % Motor Pinion Radius (m)
Lm = 0 ; % (0.18E-3 - ignored) Motor Armature Inductance (H)
Mc2 = 0.57 ; % IP02 Cart Mass (kg)
Mw = 0.37 ; % IP02 Cart Weight Mass (kg)
Jm = 3.90E-7 ; % Rotor Moment of Intertia (kgm^2)
Rm = 2.6 ; % Motor Armature Resistance (Ohms)
Beq = 4.3 ; % Equivalent Viscous Damping Coefficient as seen at the Motor Pinion (Ns/m)
Km = 0.00767 ; % Back-ElectroMotive-Force (EMF) constant (Vs/rad)

%% Simplification of the plant (Voltage to Velocity)

G2(s) = (1/(s*Mc2))/(1+(Beq*(1/(s*Mc2))));
simplify(G2(s));

G1(s) = (Eff_g*Kg/r_mp*G2(s))/(1+((Eff_g*Kg/r_mp*G2(s))*Jm*Kg*s/r_mp));
G1(s) =simplify(G1(s));

G(s) = ((G1(s))*Kt*Eff_m/(Rm+s*Lm))/(1+(((G1(s))*Kt*Eff_m/(Rm+s*Lm))*Km*Kg/r_mp));
G(s) = simplify(G(s));


%% Monic Form Conversion

G(s) =subs(G(s));

[num, denum] = (numden(G(s)));
Coeff = coeffs(denum);
Coeff = fliplr(Coeff);
G(s) = (num/Coeff(1))/(denum/(Coeff(1)))


%% Simulation for Steady State Errors   
clc;
clear all;
close;

s = tf('s');
G_ol = 2.4513/(s+17.1001); % Open Loop
G_oli = 2.4513/(s*(s+17.1001)); % Open Loop  with Integrator
G_olig = 2660*2.4513/(s*(s+17.1001)); % Open Loop  with Integrator and Gain
G_uf = G_ol/(1+G_ol); % Unity Feedback
G_ufi = (G_ol/s)/(1+(G_ol/s)); % Unity Feedback with Integrator
G_c=(6.6514*(s+80/6.6514))/(s+80*6.6514); % calculated compansator
G_olcomp = G_c*G_olig; % open loop tf with compensator
G_cl = G_olcomp/(1+G_olcomp); %Closed loop tf with compensator

P = pole(G_ol);
Z = zero(G_ol);

Gain_ol = 2.4513/17.1001;

%% Input Function
[square,t] = gensig("square",10,100,0.01);

impulse = t==0;
unitstep = t>=0;
ramp = t.*unitstep;
quad = t.^2.*unitstep;

input = unitstep;
%% Simulation
lsim(G_ol,input,t)
hold on
lsim(G_uf,input,t)
lsim(G_ufi,input,t)

syms s

[Num,Den] = tfdata(G_oli);
sys_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);

eqn =  simplify(sys_syms) ==1
solve(eqn,s)
figure;
bode(G_oli)

figure
bode(G_olig)

%%  Integrator as H(s)

lsim(G_uf/s,unitstep,t);

%%
syms s
s = tf('s');

minreal(G_olcomp);
zpk(G_olcomp)

E_step=1/(1+G_olcomp)
E_ramp=(1/s)*1/(1+G_olcomp)
E=1/(1+G_olcomp)
zpk(E_ramp)
zpk(E_step)

minreal(G_cl)
zpk(minreal(G_cl))



