%% EE406 Experiment 5 Calculations

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
Ks = 160 ;
Kt = 0.00767 ; % Motor Torque Constant (Nm/A)
Eff_m = 1; % Motor Efficiency (1/1)
Eff_g = 1; % Planetary Gearbox Efficiency (1/1)
r_mp = 6.35E-3 ; % Motor Pinion Radius (m)
Lm = 0 ; % (0.18E-3 - ignored) Motor Armature Inductance (H)
M2 = 0.5425 ; % IP02 Cart Mass (kg)
Mc = 1.1456 ; % IP02 Cart Mass (kg)
Mw = 0.37 ; % IP02 Cart Weight Mass (kg)
Jm = 3.90E-7 ; % Rotor Moment of Intertia (kgm^2)
Rm = 2.6 ; % Motor Armature Resistance (Ohms)
Beq = 5.4 ; % Equivalent Viscous Damping Coefficient as seen at the Motor Pinion (Ns/m)
Beq2 = 1.1;
Km = 0.00767 ; % Back-ElectroMotive-Force (EMF) constant (Vs/rad)
%% Symbolic definitions 

A = [0 0 1 0; 0 0 0 1; (-Ks/Mc) (Ks/Mc) ((-Beq/Mc)-((Kg^2)*Kt*Km/Mc/Rm/r_mp^2)) 0; (Ks/M2) (-Ks/M2) 0 (-Beq2/M2)];
B =[0;0;((Kg*Kt)/(Mc*Rm*r_mp)); 0];

eig(A)
eig(B)

[b,a] = ss2tf(A,B);
sys=tf(b,a);


C=[1 0 0 0 ;0 0 0 0; 0 0 0 0; 0 0 0 0]

O=[C*A ; C*A; C*A^2; C*A^3]

figure;
rlocus(sys);
title("\bf{Pole-Zero Map}",'FontSize',16);
set(gcf,'Position',[0 0 1200 520]);
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultLineLineWidth',0.6)
set(0,'defaultAxesFontSize',14)
grid on


%%

%TF = tf([25],[1 -4 4]);
TF=zpk([-7 7],[-0.1 -11 -0.17+5j -0.17-5j],-1);
rlocus(TF);
pzplot(TF);
step(TF);



figure;
step(sys);
title("\bf{Step Response}",'FontSize',16);
xlim([0 3000]);
ylim([0 500]);
set(gcf,'Position',[0 0 1200 520]);
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultLineLineWidth',0.6)
set(0,'defaultAxesFontSize',14)
grid on









Atilda = A-B*K;
sys = ss(Atilda,B,C,D);
step(sys)
stepinfo(sys)