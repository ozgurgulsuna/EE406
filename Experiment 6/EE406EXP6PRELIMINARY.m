%% EE406 Experiment 6 Calculations
clc;
clear all;
close;

set(0,'defaultTextInterpreter','latex')
set(0,'DefaultLineLineWidth',3)
set(0,'defaultAxesFontSize',15)

%% Symbolic definitions 
syms  xc alpha xc_dot alpha_dot Vm 

%% Values of the Parameters
Jm = 3.90E-7 ; % Rotor Moment of Intertia (kgm^2)
Rm = 2.6 ; % Motor Armature Resistance (Ohms)
Kt = 0.00767 ; % Motor Torque Constant (Nm/A)
Eff_m = 1 ; % Motor Efficiency (1/1)
Km = 0.00767 ; % Back-ElectroMotive-Force (EMF) constant (Vs/rad)
Kg =3.71 ; % Planetray Gear Ratio (1/1)
g = 9.81 ; % Gravitational Constant on Earth (m/s^2)
Eff_g = 1 ; % Planetary Gearbox Efficiency (1/1)
r_mp = 6.35E-3 ; % Motor Pinion Radius (m)
Beq = 4.3 ; % Equivalent Viscous Damping Coefficient as seen at the Motor Pinion (Ns/m)
Mp = 0.230 ; % Long Pendulum Mass (with T-fitting) (kg)
lp = 0.3302 ; % Long Pendulum Length from Pivot to Center Of Gravity (m)
Ip = 7.88e-3; % Long Pendulum Moment of Inertia, about its Center Of Gravity (kg*m^2)
Bp = 0.0024 ; % Viscous Damping Coefficient, as seen at the Pendulum Axis (n*m*s/rad)
Mc = 0.57 + Jm*Kg^2/(r_mp^2)  ; % or 0.38 for cart, inertia addition, + 0.230 for Mp, 
%% System ... 
xc_ddot = (-(Ip+Mp*lp^2)*Beq*(xc_dot)-(Mp^2*lp^3+Ip*Mp*lp)*sin(alpha)*alpha_dot^2 ...
    +(Ip+Mp*lp^2)*(-(Eff_g*Kg^2*Eff_m*Kt*Km*xc_dot)/(Rm*r_mp^2)+(Eff_g*Kg*Eff_m*Kt*Vm)/(Rm*r_mp)) ...
    -Mp*lp*cos(alpha)*Bp*alpha_dot+Mp^2*lp^2*g*cos(alpha)*sin(alpha))/((Mc+Mp)*Ip+Mc*Mp*lp^2+Mp^2*lp^2*sin(alpha)^2);

alpha_ddot = ((Mc+Mp)*Mp*g*lp*sin(alpha)-(Mc+Mp)*Bp*alpha_dot ...
    +(-(Eff_g*Kg^2*Eff_m*Kt*Km*xc_dot)/(Rm*r_mp^2)+(Eff_g*Kg*Eff_m*Kt*Vm)/(Rm*r_mp))*Mp*lp*cos(alpha) ...
    -Mp^2*lp^2*sin(alpha)*cos(alpha)*alpha_dot^2-Mp*lp*cos(alpha)*Beq*xc_dot)/((Mc+Mp)*Ip+Mc*Mp*lp^2+Mp^2*lp^2*sin(alpha)^2);
% Evaluate Jacobian
A = jacobian([xc_dot,alpha_dot,xc_ddot,alpha_ddot],[xc; alpha; xc_dot; alpha_dot]);
B = jacobian([xc_dot,alpha_dot,xc_ddot,alpha_ddot],[Vm]);
% Equilibrium point [0 0 0 0]
xc=0;
xc_dot = 0;
alpha = 0;
alpha_dot = 0;
Vm=0;

A = double(subs(A))
B = double(subs(B))
C = [1 0 0 0 ; 0 1 0 0];
D = [0;0];





A = [0 0 1 0; 0 0 0 1; 0 1.5216 -11.6513 0.0049; 0 -26.1093 26.8458 -0.0841];
B =[0;0;1.5304; -3.5261];
C = [1 0.6413 0 0];
D = 0;
[b,a] = ss2tf(A,B,C,D);
sys=tf(b,a);

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








K = place(A,B,[(-1.8385+1.8385i) (-1.8385-1.8385i) -3+1i -3-1i]);
K = round(K,2)
Atilda = A-B*K;
sys = ss(Atilda,B,C,D);
step(sys)
stepinfo(sys)