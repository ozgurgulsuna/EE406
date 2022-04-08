%% EE406 Experiment 2 Calculations
clc;
clear all;
close;

%% Constants
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultLineLineWidth',3)
set(0,'defaultAxesFontSize',15)

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

syms Kp Kv s omega_n zeta

%% Transfer function of the system

numerator = Kg*Kt*Eff_m*Eff_g*r_mp ;
denominator1 = [Lm*((r_mp^2*Mc2)+(Jm*Eff_g*Kg^2)),(Rm*r_mp^2*Mc2)+(Rm*Jm*Eff_g*Kg^2)+(Lm*r_mp^2*Beq),(Km*Kt*Eff_m*Eff_g*Kg^2)+(Rm*r_mp^2*Beq),0]; % Mass of the Cart Only
denominator2 =[Lm*(r_mp^2*(Mc2+Mw)+Jm*Eff_g*Kg^2),Rm*r_mp^2*(Mc2+Mw)+Rm*Jm*Eff_g*Kg^2+Lm*r_mp^2*Beq,Km*Kt*Eff_m*Eff_g*Kg^2+Rm*r_mp^2*Beq,0]; %Total Mass

sys1 = tf(numerator,denominator1);
sys1_m = minreal(sys1);

sys2 = tf(numerator,denominator2);
sys2_m = minreal(sys2);

%%  PV Controller Implementation

G(s) = 2.451/(s^2+17.1*s);  % Open-Loop Transfer Function of the System.
simplify(G(s));

T(s) = (Kp*G(s)/(1+G(s)*s*Kv))/(1+(Kp*G(s)/(1+G(s)*s*Kv))); % Closed-Loop Transfer Function of the Whole System
simplify(T(s));

%% Controller Tuning Kp, Kv

PO = 100*exp((-zeta*pi)/(sqrt(1-zeta^2))) == 10; % Max overshoot, find zeta
zeta=vpasolve(PO,zeta);

t_p = pi/(omega_n*sqrt(1-zeta^2)) == 0.15; % fixed peak time, use zeta
omega_n = vpasolve(t_p,omega_n);

Kp = vpasolve(sqrt(2.451*Kp) == omega_n,Kp);
Kv = vpasolve((17.1+2.451*Kv)/(2*(sqrt(2.451*Kp))) == zeta,Kv);







