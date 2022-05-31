%% EE406 Experiment 4 Calculations

clc;
clear all;
close;



%% Symbolic definitions 

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