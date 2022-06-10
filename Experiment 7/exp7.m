s = tf('s');
set(0,'DefaultLineLineWidth',0.6)
H = 2.45/(s*(s+17.13));
figure;
bode(H);
hold on
%% T = 0.001
T = 0.001;
Hz = c2d(H,T,'zoh');
bode(Hz);
%% T = 0.01
T = 0.01;
Hz = c2d(H,T,'zoh');
bode(Hz);
%% T = 0.1
T = 0.1;
Hz = c2d(H,T,'zoh');
bode(Hz);
%% T = 1
T = 1;
Hz = c2d(H,T,'zoh');
bode(Hz);

set(0,'defaultTextInterpreter','latex')
title("\bf{Zero Order Hold Transformation for Different Sampling Rates}",'FontSize',14);

h = findobj(gcf,'type','line');
set(h,'linewidth',2);

set(0,'defaultAxesFontSize',12)
%legend('T = 0');
%legend('T = 0.001');
%legend('T = 0','T = 0.001','T = 0.01','T = 0.1','T = 1');

set(gcf,'Position',[0 0 1200 600]);
grid on

%%
%% Tustin, Bilinear, Trapezoidal
figure;
bode(H);
hold on
%% T = 0.001
T = 0.001;
Hz = c2d(H,T,'tustin');
bode(Hz);
%% T = 0.01
T = 0.01;
Hz = c2d(H,T,'tustin');
bode(Hz);
%% T = 0.1
T = 0.1;
Hz = c2d(H,T,'tustin');
bode(Hz);
%% T = 1
T = 1;
Hz = c2d(H,T,'tustin');
bode(Hz);


set(0,'defaultTextInterpreter','latex')
title("\bf{Tustin Transformation for Different Sampling Rates}",'FontSize',14);

h = findobj(gcf,'type','line');
set(h,'linewidth',2);

set(0,'defaultAxesFontSize',12)
%legend('T = 0');
%legend('T = 0.001');
%legend('T = 0','T = 0.001','T = 0.01','T = 0.1','T = 1');

set(gcf,'Position',[0 0 1200 600]);
grid on

%%


%% Pole zero matching
figure;
bode(H);
hold on
%% T = 0.001
T = 0.001;
Hz = c2d(H,T,'matched');
bode(Hz);
%% T = 0.01
T = 0.01;
Hz = c2d(H,T,'matched');
bode(Hz);
%% T = 0.1
T = 0.1;
Hz = c2d(H,T,'matched');
bode(Hz);
%% T = 1
T = 1;
Hz = c2d(H,T,'matched');
bode(Hz);


set(0,'defaultTextInterpreter','latex')
title("\bf{Pole Zero Matching Transformation for Different Sampling Rates}",'FontSize',14);

h = findobj(gcf,'type','line');
set(h,'linewidth',2);

set(0,'defaultAxesFontSize',12)
%legend('T = 0');
%legend('T = 0.001');
%legend('T = 0','T = 0.001','T = 0.01','T = 0.1','T = 1');

set(gcf,'Position',[0 0 1200 600]);
grid on

%%
s = tf('s');

Kp=275.6;
Kd = 5.55;

Gc=Kp+Kd*s;

Gc_z=c2d(Gc,0.02,'zoh')

275.6+5.55*(1-1/s)

%%
s = tf('s');
z = tf('z');
H = 2.45/(s*(s+17.13));

for T = [1/100 1/100]
    G_z = c2d(H,T,'zoh');
    [num, den]=tfdata(G_z);
    den=cell2mat(den);
    num=cell2mat(num);
    A= den;
    B= num;
    Gc_z1=tf(A,B,T);
end

Gc_z2=tf(1,[1 -1],T);
Gc_z=Gc_z1*Gc_z2

for T = [1/100 1/100]
    G_z = c2d(H,T,'zoh');
    [num, den]=tfdata(G_z);
    den=cell2mat(den);
    num=cell2mat(num);
    A= den;
    B= -num;
    B(3)=B(3)+num(1)+num(2)+num(3);
    Gc_z=tf(A,B,T)
end



