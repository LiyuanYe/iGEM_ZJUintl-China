%Simulation of CD26 CAR-T's effect on senescent fibroblasts

f0=[2200.24,0,16.46]; % Initial [nF,nE,nM] need revision

rF=20;
rE=1.62;
lE=0.12;
lM=0.00003;
Ke=22.72;
rA=0.65;
Kf=6000;
Kp=600;
Ka=1800;

[t,f]=ode45(@Eqs_S,[0:0.1:90],f0,[], rF, rE, lE, lM, Ke, rA, Kf, Kp, Ka);

%figure;
subplot(2,2,1)
plot(t,f(:,1));
title('nF');
hold on

subplot(2,2,2)
plot(t,f(:,2));
title('nE');
hold on

subplot(2,2,3)
plot(t,f(:,3));
title('nM');




