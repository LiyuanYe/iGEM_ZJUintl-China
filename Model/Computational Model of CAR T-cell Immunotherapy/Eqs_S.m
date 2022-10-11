function df=Eqs_S(t,f, rF, rE, lE, lM, Ke, rA, Kf, Kp, Ka);

%-----ODEs-----
%f(1)=nF   (Fibrotic cells)
%f(2)=nE   (Effector CAR-T cells)
%f(3)=nM   (Memory CAR-T cells)

df=zeros(3,1);

df(1)=rF-Ke*f(1)/(f(1)+Kf)*f(2);

df(2)=rE*f(1)/(f(1)+Kp)*f(2)+rA*f(1)/(f(1)+Ka)*f(3)-lE*f(2);

df(3)=-rA*f(1)/(f(1)+Ka)*f(3)-lM*f(3);