k = 10/0.01; % N/m
m = 1; % kg
b = 0.1/0.01; % N-sec/m
T = 0.5; % sec
sigma = 1; % N
syms x1(t) x2(t) x3(t)
X = [x1;x2;x3];
F = [0,1,0;-k/m,-b/m,1/m;0,0,-1/T];
xeqns = diff(X,t)==F*X;
xcond = [x1(0)==0.01,x2(0)==0,x3(0)==0];
xsol = dsolve(xeqns,xcond);
xsol1(t) = xsol.x1;

syms p11(t) p12(t) p13(t) p21(t) p22(t) p23(t) p31(t) p32(t) p33(t)
P = [p11,p12,p13;p21,p22,p23;p31,p32,p33];
G = [0;0;1];
q = 2*sigma^2/T;
RHS = F*P + P*F' + G*q*G';
peqns = diff(P,t) == RHS;
pcond = [p11(0)==0,p12(0)==0,p13(0)==0,...
    p21(0)==0,p22(0)==0,p23(0)==0,...
    p31(0)==0,p32(0)==0,p33(0)==0,];
psol = dsolve(peqns,pcond);
psol1(t) = psol.p11;
% 
% dt = 0.01;
% timesteps = 200;
% tV = 0:dt:(timesteps-1)*dt;
% tV = tV';
% xV = zeros(1,timesteps);
% pV = zeros(1,timesteps);
% for i = 1:1:timesteps
%     xV(i) = xsol1(dt*(i-1));
% %     pV(i) = psol1(dt*(i-1));
% end
% figure(1);
% plot(tV, xV);
% xlabel('t second');
% ylabel('X m');
% legend('displacement');
% title('Displacement From the ODE');

% figure(2);
% plot(tV, pV);
% xlabel('t second');
% ylabel('variance m^2');
% legend('Variance');
% title('Variance From the ODE');

dt = 0.0001;
timesteps = 20000;
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
% Xn,Pn,VarN for numerical solution
Xn = zeros(3,timesteps);
VarN = zeros(3,3,timesteps);
Xn(1, 1) = 0.01;

for i = 2:timesteps 
    Xn(:,i) = Xn(:,i-1) + F*Xn(:,i-1)*dt;   
    VarN(:,:,i) = VarN(:,:,i-1) + (F*VarN(:,:,i-1)+VarN(:,:,i-1)*F'+G*q*G')*dt;
end

Var_XN = VarN(1,1,:);
Var_XN = Var_XN(:);

xV = Xn(1,:);
xV = xV(:);
figure(2);
plot(tV, Var_XN);
xlabel('t second');
ylabel('variance m^2');
legend('Variance');
title('Variance From the ODE');


