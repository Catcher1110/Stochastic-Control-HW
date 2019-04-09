m = 1; % kg
K = 10/0.01; % N/m
b = 0.1/0.01; % N-sec/m
q = 0.1; % noise strength N^2/sec
dt = 0.001; % time step sec
timesteps = 1000;
ensembles = 1000;
dbeta = randn(timesteps,ensembles)*sqrt(q*dt);
% x,p for monte carlo simulation
x = zeros(timesteps,ensembles);
p = zeros(timesteps,ensembles);
x0 = 0.01;
x(1,1:ensembles) = x0;

for j = 1:ensembles
    for i = 2:timesteps
        x(i,j) = x(i-1,j) + p(i-1,j)/m*dt;
        p(i,j) = p(i-1,j) + (-K*x(i-1,j)-b/m*p(i-1,j))*dt + dbeta(i-1,j);
    end
end

% Xn,Pn,VarN for numerical solution
Xn = zeros(timesteps,1);
Pn = zeros(timesteps,1);
VarN = zeros(2,2,timesteps);
Xn(1, 1) = x0;
F = [0, 1/m; -K, -b/m];
G = [0;1];
for i = 2:timesteps 
    Xn(i) = Xn(i-1) + Pn(i-1)/m*dt;
    Pn(i) = Pn(i-1) + (-K*Xn(i-1)-b/m*Pn(i-1))*dt;
    
    VarN(:,:,i) = VarN(:,:,i-1) + (F*VarN(:,:,i-1)+VarN(:,:,i-1)*F'+G*q*G')*dt;
end

Mean_x = mean(x,2);
Mean_p = mean(p,2);
Var_x = var(x,0,2);
Var_p = var(p,0,2);
Var_XN = VarN(1,1,:);
Var_XN = Var_XN(:);
Var_PN = VarN(2,2,:);
Var_PN = Var_PN(:);
%plot
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
figure(1);
plot(tV, Mean_p, 'b', tV, Pn, 'r');
xlabel('t second');
ylabel('P kg*m/s');
legend('Monte Carlo','ODE solution');
title('mean for x2 (p)');
figure(2);
plot(tV, Mean_x, 'b', tV, Xn, 'r');
xlabel('t second');
ylabel('X m');
legend('Monte Carlo','ODE solution');
title('mean for x1 (x)');
figure(3);
plot(tV, Var_x, 'b', tV, Var_XN,'r');
xlabel('t second');
ylabel('Variance X^2 m^2');
legend('Monte Carlo','ODE solution');
title('variance for x1 (x)');
figure(4);
plot(tV, Var_p, 'b', tV, Var_PN,'r');
xlabel('t second');
ylabel('Variance P^2 (kg*m/s)^2');
legend('Monte Carlo','ODE solution');
title('variance for x2 (p)');
