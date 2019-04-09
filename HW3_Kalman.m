m = 1; % kg
K = 10/0.01; % N/m
b = 0.1/0.01; % N-sec/m
q = 0.1; % noise strength N^2/sec
dt = 0.0001; % time step sec
F = [0, 1/m; -K, -b/m];
H = [1, 0];
G = [0;1];
R = 0.01^2;
%R = 0.001^2;
timesteps = 10000;
ensembles = 1000;
dbeta = randn(timesteps,1)*sqrt(q*dt);
V = randn(timesteps, 1)*sqrt(R);
x = zeros(timesteps,1);
p = zeros(timesteps,1);
x0 = 0.01;
x(1,1) = x0;
xV = zeros(2, timesteps);
xV(1,1) = x0;
x_plusV = zeros(2, timesteps);
x_minusV = zeros(2, timesteps);
x_plusV(1,1) = x0;
x_minusV(1,1) = x0;
p_plus = zeros(2, 2, timesteps);
p_minus = zeros(2, 2, timesteps);
sigma_xV = zeros(timesteps);
for i = 2:timesteps
    x(i) = x(i-1) + p(i-1)/m*dt;
    p(i) = p(i-1) + (-K*x(i-1)-b/m*p(i-1))*dt + dbeta(i-1);
    
    xV(:,i) = [x(i); p(i)];
    z = x(i) + V(i-1);
    
    x_minusV(:,i) = x_plusV(:,i-1) + F*x_plusV(:,i-1)*dt;
    p_minus(:,:,i) = p_plus(:,:,i-1) + (F*p_plus(:,:,i-1)+p_plus(:,:,i-1)*F'+G*q*G')*dt;
    KG = p_minus(:,:,i)*H'*inv(H*p_minus(:,:,i)*H'+R);
    x_plusV(:,i) = x_minusV(:,i)+KG*(z-H*x_minusV(:,i));
    p_plus(:,:,i) = p_minus(:,:,i)-KG*H*p_minus(:,:,i);
    p_minus(:,:,i) = p_plus(:,:,i);
    sigma_xV(i) = sqrt(p_plus(1,1,i));
end

abs_error = abs(xV(1,:) - x_plusV(1,:));
abs_error2 = p_plus(1,1,:);
abs_error2 = abs_error2(:);
%plot
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
figure(1);
plot(tV, x_plusV(1,:), 'b', tV, xV(1,:), 'r');
xlabel('t second');
ylabel('X m');
legend('estimation displacement','actual displacement');
title('R = 0.01^2 m^2');
figure(2);
plot(tV, abs_error);
xlabel('t second');
ylabel('X m');
title('estimation error of x');
figure(3);
plot(tV, abs_error2);
xlabel('t second');
ylabel('X m');
title('estimation error variance of x');

