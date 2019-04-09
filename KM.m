m = 1; % kg
K = 10/0.01; % N/m
b = 0.1/0.01; % N-sec/m
q = 0.1; % noise strength N^2/sec
dt = 0.001; % time step sec
F = [0, 1/m; -K, -b/m];
H = [1, 0];
G = [0;1];
%R = 0.01^2;
R = 0.001^2;
timesteps = 1000;
ensembles = 100;
dbeta = randn(timesteps,ensembles)*sqrt(q*dt);
V = randn(timesteps, ensembles)*sqrt(R);
x = zeros(timesteps,ensembles);
p = zeros(timesteps,ensembles);
x0 = 0.01;
x(1,1:ensembles) = x0;
xV = zeros(2, timesteps,ensembles);
x_plusV = zeros(2,timesteps,ensembles);
x_minusV = zeros(2,timesteps,ensembles);
x_plusV(1,1,1:ensembles) = x0;
x_minusV(1,1,1:ensembles) = x0;
p_plus = zeros(2, 2, timesteps,ensembles);
p_minus = zeros(2, 2, timesteps,ensembles);
sigma_xV = zeros(timesteps,ensembles);
for j = 1:ensembles
    for i = 2:timesteps
        x(i,j) = x(i-1,j) + p(i-1,j)/m*dt;
        p(i,j) = p(i-1,j) + (-K*x(i-1,j)-b/m*p(i-1,j))*dt + dbeta(i-1,j);       
        xV(:,i,j) = [x(i,j); p(i,j)];
        z = x(i,j) + V(i-1,j);
        
        x_minusV(:,i,j) = x_plusV(:,i-1,j) + F*x_plusV(:,i-1,j)*dt;
        p_minus(:,:,i,j) = p_plus(:,:,i,j) + (F*p_plus(:,:,i-1,j)+p_plus(:,:,i-1,j)*F'+G*q*G')*dt;
        KG = p_minus(:,:,i,j)*H'*inv(H*p_minus(:,:,i,j)*H'+R);
        x_plusV(:,i,j) = x_minusV(:,i,j)+KG*(z-H*x_minusV(:,i,j));
        p_plus(:,:,i,j) = p_minus(:,:,i,j)-KG*H*p_minus(:,:,i,j);
        p_minus(:,:,i,j) = p_plus(:,:,i,j);
        sigma_xV(i,j) = sqrt(p_plus(1,1,i,j));
    end
end

Mean_x_plusV = mean(x_plusV, 3);
Mean_xV = mean(x,2);
abs_error = abs(xV - x_plusV);
%plot
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
figure(1);
plot(tV, Mean_x_plusV(1,:), 'b', tV, Mean_xV, 'r');
xlabel('t second');
ylabel('X m');
legend('estimation displacement','actual displacement');
title('R = 0.001^2 m^2');

