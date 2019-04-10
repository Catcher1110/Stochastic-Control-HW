m = 1; % kg
K = 10/0.01; % N/m
b = 0.1/0.01; % N-sec/m
T = 0.5; % sec
sigma = 1; % N
q = 2*sigma^2/T; % noise strength N^2/sec
dt = 0.0001; % time step sec
F = [0, 1, 0; -K/m, -b/m, 1/m; 0, 0, 0];
H = [1, 0, 0];
G = [0; 0; 1];
R = 0.001^2;
timesteps = 20000;

dbeta = randn(timesteps,1)*sqrt(q*dt);
V = randn(timesteps, 1)*sqrt(R);
x1 = zeros(timesteps,1);
x2 = zeros(timesteps,1);
x3 = zeros(timesteps,1);
x0 = 0.01;
x1(1,1) = x0;
xV = zeros(3, timesteps);
xV(1,1) = x0;
x_plusV = zeros(3, timesteps);
x_minusV = zeros(3, timesteps);
x_plusV(1,1) = x0;
x_minusV(1,1) = x0;
p_plus = zeros(3, 3, timesteps);
p_minus = zeros(3, 3, timesteps);
sigma_xV = zeros(timesteps);
for i = 2:timesteps
    x1(i) = x1(i-1) + x2(i-1)*dt;
    x2(i) = x2(i-1) + (-K/m*x1(i-1)-b/m*x2(i-1)+1/m*x3(i-1))*dt;
    x3(i) = x3(i-1) + (-1/T * x3(i-1)) + dbeta(i-1);
    
    xV(:,i) = [x1(i); x2(i); x3(i)];
    z = x1(i) + V(i-1);
    
    x_minusV(:,i) = x_plusV(:,i-1) + F*x_plusV(:,i-1)*dt;
    p_minus(:,:,i) = p_plus(:,:,i-1) + (F*p_plus(:,:,i-1)+p_plus(:,:,i-1)*F'+G*q*G')*dt;
    KG = p_minus(:,:,i)*H'*inv(H*p_minus(:,:,i)*H'+R);
    x_plusV(:,i) = x_minusV(:,i)+KG*(z-H*x_minusV(:,i));
    p_plus(:,:,i) = p_minus(:,:,i)-KG*H*p_minus(:,:,i);
    p_minus(:,:,i) = p_plus(:,:,i);
    sigma_xV(i) = sqrt(p_plus(1,1,i));
end

%plot
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
figure(1);
plot(tV, x_plusV(1,:), 'b', tV, xV(1,:), 'r');
xlabel('t second');
ylabel('X m');
legend('estimation displacement','actual displacement');
title('Displacement with R = 0.001^2 m^2');

figure(2);
plot(tV, x_plusV(2,:), 'b', tV, xV(2,:), 'r');
xlabel('t second');
ylabel('V m/s');
legend('estimation velocity','actual velocity');
title('Velocity with R = 0.001^2 m^2');

figure(3);
plot(tV, x_plusV(3,:), 'b', tV, xV(3,:), 'r');
xlabel('t second');
ylabel('Force N');
legend('estimation Force','actual Force');
title('Force with R = 0.001^2 m^2');
