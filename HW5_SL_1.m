m = 1; % kg
K1 = 10/0.01; % N/m
K3 = 10/0.01^3; % N/m^3
b = 0.1/0.01; % N-sec/m
q = 0.1; % noise strength N^2/sec
dt = 0.0001; % time step sec
R = 0.001^2;% measure noise
H = [1, 0];
G = [0; 1];
timesteps = 30000;

dbeta = randn(timesteps,1)*sqrt(q*dt);
V = randn(timesteps, 1)*sqrt(R);
x1 = zeros(timesteps,1);
x2 = zeros(timesteps,1);
x0 = 0.01;
x1(1,1) = x0;
xV = zeros(2, timesteps);
xV(1,1) = x0;
x_plusV = zeros(2, timesteps);
x_minusV = zeros(2, timesteps);
x_plusV(1,1) = x0;
x_minusV(1,1) = x0;
p_plus = zeros(2, 2, timesteps);
p_minus = zeros(2, 2, timesteps);
p_minus(2,2,1) = 0.01*0.01^2;% initial P22, (m/sec)^2
sigma_xV = zeros(timesteps);
for i = 2:timesteps
    x1(i) = x1(i-1) + x2(i-1)*dt;
    x2(i) = x2(i-1) + (-K1/m*x1(i-1)-K3/m*x1(i-1)^3-b/m*x2(i-1))*dt + dbeta(i-1);
    
    xV(:,i) = [x1(i); x2(i)];
    z = x1(i) + V(i-1);
    % update the nonliner N
    nnew = 3*x_plusV(1,i-1)^2+3*p_plus(1,1,i-1);
    N = [0, 1; -K1-nnew, -b/m];
    % estimated m&P value before measurement
    x_minusV(:,i) = x_plusV(:,i-1) + N*x_plusV(:,i-1)*dt;
    p_minus(:,:,i) = p_plus(:,:,i-1) + (N*p_plus(:,:,i-1)+p_plus(:,:,i-1)*N'+G*q*G')*dt;
    % calculate the Kalman Gain Factor
    KG = p_minus(:,:,i)*H'*inv(H*p_minus(:,:,i)*H'+R);
    % update the m&P
    x_plusV(:,i) = x_minusV(:,i)+KG*(z-H*x_minusV(:,i));
    p_plus(:,:,i) = p_minus(:,:,i)-KG*H*p_minus(:,:,i);
    % store the P
    p_minus(:,:,i) = p_plus(:,:,i);
    % calculate the sigma
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
title('Displacement with x0 = 1 cm');

figure(2);
plot(tV, x_plusV(2,:), 'b', tV, xV(2,:), 'r');
xlabel('t second');
ylabel('V m/s');
legend('estimation velocity','actual velocity');
title('Velocity with x0 = 1 cm');
