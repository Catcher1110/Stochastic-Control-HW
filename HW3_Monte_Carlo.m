m = 1; % kg
K = 10/0.01; % N/m
b = 0.1/0.01; % N-sec/m
q = 0.1; % noise strength N^2/sec
dt = 0.001; % time step sec
timesteps = 1000;
ensembles = 1000;
dbeta = randn(timesteps,ensembles)*sqrt(q*dt);
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

Mean_x = mean(x,2);
Mean_p = mean(p,2);
Var_x = var(x,0,2);
Var_p = var(p,0,2);

%plot
tV = 0:dt:(timesteps-1)*dt;
tV = tV';
y = 0.01 - exp(-tV.*5).*sin(tV.*sqrt(975))/2/sqrt(975);
figure(1);
plot(tV, Var_x, 'b', tV, Var_p, 'r');
xlabel('t second');
ylabel('Variance');
figure(2);
plot(tV, Mean_x, 'b', tV, y, 'r');
legend('');
