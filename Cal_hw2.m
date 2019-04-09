% This is for (a)
% m = 60;
% zeta = 0.5;
% fn = 4;
% wn = 2*pi*fn;

% This is for (b)
m = 272;
zeta = 0.10;
fn = 26;
wn = 2*pi*fn;

% This is for (a)
Fa = @(w)((w>=0 & w<=25*2*pi)/2/pi*10^5)*1/m^2./((wn^2-w.^2).^2+(2*zeta*wn*w).^2);
% This is for (b)
Fb = @(w)((w>=0 & w<=10*2*pi)/2/pi*4*10^6 + ...
    (w>=10*2*pi & w<=20*2*pi)/2/pi.*(2.*w*10^5+2*10^6) + ...
    (w>=20*2*pi & w<=100*2*pi)/2/pi*6*10^6 + ...
    (w>=100*2*pi & w<=200*2*pi)/2/pi.*(-6.*w*10^4+12*10^6))...
    .*1/m^2./((wn^2-w.^2).^2+(2*zeta*wn*w).^2);

% This is for (a)
Ex2 = 2/(2*pi) * integral(Fb,0,inf);
p = 1 - normcdf(0.012,0,sqrt(Ex2));
% This is for (b)
% Ex2 = 2/(2*pi) * integral(Fb,0,inf);
% p = 1 - normcdf(0.012,0,sqrt(Ex2));
disp("The probablity is: ");
disp(2*p);