clc; clear; close all;

[x, y] = meshgrid(linspace(-100, 100, 20), linspace(-100, 100, 20));

k = 0.05;
d = (y - tan(pi/3).*x)./sqrt(1 + tan(pi/3)^2);
chi = pi/3 - atan(k*d); %to reverse the direction: chi = pi + atan(ky);
u = cos(chi);
v = sin(chi);

figure;
hold on;
quiver(x, y, u, v, 1);
plot([-100, 100]./tan(pi/3), [-100, 100], linewidth=1.5);
axis equal;
grid on; grid minor;
xlabel("x [m]"); ylabel("y [m]");
title("Vector Field across the line");
hold off;

r = 50;
d = sqrt(x.^2 + y.^2);
gamma = atan2(y, x);
chi_d = gamma - pi/2 - atan(k*(d - r));     %to change the direction: chi_c = gamma + pi/2 + atan(k*(d-r));


figure;
hold on;
quiver(x, y, cos(chi_d), sin(chi_d), 1);

plot(r*cos(linspace(0, 2*pi, 100)), r*sin(linspace(0, 2*pi, 100)), LineWidth=1.5);
axis equal;
grid on; grid minor;
xlabel("x [m]"); ylabel("y [m]");
title("Vector Field across the Loiter circle");
hold off;
