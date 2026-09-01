clc; clear; close all;

simTime = 50;
dt = 0.001;

%%Path Specification
o = [0, 0]';   %m
r = 150;        %m
s = 1; %-1: clockwise; 1: anticlockwise

%%UAV Specifications
va = 25; %m/s

p0 = [-350, -555]';   %m
psi0 = -3*pi/4;        %radians

k1 = 1;
k2 = 0.2;

g = 9.81;
phiMax = 45*pi/180; %radians
dpsiMax = g*tan(phiMax)/va; %radians/sec

state = zeros([6, simTime/dt + 1]);
state(1:3, 1) = [p0', psi0]';

i = 2;
for t=dt:dt:simTime
    p = state(1:2, i - 1); x = p(1); y = p(2);
    psi = state(3, i - 1);
    
    d = sqrt((x - o(1))^2 + (y - o(2))^2) - r;
    theta = atan2(y - o(2), x - o(1));
    thetad = theta + s*(pi/2);
    
    dpsi = k1*wrapToPi(thetad - psi) + s*atan(k2*d);

    psi = wrapToPi(psi + dpsi*dt);
    x = x + va*cos(psi)*dt;
    y = y + va*sin(psi)*dt;


    state(:, i) = [x, y, psi, dpsi, wrapToPi(thetad - psi), d]';

    i = i + 1;
end


figure;
hold on;
grid on; grid minor;
axis equal;
x = (o(1) - r):0.001:(o(1) + r);
y1 = o(2) - sqrt(r^2 - (x - o(1)).^2);
y2 = o(2) + sqrt(r^2 - (x - o(1)).^2);
x = [x(1:end-1), flip(x)]';
y = [y1(1:end-1), flip(y2)]';
plot(x, y,"--", lineWidth=1, Color="#ff0000");
plot(p0(1), p0(2), "Marker","o", "MarkerSize", 8, 'LineWidth', 1.5, "DisplayName", "Inital Position");
plot(state(1, :), state(2, :), linewidth = 1.5, color='#0099aa');
plot(state(1, end), state(2, end), '^', markersize = 8, color = '#0099ff', linewidth=1.5);
xlabel("x-axis");
ylabel("y-axis");
title("Loitering with NLGL");
hold off;

figure; 
subplot(2, 2, 1);
plot(0:dt:simTime, state(3, :)*180/pi,  "LineWidth", 1.5);
hold on; grid on; grid minor;
xlabel("Time [s]");
ylabel("Heading Angle, \psi [Degree]");
hold off;

subplot(2, 2, 2);
plot(0:dt:simTime, state(4, :),  "LineWidth", 1.5);
hold on; grid on; grid minor;
xlabel("Time [s]");
ylabel("Heading Angle Rate, \psi_D_O_T [rad/s]");
hold off;

subplot(2, 2, 3);
plot(0:dt:simTime, state(5, :)*180/pi,  "LineWidth", 1.5);
hold on; grid on; grid minor;
xlabel("Time [s]");
ylabel("Angle Error, \eta [Degree]");
hold off;

subplot(2, 2, 4);
plot(0:dt:simTime, state(6, :),  "LineWidth", 1.5);
hold on; grid on; grid minor;
xlabel("Time [s]");
ylabel("Cross Track Error");
hold off;