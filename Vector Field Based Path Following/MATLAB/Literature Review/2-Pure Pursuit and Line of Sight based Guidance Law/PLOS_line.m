%this code does not deal with the vertical line following edge case
clear;
clc;
close all;

simTime = 30;
dt = 0.001;

w1 = [-200, -200];
w2 = [150, 400];

lineSlope = (w2(2) - w1(2))/(w2(1) - w1(1));
lineIntercept = w1(2) - lineSlope * w1(1);

p0 = [400, -200];   %initial Position of UAV (can change)
psi0 = pi/4;       %initial heading angle of the UAV (can change)

va = 25;            %UAV velocity (assuming constant throughout) (can change)

k1 = 100;
k2 = 1;
kappa = 0.5;

g = 9.81;
phiMax = 45*pi/180; %radians
psiDotMax = g*tan(phiMax)/va; %radians/sec

UAVState = zeros(simTime/dt, 5);
UAVState(1, 1:3) = [p0, psi0];

i = 2;
for t = dt:dt:simTime
    p = UAVState(i-1, 1:2); psi = UAVState(i-1, 3);
    x = p(1); y = p(2);

    thetad = atan2(w2(2) - y, w2(1) - x);
    % d = ((w1(2) - w2(2))*x + (w2(1) - w1(1))*y + w1(1)*w2(2) - w2(1)*w1(2))/sqrt((w2(1) - w1(1))^2 + (w2(2) - w1(2))^2);
    
    d = -(y - lineSlope*x - lineIntercept)/sqrt(1 + lineSlope^2);

    %PLOS Guidance Law
    dpsi = k1*wrapToPi(thetad - psi) + atan(k2*d);

    
    dpsi = max(min(dpsi, psiDotMax), -psiDotMax); %1st order Heading Control: psiDot = kappa*(psid - psi);
    psi = wrapToPi(psi + dpsi*dt);  %psi

    x = x + va*cos(psi)*dt;
    y = y + va*sin(psi)*dt;

    UAVState(i, :) = [x, y, psi, dpsi, d];

    i = i + 1;
end

figure;
hold on; grid on; grid minor;
plot([w1(1), w2(1)], [w1(2), w2(2)], "LineWidth", 1.5, "LineStyle", ":");
plot(p0(1), p0(2), "o", "MarkerSize", 8, "LineWidth", 1.5);
plot(UAVState(:, 1), UAVState(:, 2), "LineWidth", 1.5);
xlabel("x [m]"); ylabel("y [m]");
hold off;

figure;
subplot(2, 1, 1);
hold on; grid on; grid minor;
plot(0:dt:simTime, UAVState(:, 3)*180/pi, "LineWidth", 1.5);
xlabel("Time [s]"); ylabel("\psi [degrees]")
hold off;

subplot(2, 1, 2);
hold on; grid on; grid minor;
plot(0:dt:simTime, UAVState(:, 4), "LineWidth", 1.5);
xlabel("Time [s]"); ylabel("\psi_d_o_t [rad/s]")
hold off;

subplot(3, 1, 3);
hold on; grid on; grid minor;
plot(0:dt:simTime, UAVState(:, 5), "LineWidth", 1.5);
xlabel("Time [s]"); ylabel("Cross Track Error [m]")
hold off;

