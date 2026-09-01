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

p0 = [0, -200];   %initial Position of UAV (can change)
psi0 = -3*pi/4;        %initial heading angle of the UAV (can change)

va = 25;            %UAV velocity (assuming constant throughout) (can change)

k1 = [1, 1.5, 2, 2.5];
k2 = 0.1;
kappa = 0.5;

g = 9.81;
phiMax = 45*pi/180; %radians
psiDotMax = g*tan(phiMax)/va; %radians/sec



Ne = size(k1, 2);
UAVState = cell(Ne);
for e = 1:Ne
    UAVState{e} = zeros(simTime/dt, 5);
    UAVState{e}(1, 1:3) = [p0, psi0];
end




parfor e = 1:Ne
    i = 2;
    for t = dt:dt:simTime
        p = UAVState{e}(i-1, 1:2); psi = UAVState{e}(i-1, 3);
        x = p(1); y = p(2);

        thetad = atan2(w2(2) - y, w2(1) - x);
        % d = ((w1(2) - w2(2))*x + (w2(1) - w1(1))*y + w1(1)*w2(2) - w2(1)*w1(2))/sqrt((w2(1) - w1(1))^2 + (w2(2) - w1(2))^2);

        d = -(y - lineSlope*x - lineIntercept)/sqrt(1 + lineSlope^2);

        %PLOS Guidance Law
        dpsi = k1(e)*wrapToPi(thetad - psi) + atan(k2*d);

        dpsi = max(min(dpsi, psiDotMax), -psiDotMax); %1st order Heading Control: psiDot = kappa*(psid - psi);
        psi = wrapToPi(psi + dpsi*dt);  %psi

        x = x + va*cos(psi)*dt;
        y = y + va*sin(psi)*dt;

        UAVState{e}(i, :) = [x, y, psi, dpsi, d];

        i = i + 1;
    end
end


figure;
hold on;
grid on; grid minor;
axis equal;
h1 = plot([w1(1), w2(1)], [w1(2), w2(2)],"--", "LineWidth", 1.5, "Color", "#ff0000" );

h2 = plot(p0(1), p0(2), "o", "LineWidth", 1.5, "MarkerSize", 8, "Color", "#00aa99");
clr = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
leg = ["Line to Follow",  "Initial Position"];

for e = 1:Ne
    leg(e+2) = "k_1 = " + k1(e);
    h(e) = plot(UAVState{e}(:,1), UAVState{e}(:,2), "Color", clr(e), "LineWidth", 1.5);
    plot(UAVState{e}(end,1), UAVState{e}(end,2), '^', "MarkerSize", 8,"Color", clr(e), "LineWidth", 1.5);
end

xlabel("x [m]");
ylabel("y [m]");
title("PLOS algorithm for line following (variation k_1)");

legend([h1, h2, h], leg);
hold off;