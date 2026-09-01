clc; clear; close all;

simTime = 60;
dt = 0.001;

%%Path Specification
o = [0, 0]';   %m
r = 150;        %m
s = 1; %-1: clockwise; 1: anticlockwise

%%UAV Specifications
va = 25; %m/s

p0 = [-100, -300]';   %m
psi0 = -3*pi/4;       %radians

k1 = 1;
k2 = [0.1, 0.07, 0.04, 0.01];

g = 9.81;
phiMax = 45*pi/180; %radians
dpsiMax = g*tan(phiMax)/va; %radians/sec

Ne = size(k2, 2);
state = cell(Ne);
for e = 1:Ne
    state{e} = zeros([6, simTime/dt + 1]);
    state{e}(1:3, 1) = [p0', psi0]';
end

parfor e = 1:Ne
    i = 2;
    for t=dt:dt:simTime
        p = state{e}(1:2, i - 1); x = p(1); y = p(2);
        psi = state{e}(3, i - 1);

        d = sqrt((x - o(1))^2 + (y - o(2))^2) - r;
        theta = atan2(y - o(2), x - o(1));
        thetad = theta + s*(pi/2);

        dpsi = k1*wrapToPi(thetad - psi) + s*atan(k2(e)*d);

        psi = wrapToPi(psi + dpsi*dt);
        x = x + va*cos(psi)*dt;
        y = y + va*sin(psi)*dt;


        state{e}(:, i) = [x, y, psi, dpsi, wrapToPi(thetad - psi), d]';

        i = i + 1;
    end
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
h1 = plot(x, y,"--", lineWidth=1, Color="#ff0000");
h2 = plot(p0(1), p0(2), "Marker","o", "MarkerSize", 8, 'LineWidth', 1.5);

clr = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
leg = ["Loiter Circle",  "Initial Position"];

for e = 1:Ne
    leg(e+2) = "k_2 = " + k2(e);
    h(e) = plot(state{e}(1,:), state{e}(2, :), "Color", clr(e), 'LineWidth', 1.5);
    plot(state{e}(1, end), state{e}(2, end), '^', "MarkerSize", 8,"Color", clr(e), "LineWidth", 1.5);
end

xlabel("x [m]");
ylabel("y [m]");
title("PLOS Algorithm for loiter (Variation k_2)");

legend([h1, h2, h], leg);
hold off;