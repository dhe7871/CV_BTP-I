%this code does not deal with the vertical line following edge case
clear;
clc;
close all;

simTime = 40;
dt = 0.001;

w1 = [-200, -200];
w2 = [150, 400];

lineSlope = (w2(2) - w1(2))/(w2(1) - w1(1));
lineIntercept = w1(2) - lineSlope * w1(1);

p0 = [200, -100];   %initial Position of UAV (can change)
psi0 = [0, pi/2, pi, -pi/2];        %initial heading angle of the UAV (can change)

L = 100;            %finding radius L (can change)
va = 25;            %UAV velocity (assuming constant throughout) (cna change)
back = false;        %Tracing the line backward and forward (can change)

g = 9.81;
phiMax = 45*pi/180; %radians
psiDotMax = g*tan(phiMax)/va; %radians/sec


% lineParams = [m, c]
lineParams = [lineSlope, lineIntercept];


Ne = size(psi0, 2);
UAVState{1} = zeros(simTime/dt, 5);
UAVState{2} = zeros(simTime/dt, 5);
UAVState{3} = zeros(simTime/dt, 5);
UAVState{4} = zeros(simTime/dt, 5);

UAVState{1}(1, 1:3) = [p0, psi0(1)];
UAVState{2}(1, 1:3) = [p0, psi0(2)];
UAVState{3}(1, 1:3) = [p0, psi0(3)];
UAVState{4}(1, 1:3) = [p0, psi0(4)];



parfor e = 1:Ne
    i = 2;
    for t = dt:dt:simTime
        % circleParams = [x0, y0, r]
        circleParams = [UAVState{e}(i-1,1), UAVState{e}(i-1,2), L];
        [ptsExists, VTP] = getIntersectionPts(circleParams, lineParams, back);
        if ptsExists
            psid = atan2(VTP(2) - UAVState{e}(i-1,2), VTP(1) - UAVState{e}(i-1,1));
        else
            psid = atan(-1/lineSlope);
            if  UAVState{e}(i-1, 2) - lineSlope*UAVState{e}(i-1,1) - lineIntercept < 0
                if psid < 0
                    psid = pi + psid;
                end
            else
                if psid >= 0
                    psid = -pi + psid;
                end
            end
        end

        UAVState{e}(i-1, 5) = psid - UAVState{e}(i-1, 3);
        UAVState{e}(i-1, 4) = max(min(2*va*sin(UAVState{e}(i-1, 5))/L, psiDotMax), -psiDotMax); %Non-Linear guidance law: psiDot = 2*va*sin(psid - psi)/L
        UAVState{e}(i, 3) = wrapToPi(UAVState{e}(i-1, 3) + UAVState{e}(i-1, 4)*dt);  %psi

        %Using adams-moultan method for solving ODE
        UAVState{e}(i,1) = UAVState{e}(i-1, 1) + (va*cos(UAVState{e}(i-1, 3)) + va*cos(UAVState{e}(i, 3)))*dt/2;
        UAVState{e}(i,2) = UAVState{e}(i-1, 2) + (va*sin(UAVState{e}(i-1, 3)) + va*sin(UAVState{e}(i, 3)))*dt/2;

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
leg = ["Line to Follow",  "Initial Position", "\psi_0 = 0", "\psi_0 = 90^o", "\psi_0 = 180^o", '\psi_0 = 270^o'];

for e = 1:Ne
    p(e) = plot(UAVState{e}(:,1), UAVState{e}(:,2), "Color", clr(e), 'LineWidth', 1.5);
    plot(UAVState{e}(end,1), UAVState{e}(end,2), '^', "MarkerSize", 8,"Color", clr(e), "LineWidth", 1.5);
end

xlabel("x [m]");
ylabel("y [m]");
title("NLGL algorithm for line following (variation \psi)");

legend([h1, h2, p], leg);
hold off;


function [ptsExists, sol] = getIntersectionPts(circleParams, lineParams, back)
    a = 1 + lineParams(1)^2;
    b = 2*lineParams(1)*(lineParams(2) - circleParams(2)) - 2*circleParams(1);
    c = circleParams(1)^2 + (lineParams(2) - circleParams(2))^2 - circleParams(3)^2;

    D = b^2 - 4*a*c;

    if D < 0
        sol = [0, inf];
        ptsExists = false;
    else
        sgn = -1*(back) + 1*(~back); 
        x = (-b + sgn*sqrt(D))/(2*a);
        y = lineParams(1)*x + lineParams(2);
        sol = [x, y];
        ptsExists = true;
    end
end
