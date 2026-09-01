%this code does not deal with the vertical line following edge case
clear;
clc;
close all;

simTime = 50;
dt = 0.001;

w1 = [-200, -200];
w2 = [350, 700];
figure;
hold on;
grid on; grid minor;
axis equal;
plot([w1(1), w2(1)], [w1(2), w2(2)],"--", lineWidth=1, Color="#ff0000");
xlabel("x-axis");
ylabel("y-axis");
title("Line following with NLGL");



lineSlope = (w2(2) - w1(2))/(w2(1) - w1(1));
lineIntercept = w1(2) - lineSlope * w1(1);

p0 = [150, -225];   %initial Position of UAV (can change)
psi = -pi/2;        %initial heading angle of the UAV (can change)

L = 100;            %finding radius L (can change)
va = 25;            %UAV velocity (assuming constant throughout) (cna change)
back = false;        %Tracing the line backward and forward (can change)

g = 9.81;
phiMax = 45*pi/180; %radians
psiDotMax = g*tan(phiMax)/va; %radians/sec


% lineParams = [m, c]
lineParams = [lineSlope, lineIntercept];

plot(p0(1), p0(2), "o", "LineWidth", 1.5, "MarkerSize", 8, "Color", "#00aa99", "DisplayName", "Initial Position");

UAVState = zeros([simTime/dt, 5]);
UAVState(1,:) = [p0, psi, 0, 0];
i = 2;
for t = dt:dt:simTime
    % circleParams = [x0, y0, r]
    circleParams = [UAVState(i-1,1), UAVState(i-1,2), L];
    [ptsExists, VTP] = getIntersectionPts(circleParams, lineParams, back);
    if ptsExists
        psid = atan2(VTP(2) - UAVState(i-1,2), VTP(1) - UAVState(i-1,1));
    else
        psid = atan(-1/lineSlope);
        if  UAVState(i-1, 2) - lineSlope*UAVState(i-1,1) - lineIntercept < 0
            if psid < 0
                psid = pi + psid;
            end
        else
            if psid >= 0
                psid = -pi + psid;
            end
        end
    end
    
    UAVState(i-1, 5) = psid - UAVState(i-1, 3);
    UAVState(i-1, 4) = max(min(2*va*sin(UAVState(i-1, 5))/L, psiDotMax), -psiDotMax); %Non-Linear guidance law: psiDot = 2*va*sin(psid - psi)/L
    UAVState(i, 3) = wrapToPi(UAVState(i-1, 3) + UAVState(i-1, 4)*dt);  %psi

    %Using adams-moultan method for solving ODE
    UAVState(i,1) = UAVState(i-1, 1) + (va*cos(UAVState(i-1, 3)) + va*cos(UAVState(i, 3)))*dt/2;
    UAVState(i,2) = UAVState(i-1, 2) + (va*sin(UAVState(i-1, 3)) + va*sin(UAVState(i, 3)))*dt/2;


    if ~mod(i, 100)
        if i ~= 100
            delete(h1);
            delete(h2);
        end
        h1 = plot(UAVState(1:i,1), UAVState(1:i,2), linewidth = 1.5, color='#0099aa');
        h2 = plot(UAVState(i,1), UAVState(i,2), '^', markersize = 8, color = '#0099ff', linewidth=1.5);
        
        pause(dt/1000);
    end

    i = i + 1;
end

legend("Target Line", "Actual traced line", "current position");
hold off;

figure;
subplot(1, 3, 1);
plot( 0:dt:simTime,UAVState(:,3)*180/pi, LineWidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("Heading Angle \psi [degrees]");
% title("Variation in UAV Heading Angle with time");
hold off;

subplot(1, 3, 2);
plot( 0:dt:simTime, UAVState(:,4), LineWidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time (in s)");
ylabel("Heading Angle rate \psi_d_o_t [rad/s]");
% title("Variation in UAV Heading Angle rate with time");
hold off;

subplot(1, 3, 3);
plot( 0:dt:simTime, UAVState(:,5)*180/pi, LineWidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("Heading Angle Error ,(\eta = \psi_d - \psi) [degrees]");
% title("Variation in UAV Heading Error with time");
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
