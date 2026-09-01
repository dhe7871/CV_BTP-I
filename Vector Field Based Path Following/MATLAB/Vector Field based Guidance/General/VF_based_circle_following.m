clear; clc; close all;

%Time definition
sim_time = 50;
dt = 0.001;

%circle definition
r = 150;
O = [10, 10];

theta = linspace(0, 2*pi, 500);


writerObj = VideoWriter('VFLoiterAnimation.avi');
open(writerObj);

figure;
hold on;
grid on; grid minor;
axis equal;
xlabel("x [m]");
ylabel("y [m]");
title("Loiter Path Following with VF based guidance");
plot(O(1) + r*cos(theta), O(2) + r*sin(theta),"--", lineWidth=1.5, color="#ff0022");

%curvature constraint definition
g = 9.81;
phi_max = pi/4;

%UAV definition
vg = 25;
p0 = [300, 0];
chi = 3*pi/4;
chi_inf = pi/2;

reverse = true;

% gains
alpha = 1.65; %i guess will be predefined

k = 0.05;
kappa = pi/2;

epsilon = 0.1;

state = zeros([sim_time/dt + 1, 7]);
state(1, :) = [p0,chi, 0, 0, 0, 0];

plot(state(1,1), state(1, 2), 'o', markersize=8, MarkerEdgeColor='#0088ff', LineWidth=1.5);

%curvature constraint definition
chi_dot_max = (g*tan(phi_max))/vg;

i = 2;
for t=dt:dt:sim_time
    %calculation for previous step
    chi = state(i-1, 3);
    d = norm(state(i-1, 1:2) - O);
    gamma = atan2(state(i-1, 2) - O(2), state(i-1, 1) - O(1));
    
    chi_d = wrapToPi(gamma - (reverse*(-1) + (~reverse))*pi/2 - (reverse*(-1) + (~reverse))*atan(k*(d-r)));
    chi_tilda = wrapToPi(chi - chi_d);

    beta = k/(1+(k*(d-r))^2);

    chi_c = wrapToPi(chi - (vg/(alpha*d))*sin(chi - gamma) - (beta/alpha)*vg*cos(chi - gamma) - (kappa/alpha)*sat(chi_tilda/epsilon));
    
    chi_dot = max(min(alpha*wrapToPi(chi_c - chi), chi_dot_max), -chi_dot_max);
    
    state(i-1, 4:7) = [chi_dot, chi_c, chi_d, beta];

    %calculation for this step
    chi = chi + chi_dot * dt;
    
    state(i, 1) = state(i-1, 1) + vg*cos(chi)*dt;
    state(i, 2) = state(i-1, 2) + vg*sin(chi)*dt;
    state(i, 3) = chi;

    if ~mod(i, 100)
        if i ~= 100
            delete(h1);
            delete(h2);

        end
      
        h1 = plot(state(1:i, 1), state(1:i, 2), lineWidth = 1.5, color = "#0088ff");
        h2 = plot(state(i,1), state(i,2), marker="^", markersize = 8, color = '#0088ff', linewidth=1.5);
        
        frame = getframe(gcf);
        writeVideo(writerObj, frame);

        pause(dt);
    end

    i = i + 1;
end

legend("Orbit to follow", "Initial Position", "Actual trajectory", "Current Position");
hold off;

frame = getframe(gcf);
writeVideo(writerObj, frame);
close(writerObj);

figure;
subplot(3, 1, 1);
hold on;
plot(0:dt:(sim_time-dt), state(1:end-1, 6)*180/pi,"--", lineWidth=1.5);
plot(0:dt:(sim_time-dt), state(1:end-1, 5)*180/pi,lineStyle=":", lineWidth=1.5);
plot(0:dt:sim_time, wrapToPi(state(:, 3))*180/pi, lineWidth=1.5);
grid on; grid minor;
xlabel("Time [s]");
ylabel("\chi [deg]");
title("Course Angle (\chi) Variation");
legend("Desired, \chi^d", "Commanded, \chi^c", "Actual, \chi");
hold off;

subplot(3, 1, 2);
hold on;
grid on; grid minor;
plot(0:dt:(sim_time-dt), state(1:end-1, 4), lineWidth=1.5);
title("Course rate variation with time");
xlabel("Time [s]");
ylabel("\chi_d_o_t [rad/s]");
hold off;

subplot(3, 1, 3);
hold on;
grid on; grid minor;
plot(0:dt:(sim_time-dt), state(1:end-1, 7), lineWidth=1.5);
% title("Cross track error variation with time");
xlabel("Time [s]");
ylabel("\beta [m]");
hold off;


function val = sat(x)
    if abs(x) <= 1
        val = x;
    else
        val = sign(x);
    end
end
