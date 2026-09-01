clear; clc; close all;

%Time definition
sim_time = 30;
dt = 0.001;


%Path definition
length = 500;
d0 = [100, 200];
dv = [0.5, 1];
df = d0 + length * dv;
dv = dv/norm(dv);

writerObj = VideoWriter('VFLineAnimation.avi');
open(writerObj);

figure(1);
hold on;
grid on; grid minor;
axis equal;
xlabel("x [m]");
ylabel("y [m]");
title("Line Path Following with VF based guidance");
plot([d0(1), df(1)], [d0(2), df(2)],":", lineWidth=1.5, color="#ff0022");

%curvature constraint definition
g = 9.81;
phi_max = pi/4;


%UAV definition
vg = 25;
p0 = [-200, 250];
chi = 3*pi/4;
chi_inf = pi/2 + atan2(dv(2), dv(1));

% gains
alpha = 1.65; %i guess will be predefined

k = 0.02;
kappa = pi/2;

epsilon = 0.1;


state = zeros([sim_time/dt + 1, 7]);
state(1, :) = [p0, chi, 0, 0, 0, 0];


plot(state(1,1), state(1, 2), 'o', markersize=8, MarkerEdgeColor='#0099aa', LineWidth=1.5);


%curvature constraint definition
chi_dot_max = (g*tan(phi_max))/vg;
i = 2;
for t=dt:dt:sim_time
    %calculation for previous step
    chi = state(i-1, 3);
    dist = cross([state(i-1, 1:2) - d0, 0], [dv, 0]);
    y = -dist(3);

    chi_d = wrapToPi(atan2(dv(2),dv(1)) - chi_inf*(2/pi)*tanh(k*y));
    chi_tilda = wrapToPi(chi - chi_d);
    chi_c = wrapToPi(chi - (1/alpha)*chi_inf*(2/pi)*(k/(1+ (k*y)^2))*vg*sin(chi) - (kappa/alpha)*sat(chi_tilda/epsilon));
    
    chi_dot = max(min(alpha*wrapToPi(chi_c - chi), chi_dot_max), -chi_dot_max);
    
    state(i-1, 4:7) = [chi_dot, chi_c, chi_d, y];

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
      
        h1 = plot(state(1:i, 1), state(1:i, 2), lineWidth = 1.5, color = "#0099aa");
        h2 = plot(state(i,1), state(i,2), marker="^", markersize = 8, color = '#0099aa', linewidth=1.5);
        
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
        pause(dt);
    end

    i = i + 1;
end

legend("Line to follow", "Initial Position", "Actual trajectory", "Current Position");
hold off;

frame = getframe(gcf);
writeVideo(writerObj, frame);
close(writerObj);


figure(2);
subplot(3, 1, 1);
hold on;
plot(0:dt:(sim_time-dt), state(1:end-1, 6)*180/pi,"--", lineWidth=1.5);
plot(0:dt:(sim_time-dt), state(1:end-1, 5)*180/pi,lineStyle=":", lineWidth=1.5);
plot(0:dt:sim_time, state(:, 3)*180/pi, lineWidth=1.5);
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
title("Cross track error variation with time");
xlabel("Time [s]");
ylabel("d [m]");
hold off;



function val = sat(x)
    if abs(x) <= 1
        val = x;
    else
        val = sign(x);
    end
end
