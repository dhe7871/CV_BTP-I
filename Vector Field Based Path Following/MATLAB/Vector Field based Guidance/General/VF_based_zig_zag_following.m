clear; clc; close all;

%Time definition
sim_time = 100;
dt = 0.001;


%Path definition
d0 = [0, 500];
d1 = [150, 500];
d2 = [300, 300];
d3 = [500, 450];
d4 = [500, 0];
d5 = [200, -100];
d6 = [-100, 0];
d7 = [0, 500];

lines = [d0,(d1-d0)/norm(d1-d0), d1;
      d1,(d2-d1)/norm(d2-d1), d2;
      d2,(d3-d2)/norm(d3-d2), d3;
      d3,(d4-d3)/norm(d4-d3), d4;
      d4,(d5-d4)/norm(d5-d4), d5;
      d5,(d6-d5)/norm(d6-d5), d6;
      d6,(d7-d6)/norm(d7-d6), d7;
      ];

writerObj = VideoWriter('VFLineAnimation.avi');
open(writerObj);

figure;
hold on;
grid on; grid minor;
axis equal;
xlabel("x [m]");
ylabel("y [m]");
title("Complex Path (Quadrilateral) Path Following with VF based guidance");

for line=lines'
    plot([line(1), line(5)], [line(2), line(6)],"--", lineWidth=1.5, color="#ff0022");
end

%curvature constraint definition
g = 9.81;
phi_max = pi/4;


%UAV definition
vg = 25;
p0 = [-200, 300];
chi = -3*pi/4;
chi_inf = pi/2;


% gains
alpha = 1.65; %i guess will be predefined

k = 0.03;
kappa = pi/2;

epsilon = 0.1;


state = zeros([sim_time/dt + 1, 7]);
state(1, :) = [p0,chi, 0, 0, 0, 0];

plot(state(1,1), state(1, 2), 'o', markersize=8, MarkerEdgeColor='#0099aa', LineWidth=1.5);


%curvature constraint definition
chi_dot_max = (g*tan(phi_max))/vg;

Rmax = vg/chi_dot_max;

%first line selection
min_dist = inf;
c_line = 1;
i = 1;
for line=lines'
    dist = norm(cross([p0-line(1:2)', 0], [line(3:4)', 0]));
    if  dist < min_dist
        min_dist = dist;
        c_line = i;
    end
    i = i + 1;
end


i = 2;
for t=dt:dt:sim_time
    %calculation for previous step
    chi = state(i-1, 3);
    dist = cross([state(i-1, 1:2) - lines(c_line,1:2), 0], [lines(c_line, 3:4), 0]);
    y = -dist(3);

    chi_d = wrapToPi(atan2(lines(c_line, 4), lines(c_line, 3)) - chi_inf*(2/pi)*atan(k*y));
    chi_tilda = wrapToPi(chi - chi_d);
    chi_c = wrapToPi(chi - (1/alpha)*chi_inf*(2/pi)*(k/(1+ (k*y)^2))*vg*sin(chi) - (kappa/alpha)*sat(chi_tilda/epsilon));
    
    chi_dot = max(min(alpha*wrapToPi(chi_c - chi), chi_dot_max), -chi_dot_max);
    
    state(i-1, 4:7) = [chi_dot, chi_c, chi_d, y];

    %calculation for this step
    chi = chi + chi_dot * dt;
    
    state(i, 1) = state(i-1, 1) + vg*cos(chi)*dt;
    state(i, 2) = state(i-1, 2) + vg*sin(chi)*dt;
    state(i, 3) = chi;

    ltheta = pi - acos(lines(c_line, 3:4)*lines((c_line+1)*(size(lines, 1) >= c_line+1)...
        + 1*(size(lines, 1) < c_line+1), 3:4)');

    if norm(lines(c_line, 5:6) - state(i, 1:2)) <= Rmax/tan(ltheta/2)
        c_line = c_line + 1;

        disp(Rmax); disp(Rmax/tan(ltheta/2))

        if c_line > size(lines, 1)
            c_line = 1;
        end
    end
        
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

hold off;

close(writerObj);

figure;
subplot(3, 1, 1);
hold on;
plot(0:dt:(sim_time-dt), state(1:end-1, 6)*180/pi,"--", lineWidth=1.5);
plot(0:dt:(sim_time-dt), state(1:end-1, 5)*180/pi,lineStyle=":", lineWidth=1.5);
plot(0:dt:sim_time, state(:, 3)*180/pi, lineWidth=1.5);
grid on; grid minor;
xlabel("Time [s]");
ylabel("\chi [deg]");
title("Course Angle (\chi) Variation with time");
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

