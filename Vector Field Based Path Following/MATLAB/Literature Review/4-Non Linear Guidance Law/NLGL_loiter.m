clc; clear; close all;

simTime = 50;
dt = 0.001;

%%Path Specification
o = [0, 0]';   %m
r = 150;        %m
direction = -1; %1: clockwise; -1: anticlockwise

%%UAV Specifications
va = 25; %m/s
L = 50; %m (should always be greater than r; if L > r, then the algorithm can be dynamically unstable i.e. oscillation never dies)

p0 = [-350, -555]';   %m
psi0 = -3*pi/4;        %radians

g = 9.81;
phiMax = 45*pi/180; %radians
dpsiMax = g*tan(phiMax)/va; %radians/sec

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
xlabel("x-axis");
ylabel("y-axis");
title("Loitering with NLGL");

state = zeros([4, simTime/dt + 1]);
state(1:3, 1) = [p0', psi0]';


i = 2;
for t=dt:dt:simTime
    p = state(1:2, i - 1); x = p(1); y = p(2);
    psi = state(3, i - 1);
    
    if(i == 2)
        cte = (x - o(1))^2 + (y - o(2))^2 - r^2;
        state(6, 1) = cte;
    end

    s = getIntersectionPt(o, r, p, L, direction);
    
    normOP = (p - o)/sqrt((p - o)'*(p - o));
    if(isnan(s))
        s = o + r*normOP;
    end

    psit = atan2(s(2) - p(2), s(1) - p(1));
    eta = wrapToPi(psit - psi);

    dpsi = max(min(2*va*sin(eta)/L, dpsiMax), -dpsiMax);
    psi = wrapToPi(psi + dpsi*dt);
    x = x + va*cos(psi)*dt;
    y = y + va*sin(psi)*dt;

    cte = (x - o(1))^2 + (y - o(2))^2 - r^2;

    state(:, i) = [x, y, psi, dpsi, eta, cte]';

    if ~mod(i, 100)
        if i ~= 100
            delete(h1);
            delete(h2);
            delete(h3);
            delete(h4);
            delete(h5);
            delete(h6);
        end
        h1 = plot(state(1, 1:i), state(2, 1:i), linewidth = 1.5, color='#0099aa');
        h2 = plot(state(1, i), state(2, i), '^', markersize = 8, color = '#0099ff', linewidth=1.5);
        
        tx = (x - L):0.01:(x + L);
        ty1 = y - sqrt(L^2 - (tx - x).^2 + 0.01);
        ty2 = y + sqrt(L^2 - (tx - x).^2 + 0.01);
        tx = [tx(1:end-1), flip(tx)]';
        ty = [ty1(1:end-1), flip(ty2)]';
        h3 = plot(tx, ty,":", lineWidth=1, Color="#00ff22");

        h4 = plot(s(1), s(2), "Marker", "square", "MarkerSize", 8, "Color","#ffff00");
        h5 = plot([x, s(1)], [y, s(2)], "LineStyle","-.", "Color", "#ffff00");

        w = [x, y]' + [cos(psi), sin(psi)]'.*L;
        h6 = plot([x, w(1)], [y, w(2)], "Color", "#00ff55");

        pause(dt);
    end

    i = i + 1;
end


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

%%function which gives the one of the two intersection pts according to the
%%direction specified,
%%If the circles do not intersect in that case "NaN" is sent instead of a
%%intersection pt.
function s = getIntersectionPt(o, r, p, L, direction)
    d = sqrt((p - o)'*(p - o));

    if(abs(d - r - L) < 1e-2)
        s = o + (p - o).*(r/d);
    elseif (abs(d + r - L) < 1e-2)
        s = o - (p - o).*(r/d);
    elseif(d > r + L || L > d + r || r > d + L)
        %1- no intersection of between the circles (both circles are distant)
        %2- one circle is inside other circle
        s = NaN;
    else
        %%all small variables are "scalar" and all caps are "vector" and norm__
        %%are nothing but the "normalized vectors"
        OP = (p - o);
        op = sqrt((p - o)'*(p - o));
        normOP = OP./op;

        a = (r^2 - L^2 + op^2)/(2*op);
        OR = a.*normOP;

        normRQ = [-normOP(2), normOP(1)]';
        normRS = -normRQ;
        rq = sqrt(r^2 - a^2);
        rs = rq;

        RQ = rq*normRQ;
        RS = rs*normRS;

        R = o + OR;

        Q = R + RQ;
        S = R + RS;
        pts = [Q, S];

        PO = -OP;
        PSSdash = pts - p;
        sgn = sign([ ...
            sum(cross([PO', 0]', [PSSdash(:, 1)', 0]'), "all"), ...
            sum(cross([PO', 0]', [PSSdash(:, 2)', 0]'), "all") ...
            ]');
        mask = direction == sgn;
        s = pts*mask;
    end
end