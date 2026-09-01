clc; clear; close all;

simTime = 60;
dt = 0.001;

%%Path Specification
o = [0, 0]';   %m
r = 150;        %m
direction = -1; %1: clockwise; -1: anticlockwise

%%UAV Specifications
va = 25; %m/s
L = 50; %m (should always be greater than r; if L > r, then the algorithm can be dynamically unstable i.e. oscillation never dies)

p0 = [-300, -200]';   %m
psi0 = [pi/4, 3*pi/4, -3*pi/4, -pi/4];        %radians

g = 9.81;
phiMax = 45*pi/180; %radians
dpsiMax = g*tan(phiMax)/va; %radians/sec


Ne = size(psi0, 2);
state = cell(Ne);
for e = 1:Ne
    state{e} = zeros([4, simTime/dt + 1]);
    state{e}(1:3, 1) = [p0', psi0(e)]';
end


parfor e = 1:Ne
    i = 2;
    for t=dt:dt:simTime
        p = state{e}(1:2, i - 1); x = p(1); y = p(2);
        psi = state{e}(3, i - 1);
        
        if(i == 2)
            cte = (x - o(1))^2 + (y - o(2))^2 - r^2;
            state{e}(6, 1) = cte;
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
    
        state{e}(:, i) = [x, y, psi, dpsi, eta, cte]';
    
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
h2 = plot(p0(1), p0(2), "Marker","o", "MarkerSize", 8, 'LineWidth', 1.5, "Color", "#0099aa");

clr = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
leg = ["Loiter Circle",  "Initial Position", "\psi_0 = 45^o", "\psi_0 = 135^0", "\psi_0 = -135^o", "\psi_0 = -45^o"];

for e = 1:Ne
    h(e) = plot(state{e}(1, :), state{e}(2, :), "LineWidth", 1.5, "Color", clr(e));
    plot(state{e}(1, end), state{e}(2, end), '^', "LineWidth", 1.5, "MarkerSize", 8, "Color", clr(e));
end

xlabel("x [m]");
ylabel("y [m]");
title("Non-Linear Guidance Law for Loiter (Variation \psi)");
legend([h1, h2, h], leg);
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