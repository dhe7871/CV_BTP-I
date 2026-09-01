clc; clear; close all;

simTime = 100; %s
dt= 0.001; %s

syms x y;
a0 = 200;
b0 = 10;
% fn = y - 3*x;
% fn = (x - 2)^2 + (y + 10)^2 - 9775;
% fn = (y - 200*sin(x/50))/100;
% fn =(1/((x^2 + y^2)^(9/2) + 1))*...
%     (-a0*x^8 - 8*b0*x^7*y - 4*a0*x^6*y^2 + 56*b0*x^5*y^3 - 6*a0*x^4*y^4 - 56*b0*x^3*y^5 - 4*a0*x^2*y^6 + 8*b0*x*y^7 - a0*y^8 ...
%     + (x^8 + x^6*y^2 + x^4*y^4 + x^2*y^6 + y^8)*sqrt(x^2 + y^2));
% fn = (1/((5*x^2 + x*y + 5*y^2)^4.5 + 1))*(-(5*x^2 + x*y + 5*y^2)^4.5 ...
%     + 300*sqrt(2)*(5*x^2 + x*y + 5*y^2)^4 ...
%     + 45*sqrt(22)*(x^2 - y^2)*(x^2 + 20*x*y + y^2)*(49*x^4 - 20*x^3*y - 300*x^2*y^2 - 20*x*y^3 + 49*y^4));
S = 100;

fn = (x./S).^2 + (y./S).^2 ...
   + 0.6*exp(-((x./S - 0.5).^2 + (y./S + 0.3).^2)) ...
   + 0.4*exp(-((x./S + 0.7).^2 + (y./S - 0.6).^2)) ...
   + 0.35*exp(-((x./S + 0.2).^2 + (y./S + 0.8).^2)) ...
   + 0.25*exp(-((x./S - 0.9).^2 + (y./S - 0.1).^2)) ...
   + 0.15*sin(3*x./S).*cos(2*y./S) ...
   + 0.12*sin(5*y./S + 0.7) ...
   + 0.10*cos(4*x./S - 1.1) ...
   + 0.08*(x./S).*(y./S) - 1.8;


%differentiated matlab functions
fm = matlabFuncTillDoubleDiff(fn);

%UAV Specifications
s = 1; %direction
va = 15; %m/s
w = [6, 8]; %m/s

p0 = [-200, -200]; %m
chi0 = -3*pi/4; %rad

g = 9.81;
phiMax = pi/4;
chiDotMax = g*tan(phiMax)/va; %rads/s

kappa0 = 10;
kChi = 1.5;

%State Initialization
state = zeros([simTime/dt + 1, 10]);
state(1, :) = [p0, zeros([1, 2]), 0, chi0, 0, zeros([1, 3])];

writerObj = VideoWriter('CVF2DPathFollowingAnimation.avi');
open(writerObj);

figure;
hold on;
grid on; grid minor;
axis equal;
% xlim([-400, 200]); ylim([-400, 200]);
fimplicit(fn, [-300 300 -300 200], "LineStyle", "--","LineWidth",1.5, "Color","#ff0055");
scatter(p0(1), p0(2), "LineWidth", 1.5, "Marker", "o", "SizeData", 100);
h1 = plot(NaN, NaN, "LineWidth", 1.5, "Color", "#0099aa");
h2 = plot(NaN, NaN, "Marker","^", "MarkerSize", 8, "Color",'#0099aa', "LineWidth", 1.5);
legend("Curve to follow", "Initial Position", "Actual trajectory", "Current Position");
title("Path Following with CVF based guidance");


i = 2;
for t=dt:dt:simTime
    df = calculateDerivatives(fm, state(i-1, 1:2));

    chi = state(i-1, 6);

    vg = (w(1)*cos(chi) + w(2)*sin(chi)) + sqrt(va^2 - (w(1)*sin(chi) - w(2)*cos(chi))^2);
    r = df(1);
    state(i-1, 10) = r;
    if i==2    
        kappa = kappa0;
        if r > 0.1
            kappa = kappa0/abs(r);
        end
    end

    vc = [df(2);df(3)];
    vc = vc./sqrt(vc'*vc);
    vs = [df(3);-df(2)];
    vs = vs./sqrt(vs'*vs);
    vd = -vg*tanh(kappa*r)*vc + s*vg*sech(kappa*r)*vs;

    chid = atan2(vd(2), vd(1));
    chie = wrapToPi(chi - chid);

    %coeffiecient functions vector of time derivative of desired course angle
    % chi_d_dot = A_x*x_dot + A_y*y_dot
    Ax = -s*kappa*df(2)*sech(kappa*df(1)) - (df(4)*df(3) - df(2)*df(6))/(df(2)^2 + df(3)^2);
    Ay = -s*kappa*df(3)*sech(kappa*df(1)) + (df(5)*df(2) - df(3)*df(6))/(df(2)^2 + df(3)^2); 

    chidDot = vg*(Ax*cos(chi) + Ay*sin(chi));

    omega1 = -kChi*chie + chidDot;
    chicDot = max(min(omega1, chiDotMax), -chiDotMax);

    state(i-1, 5) = vg;
    state(i-1, 7:9) = [chid, chidDot,chicDot];

    chi = wrapToPi(chi + chicDot*dt);

    xDot = vg*cos(chi);
    yDot = vg*sin(chi);

    state(i, 1:2) = state(i-1, 1:2) + dt*[xDot, yDot];
    state(i, 3:4) = [xDot, yDot];
    state(i, 6) = chi;

    if ~mod(i, 100)
        set(h1, "XData", state(1:i, 1), "YData",state(1:i, 2));
        set(h2, "XData", state(i, 1), "YData", state(i, 2));

        drawnow limitrate;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end


    i = i + 1;
end

hold off;

frame = getframe(gcf);
writeVideo(writerObj, frame);
close(writerObj);

figure;
subplot(2, 2, 1);
plot(0:dt:simTime, state(:, 7)*180/pi,"--", linewidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 6)*180/pi, lineWidth=1.5);
xlabel("Time [s]");
ylabel("\chi, [deg]");
legend("Desired","Actual");
title("Course angle variation with time");
hold off;

subplot(2, 2, 2);
plot(0:dt:simTime, state(:, 10), lineWidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("Cross Track Error");
title("Cross Track Error variation with time");
hold off;

subplot(2, 2, 3);
plot(0:dt:simTime, state(:, 5), linewidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("V_g, [m/s]");
title("Ground speed variation with time");
hold off;

subplot(2, 2, 4);
plot(0:dt:simTime, state(:, 8),"--", linewidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 9), linewidth=1.5);
xlabel("Time [s]");
ylabel("d\chi/dt, [rad/s]");
legend("Desired", "commanded");
title("Course rate variation with time");
hold off;


function fnm = matlabFuncTillDoubleDiff(fn)
    syms x y;
    fx = diff(fn, x);
    fy = diff(fn , y);
    fxx = diff(fx, x);
    fyy = diff(fy, y);
    fxy = diff(fx, y);

    fm = matlabFunction(fn, "Vars", [x, y]);
    fxm = matlabFunction(fx, "Vars", [x, y]);
    fym = matlabFunction(fy, "Vars", [x, y]);
    fxxm = matlabFunction(fxx, "Vars", [x, y]);
    fyym = matlabFunction(fyy, "Vars", [x, y]);
    fxym = matlabFunction(fxy, "Vars", [x, y]);

    fnm = {fm, fxm, fym, fxxm, fyym, fxym};
end

function df = calculateDerivatives(fm, p)
    a = p(1);
    b = p(2);

    f = fm{1};
    fx = fm{2};
    fy = fm{3};
    fxx = fm{4};
    fyy = fm{5};
    fxy = fm{6};

    df = [f(a, b), fx(a, b), fy(a, b), fxx(a, b), fyy(a, b), fxy(a, b)];
end
