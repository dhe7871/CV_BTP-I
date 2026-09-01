clc; clear; close all;
addpath(fullfile(pwd, "v1/lib"));

try
    load(fullfile(pwd, "curve.mat"));
catch exception
    disp("Curve.mat doesn't exist...");
end

simTime = 100; %s
dt= 0.001; %s

syms x y z st vgt kappat;
%true: for using previously generated curve
prevCurve = false;
% fn = ((x - 2)^2 + (y + 10)^2 + (z - 15)^2 - 10000)/100;   %kappa0 = 10;
% gn = 2*x + 3*y - 2*z - 10;
% fn = x - 200*cos(z/20);     %kappa0 = 8;
% gn = y - 150*sin(z/20);
% fn = 100*(x + y + z)/(sqrt(5*x^2 + x*y + 5*y^2) + 1);       %kappa0 = 5;
% gn = 100*(1/((5*x^2 + x*y + 5*y^2)^4.5 + 1))*(-(5*x^2 + x*y + 5*y^2)^4.5 ...
%     + 300*sqrt(2)*(5*x^2 + x*y + 5*y^2)^4 ...
%     + 45*sqrt(22)*(x^2 - y^2)*(x^2 + 20*x*y + y^2)*(49*x^4 - 20*x^3*y - 300*x^2*y^2 - 20*x*y^3 + 49*y^4));

fn = (x/100)^2 + (y/100)^2 + (z/50)^2 - 1;
gn = (x + y + 3*z - 1)/100;

%differentiated matlab functions
fm = matlabFuncTillDoubleDiff(fn);
gm = matlabFuncTillDoubleDiff(gn);

%UAV Specifications
s = 1; %direction
va = 20; %m/s
w = {@(t) 6, @(t) 8, @(t) 0}; %m/s

p0 = [200, 200, -100]; %m
chi0 = pi/3; %rad
gamma0 = pi/4; %rad

gammaMax = pi/3;
chiDotMax = 0.5; %rads/s
gammaDotMax = 0.3; %rad/s

kappa0 = 8;
kChi = 1;
kGamma = 1;

%coefficient functions vector of time derivative of desired course and climb angle
% chi_d_dot = B_x*x_dot + B_y*y_dot + B_Z*z_dot
% gamma_d_dot = C_x*x_dot + C_y*y_dot + C_z*z_dot
coeffs = getGradCoeffsMatlabFunc(fn, gn);
BC = cell(6);
parfor i = 1:6
    row = ceil(i/3);
    col = i + 3 - 3*row;
    BC{i} = matlabFunction(coeffs(row, col), "Vars", [x, y, z, st, vgt, kappat]);
end
Bx = BC{1}; By = BC{2}; Bz = BC{3};
Cx = BC{4}; Cy = BC{5}; Cz = BC{6};


state = zeros([simTime/dt + 1, 16]);
state(1, :) = [p0, zeros([1, 3]), 0, chi0, 0, gamma0, 0, zeros([1, 5])];

writerObj = VideoWriter('VF3DHelix.avi');
open(writerObj);

fig = figure();
ax = axes(fig);
hold(ax, "on");
grid(ax, "on"); grid(ax, "minor");
view(ax, 45, 25);
axis(ax, "equal");
camproj(ax, "perspective");
xlabel(ax, "x [m]"); ylabel(ax, "y [m]"); zlabel(ax, "z [m]");

if(~prevCurve || ~exist("curve", "var"))
    curve = generate3DCurvePoints(fm{1}, gm{1}, p0, 50000);
    save("curve.mat", "curve");
end
plot3(curve(:, 1), curve(:, 2), curve(:, 3), "LineWidth",1.5, "Color","#ff5500", "LineStyle","--");
iniPt = plot3(ax, state(1, 1), state(1, 2), state(1,3), "LineWidth", 1.5, "MarkerEdgeColor","#0099aa","Marker", "o", "MarkerSize",8);
line = plot3(ax, NaN, NaN, NaN, "LineWidth",1.5, "Color","#0099aa");
pt = plot3(ax, NaN, NaN, NaN, "^", "Color","#0099aa", "MarkerSize",8, "LineWidth",1.5);
title(ax, "Inclined Circle Path following with CVF based guidance");

clear textprogressbar;
textprogressbar('Progress: ');
i = 2;
for t=dt:dt:simTime
    w1 = w{1}; w2 = w{2};

    df = calculateDerivatives(fm, state(i-1, 1:3));
    dg = calculateDerivatives(gm, state(i-1, 1:3));

    chi = state(i-1, 8);
    gamma = state(i-1, 10);

    vg = (w1(t)*cos(chi) + w2(t)*sin(chi))*cos(gamma) ...
        + 0.5*sqrt(4*va^2 - 3*(w1(t)^2 + w2(t)^2) + (1 + 2*cos(2*gamma))*(w1(t)*cos(chi) + w2(t)*sin(chi))^2 - (w1(t)*sin(chi) - w2(t)*cos(chi))^2);
    r = sqrt(df(1)^2 + dg(1)^2);
    state(i-1, 16) = r;

    if i == 2
        kappa = kappa0;
        if r > 0.001
            kappa = kappa0/abs(r);
        end
    end


    vc = [df(1)*df(2) + dg(1)*dg(2); df(1)*df(3) + dg(1)*dg(3); df(1)*df(4) + dg(1)*dg(4);];
    vc = vc./sqrt(vc'*vc);
    vs = [df(3)*dg(4) - df(4)*dg(3); df(4)*dg(2) - df(2)*dg(4); df(2)*dg(3) - df(3)*dg(2);];
    vs = vs./sqrt(vs'*vs);
    vd = -vg*tanh(kappa*r)*vc + s*vg*sech(kappa*r)*vs;

    chid = atan2(vd(2), vd(1));
    gammad = atan2(vd(3), sqrt(vd(1)^2 + vd(2)^2));

    chie = wrapToPi(chi - chid);
    gammae = wrapToPi(gamma - gammad);

    chidDot = vg*(Bx(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*cos(gamma)*cos(chi) ...
        + By(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*cos(gamma)*sin(chi) ...
        + Bz(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*sin(gamma));
    gammadDot = vg*(Cx(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*cos(gamma)*cos(chi) ...
        + Cy(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*cos(gamma)*sin(chi) ...
        + Cz(state(i-1, 1), state(i-1, 2), state(i-1, 3), s, vg, kappa)*sin(gamma));



    omega1 = -kChi*chie + chidDot;
    chicDot = max(min(omega1, chiDotMax), -chiDotMax);

    omega2 = -kGamma*gammae + gammadDot;
    gammacDot = max(min(omega2, gammaDotMax), -gammaDotMax);

    state(i-1, 7) = vg;
    state(i-1, 9) = chid;
    state(i-1, 11:15) = [gammad, chidDot,chicDot, gammadDot, gammacDot];

    chi = wrapToPi(chi + chicDot*dt);
    gamma = wrapToPi(gamma + gammacDot*dt);
    gamma = max(min(gamma, gammaMax), -gammaMax);

    xDot = vg*cos(gamma)*cos(chi);
    yDot = vg*cos(gamma)*sin(chi);
    zDot = vg*sin(gamma);

    state(i, 1:3) = state(i-1, 1:3) + dt*[xDot, yDot, zDot];
    state(i, 4:6) = [xDot, yDot, zDot];
    state(i, 8) = chi;
    state(i, 10) = gamma;

    if mod(i, 100) == 0 || i == (simTime/dt + 1)
        set(line, "XData", state(1:i, 1), "YData", state(1:i, 2), "ZData", state(1:i, 3));
        set(pt, "XData", state(i, 1), "YData", state(i, 2), "ZData", state(i, 3));

    
        drawnow limitrate;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end

    i = i + 1;
    textprogressbar(t*100/simTime);
end

legend(ax, "Curve to follow", "Initial Position", "Actual Trajectory", "Current Position");
hold(ax, "off");

frame = getframe(gcf);
writeVideo(writerObj, frame);
close(writerObj);

textprogressbar('done');
disp("Simulation Complete...");

figure;
subplot(3, 2, 1);
plot(0:dt:simTime, state(:, 9)*180/pi,":", lineWidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 8)*180/pi, lineWidth=1.5);

xlabel("Time [s]");
ylabel("\chi [deg]");
legend("Desired Course Angle", "Actual Course Angle");
title("Course angle variation");
hold off;

subplot(3, 2, 2);
plot(0:dt:simTime, state(:, 11)*180/pi,":", lineWidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 10)*180/pi, lineWidth=1.5);
xlabel("Time [s]");
ylabel("\gamma [deg]");
legend("Desired Climb Angle", "Actual Climb Angle");
title("Climb angle variation");
hold off;

subplot(3, 2, 3);
plot(0:dt:simTime, state(:, 12),":", linewidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 13), linewidth=1.5);
xlabel("Time [s]");
ylabel("d\chi/dt [rad/s]");
legend("desired", "commanded");
title("Course rate variation");
hold off;

subplot(3, 2, 4);
plot(0:dt:simTime, state(:, 14),":", linewidth=1.5);
hold on;
grid on; grid minor;
plot(0:dt:simTime, state(:, 15), linewidth=1.5);
xlabel("Time [s]");
ylabel("d\gamma/dt [rad/s]");
legend("desired", "commanded");
title("Climb rate variation");
hold off;

subplot(3, 2, 5);
plot(0:dt:simTime, state(:, 7), linewidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("V_g [m/s]");
title("Ground speed variation");
hold off;

subplot(3, 2, 6);
plot(0:dt:simTime, state(:, 16), lineWidth=1.5);
hold on;
grid on; grid minor;
xlabel("Time [s]");
ylabel("Cross track error, |r|");
title("Absolute cross track error variation");
hold off;
