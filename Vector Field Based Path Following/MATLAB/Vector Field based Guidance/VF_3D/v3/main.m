clc; clear; close all; clear kappa;
addpath(fullfile(pwd, "lib"));

simTime = 300;
dt = 0.001;

%%%%%Quadcopter Parameters [Q250]%%%%%
params.dt = dt;
params.simTime = simTime;
params.g = 9.81;    %[m/s^2]
params.m = 0.444;   %[Kg]
params.l = 0.177;   %[m]
params.I = diag([0.001758, 0.001758, 0.003493]);
params.invI = diag(1./[0.001758, 0.001758, 0.003493]);
% params.Jr = 3e-6;
params.b = 9.92e-7;
params.d = 1.49e-8;
params.omegaMax = 18510*2*pi/60;
params.phiMax = 60*pi/180;
params.thetaMax = 60*pi/180;
    
%   k, lambda    
params.gains = [
    % 5, 0.5;             %[z]  {for ddx, ddy, ddz ~= 0}
    % 0.04, 0.35;         %[x]
    % 0.04, 0.35;         %[y]
    % 1.5, 2.5;           %[phi]
    % 1.5, 2.5;           %[theta]
    % 1.2, 0.7;           %[psi]
     5, 0.5;             %[z]
    0.35, 1.5;         %[x]
    0.35, 1.5;         %[y]
    2.5, 5;           %[phi]
    2.5, 5;           %[theta]
    1.2, 0.7;           %[psi]
];

params.delta = 0.3;

pathChoice = input("Path Choice: ");
%%%%%Path parameters [pathParams]%%%%%
syms x y z;

%%3D Path Definition
pathParams.pathChoice = pathChoice;
pathParams.prevCurve = false;

% PathChoice = 3 => Combined VF2D
% pathParams.fn = ((x - 2)^2 + (y + 10)^2 - 9775)*10; %kappa = 10
% pathParams.fn = 500*(y - 3*x + 10);     %kappa0 = 8
% pathParams.fn = (y - 200*sin(x/50))*150; %kappa0 = 5
% pathParams.fn = (1/((5*x^2 + x*y + 5*y^2)^4.5 + 1))*(-(5*x^2 + x*y + 5*y^2)^4.5 ...
%     + 300*sqrt(2)*(5*x^2 + x*y + 5*y^2)^4 ...
%     + 45*sqrt(22)*(x^2 - y^2)*(x^2 + 20*x*y + y^2)*(49*x^4 - 20*x^3*y - 300*x^2*y^2 - 20*x*y^3 + 49*y^4));



%PathChoice = 2 => Combined VF3D
% pathParams.fn = ((x - 2)^2 + (y + 10)^2 + (z - 15)^2 - 10000)/100; %kappa0 = 10
% pathParams.gn = 3*x + y + 5*z - 350;
% pathParams.fn = y + z - 100;
% pathParams.gn = 200*sin((x-y)/100) - z;
% pathParams.fn = x - 200*cos(z/20);                              %kappa0 = 8
% pathParams.gn  = y - 150*sin(z/20);
% pathParams.fn = 100*(x + y + z)/(sqrt(5*x^2 + x*y + 5*y^2) + 1); %kappa0 = 5
% pathParams.gn = 100*(1/((5*x^2 + x*y + 5*y^2)^4.5 + 1))*(-(5*x^2 + x*y + 5*y^2)^4.5 ...
%     + 300*sqrt(2)*(5*x^2 + x*y + 5*y^2)^4 ...
%     + 45*sqrt(22)*(x^2 - y^2)*(x^2 + 20*x*y + y^2)*(49*x^4 - 20*x^3*y - ...
%     300*x^2*y^2 - 20*x*y^3 + 49*y^4));

pathParams.fn = (x/10)^2 + (y/10)^2 + (z/5)^2 - 1;
pathParams.gn = x + y + z - 1;


%%UAV condition and Other settings
pathParams.s = 1; %direction
pathParams.va = 10; %m/s
pathParams.w = [3, 4, 0]'; %m/s (can only deal with wind in x and y direction => wz = 0)
pathParams.gammaMax = pi/3;
pathParams.chiDotMax = 0.5; %rads/s
pathParams.gammaDotMax = 0.3; %rad/s

%%Path convergence Parameters for Path generator
pathParams.kappa0 = 1000;
pathParams.kChi = 1;
pathParams.kGamma = 1;

%%Initial State of the UAV
pathParams.gamma0 = -pi/4; %rad
pathParams.chi0 = -3*pi/4; %rad
% pathParams.p0 = [400, -200, 100]';
pathParams.p0 = [200, 200, 150]';

%PathChoice = 1 => NLGL Line Path Parameters (Uncomment when using NLGL 2D line path
%generator)
pathParams.w1 = [0, 0]';
pathParams.w2 = [100, 1000]';


%state: [x, dx, y, dy, z, dz, phi, dphi, theta, dtheta, psi, dpsi];
state = zeros(12, simTime/dt + 1);
state([1, 3, 5], 1) = pathParams.p0;
state(11, 1) = pathParams.chi0;

miscState = zeros(7, simTime/dt + 1);

if(pathChoice == 1) %NLGL 2D line
    %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
    kineState = zeros(12, simTime/dt + 1);
    kineState([1, 4, 7], 1) = pathParams.p0;
    kineState(11, 1) = pathParams.chi0;
elseif(pathChoice == 3) %Combined VF 2D
    %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
    kineState = zeros(12, simTime/dt + 1);
    kineState([1, 4, 7], 1) = pathParams.p0;
    kineState(11, 1) = pathParams.chi0;

    %CombinedVF2DPath preprocessing

    %differentiated matlab functions
    fm = matlabFuncTillDoubleDiff2D(pathParams.fn);
    pathParams.fm = fm;

elseif(pathChoice == 2) %Combined VF 3D
    %kineState: [x, y, z, dx, dy, dz, vg, chi, chid, gamma, gammad, dchid, dchic, dgammad, dgammac, r];
    kineState = zeros(16, simTime/dt + 1);
    kineState(1:3, 1) = pathParams.p0;
    kineState(8, 1) = pathParams.chi0;
    kineState(10, 1) = pathParams.gamma0;

    %CombinedVF3DPath pre processing

    %%coeffiecient functions vector of time derivative of desired course and climb angle
    %%chi_d_dot = B_x*x_dot + B_y*y_dot + B_Z*z_dot
    %%gamma_d_dot = C_x*x_dot + C_y*y_dot + C_z*z_dot
    disp("Calculating Gradient Coeffients for the Path Generator...");
    [B, C] = getGradCoeffsMatlabFunc(pathParams.fn, pathParams.gn);
    disp("Gradient Coefficient calculation complete...");

    %%Differentiated Matlab functions
    fm = matlabFuncTillDoubleDiff(pathParams.fn);
    gm = matlabFuncTillDoubleDiff(pathParams.gn);

    if(pathParams.prevCurve && exist(fullfile(pwd, "temp/curve.mat"), "file"))
        load(fullfile(pwd, "temp/curve.mat"), "curve");
    else
        disp("Generating 3D Curve points, might take some time...");
        curve = generate3DCurvePoints(fm{1}, gm{1}, pathParams.p0, 50000);
        save(fullfile(pwd, "temp/curve.mat"), "curve");
        disp("3D Curve point generation complete...");
    end

    pathParams.B = B; pathParams.C = C;
    pathParams.fm = fm; pathParams.gm = gm;
    pathParams.curve = curve;
end

disp("Quadcopter Simulation Started...");
clear textprogressbar;
textprogressbar('Progress: ');
i = 2;
for t=dt:dt:simTime
    if(pathChoice == 1)
        %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid]';
        kineState(:, i) = NLGL2DPathGenerator(params, pathParams, state(: , i-1));
        %%desKineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid]';
        desKineState = kineState(:, i);
    elseif(pathChoice == 3)
        %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid]';
        kineState(:, i) = combinedVF2DPathGenerator(params, pathParams, state(:, i-1));
        %%desKineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid]';
        desKineState = kineState(:, i);
    elseif(pathChoice == 2)
        %kineState: [x, y, z, dx, dy, dz, vg, chi, chid, gamma, gammad, dchid, dchic, dgammad, dgammac, r]';
        kineState(:, i) = combinedVF3DPathGenerator(params, pathParams, kineState(: , i-1));
        %%desKineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid]';
        desKineState = generateDesKineStateVF3D(pathParams, kineState(:, i));
    end

    [wsq, desEuler] = quadCopterSMC(params, state(:, i-1), desKineState);

    U1 = params.b*(wsq(1) + wsq(2) + wsq(3) + wsq(4));    %temp
    U2 = params.b*(-wsq(2) + wsq(4));
    U3 = params.b*(-wsq(1) + wsq(3));
    U4 = params.d*(-wsq(1) + wsq(2) - wsq(3) + wsq(4));
    
    miscState(1:7, i) = [sqrt(wsq')*60/(2*pi), desEuler']';

    state(:, i) = quadCopterDynamics(params, state(:, i-1), [U1, U2, U3, U4]);
    
    if any(~isfinite(state(:,i))) || any(~isfinite(kineState(:,i))) || any(~isfinite(wsq))
        error('NaN/Inf at step %d (t=%.3f s)', i, t);
    end

    textprogressbar(t*100/params.simTime);
    i = i + 1;
end
textprogressbar('done');
disp("Quadcopter Simulation Complete...");

disp("Plotting figures...");
plotFigures(params, pathParams, kineState, state, miscState);
disp("Figure plotting complete...");