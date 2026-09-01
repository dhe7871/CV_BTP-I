function kineState = NLGL2DPathGenerator(params, pathParams, state)
% %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
%state: [x, dx, y, dy, z, dz, phi, dphi, theta, dtheta, psi, dpsi];

    %NLGL 2D line
    w1 = pathParams.w1;
    w2 = pathParams.w2;
    
    lineSlope = (w2(2) - w1(2))/(w2(1) - w1(1));
    lineIntercept = w1(2) - lineSlope * w1(1);
    
    x = state(1); y = state(3);
    z = state(5); dz = state(6);
    psi = state(11);
    
    L = 100;            %finding radius L (can change)
    va = 8;            %UAV velocity (assuming constant throughout) (can change)
    back = false;        %Tracing the line backward and forward (can change)
    
    phiMax = params.phiMax; %radians
    psiDotMax = params.g*tan(phiMax)/va; %radians/sec
    
    % lineParams = [m, c]
    lineParams = [lineSlope, lineIntercept];
    
    % circleParams = [x0, y0, r]
    circleParams = [x, y, L];
    [ptsExists, VTP] = getIntersectionPts(circleParams, lineParams, back);
    if ptsExists
        psid = atan2(VTP(2) - y, VTP(1) - x);
    else
        psid = atan(-1/lineSlope);
        if  y - lineSlope*x - lineIntercept < 0
            if psid < 0
                psid = pi + psid;
            end
        else
            if psid >= 0
                psid = -pi + psid;
            end
        end
    end
    
    de = wrapToPi(psid - psi);
    dpsi = max(min(2*va*sin(de)/L, psiDotMax), -psiDotMax); %Non-Linear guidance law: psiDot = 2*va*sin(psid - psi)/L
    psi = wrapToPi(psi + dpsi*params.dt);  %psi

    dx = va*cos(psi);
    dy = va*sin(psi);
    ddx = -va*sin(psi)*dpsi;
    ddy = va*cos(psi)*dpsi;

    x = x + dx*params.dt;
    y = y + dy*params.dt;


    kineState = [x, dx, ddx, y, dy, ddy, z, dz, 0, psi, dpsi, 0]';
end

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