function kineState = combinedVF2DPathGenerator(params, pathParams, state)
% %kineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
%state: [x, dx, y, dy, z, dz, phi, dphi, theta, dtheta, psi, dpsi];
    dt= params.dt; %s
    
    %UAV Specifications
    s = pathParams.s; %direction
    va = pathParams.va; %m/s
    wx = pathParams.w(1); wy = pathParams.w(2);
    
    chiDotMax = pathParams.chiDotMax; %rads/s
    
    kappa0 = pathParams.kappa0;
    kChi = pathParams.kChi;

    x = state(1); y = state(3); z = state(5);
    chi = state(11);
    
    %State Initialization
    % state(1, :) = [p0, zeros([1, 2]), 0, chi0, 0, zeros([1, 3])];
    
 
    df = calculateDerivatives2D(pathParams.fm, [x, y]');

    vg = (wx*cos(chi) + wy*sin(chi)) + sqrt(va^2 - (wx*sin(chi) - wy*cos(chi))^2);
    r = df(1);
    % state(i-1, 10) = r;
    persistent kappa;
    if(isempty(kappa))    
        kappa = kappa0;
        if r > 0.001
            kappa = kappa0/abs(r);
        end
    end

    vc = [df(2);df(3);];
    vc = vc./sqrt(vc'*vc);
    vs = [df(3);-df(2);];
    vs = vs./sqrt(vs'*vs);
    vd = -vg*tanh(kappa*r)*vc + s*vg*sech(kappa*r)*vs;

    chid = atan2(vd(2), vd(1));
    chie = wrapToPi(chi - chid);

    %coefficient functions vector of time derivative of desired course angle
    % chi_d_dot = A_x*x_dot + A_y*y_dot
    Ax = -s*kappa*df(2)*sech(kappa*df(1)) - (df(4)*df(3) - df(2)*df(6))/(df(2)^2 + df(3)^2);
    Ay = -s*kappa*df(3)*sech(kappa*df(1)) + (df(5)*df(2) - df(3)*df(6))/(df(2)^2 + df(3)^2); 

    chidDot = vg*(Ax*cos(chi) + Ay*sin(chi));

    omega1 = -kChi*chie + chidDot;
    chicDot = max(min(omega1, chiDotMax), -chiDotMax);

    chi = wrapToPi(chi + chicDot*dt);

    dx = vg*cos(chi);
    dy = vg*sin(chi);

    x = x + dx*dt;
    y = y + dy*dt;
    
    kineState = [x, dx, 0, y, dy, 0, z, 0, 0, chi, chicDot, 0]';
end



