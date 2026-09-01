function kineState = combinedVF3DPathGenerator(params, pathParams, kineState)
    %kineState: [x, y, z, dx, dy, dz, vg, chi, chid, gamma, gammad, dchid, dchic, dgammad, dgammac, r]';
    %state: [x, dx, y, dy, z, dz, phi, dphi, theta, dtheta, psi, dpsi];

    wx = pathParams.w(1); wy = pathParams.w(2);
    va = pathParams.va; s = pathParams.s;
    kappa0 = pathParams.kappa0;
    kChi = pathParams.kChi; kGamma = pathParams.kGamma;
    chiDotMax = pathParams.chiDotMax; gammaDotMax = pathParams.gammaDotMax;
    Bx = pathParams.B{1}; By = pathParams.B{2}; Bz = pathParams.B{3};
    Cx = pathParams.C{1}; Cy = pathParams.C{2}; Cz = pathParams.C{3};

    x = kineState(1); y = kineState(2); z = kineState(3);
    chi = kineState(8); gamma = kineState(10);
    
    % x = state(1); y = state(3); z = state(5);
    % dx = state(2); dy = state(4); dz = state(6);
    % chi = state(11); gamma = atan(dz/sqrt(dx^2 + dy^2));

    
    df = calculateDerivatives(pathParams.fm, [x, y, z]');
    dg = calculateDerivatives(pathParams.gm, [x, y, z]');

    a = wx*cos(chi) + wy*sin(chi);
    b = wx*sin(chi) - wy*cos(chi);
    vg = a*cos(gamma) ...
        + 0.5*sqrt(4*va^2 - 3*(wx^2 + wy^2) + (1 + 2*cos(2*gamma))*a^2 - b^2);
    r = sqrt(df(1)^2 + dg(1)^2);
    
    persistent kappa;
    if(isempty(kappa))
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

    chidDot = vg*(Bx(x, y, z, s, vg, kappa)*cos(gamma)*cos(chi) ...
        + By(x, y, z, s, vg, kappa)*cos(gamma)*sin(chi) ...
        + Bz(x, y, z, s, vg, kappa)*sin(gamma));
    gammadDot = vg*(Cx(x, y, z, s, vg, kappa)*cos(gamma)*cos(chi) ...
        + Cy(x, y, z, s, vg, kappa)*cos(gamma)*sin(chi) ...
        + Cz(x, y, z, s, vg, kappa)*sin(gamma));

    omega1 = -kChi*chie + chidDot;
    chicDot = max(min(omega1, chiDotMax), -chiDotMax);

    omega2 = -kGamma*gammae + gammadDot;
    gammacDot = max(min(omega2, gammaDotMax), -gammaDotMax);

    chi = wrapToPi(chi + chicDot*params.dt);
    gamma = wrapToPi(gamma + gammacDot*params.dt);
    gamma = max(min(gamma, pathParams.gammaMax), -pathParams.gammaMax);

    dx = vg*cos(gamma)*cos(chi);
    dy = vg*cos(gamma)*sin(chi);
    dz = vg*sin(gamma);

    kineState(1:3) = [x, y, z]' + params.dt*[dx, dy, dz]';
    kineState(4:6) = [dx, dy, dz]'; kineState(7) = vg;
    kineState([8, 9]) = [chi, chid]';
    kineState([10, 11]) = [gamma, gammad]';
    kineState(12:15) = [chidDot, chicDot, gammadDot, gammacDot]';
    kineState(16) = r;
end