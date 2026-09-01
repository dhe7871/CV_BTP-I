function [wsq, desEuler] = quadCopterSMC(params, state, desKineState)
    %desKineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
    xd = desKineState(1); dxd = desKineState(2); ddxd = desKineState(3);
    yd = desKineState(4); dyd = desKineState(5); ddyd = desKineState(6);
    zd = desKineState(7); dzd = desKineState(8); ddzd = desKineState(9);
    psid = wrapToPi(desKineState(10)); dpsid = desKineState(11); ddpsid = desKineState(12);

    %state: [x, dx, y, dy, z, dz, phi, dphi, theta, dtheta, psi, dpsi];
    x = state(1); dx = state(2);
    y = state(3); dy = state(4);
    z = state(5); dz = state(6);
    phi = state(7); theta = state(9); psi = state(11);
    dphi = state(8); dtheta = state(10); dpsi = state(12);

    kz = params.gains(1, 1); lambdaz = params.gains(1, 2);
    kx = params.gains(2, 1); lambdax = params.gains(2, 2);
    ky = params.gains(3, 1); lambday = params.gains(3, 2);
    kphi = params.gains(4, 1); lambdaphi = params.gains(4, 2);
    ktheta = params.gains(5, 1); lambdatheta = params.gains(5, 2);
    kpsi = params.gains(6, 1); lambdapsi = params.gains(6, 2);

    a2 = (params.I(2,2) - params.I(3, 3))/params.I(1, 1); b2 = params.I(1, 1)/params.l;
    a3 = (params.I(3, 3) - params.I(1, 1))/params.I(2, 2); b3 = params.I(2, 2)/params.l;
    a4 = (params.I(1, 1) - params.I(2, 2))/params.I(3, 3); b4 = params.I(3, 3);

    %Position Controller
    %z (Altitude)
    e = zd - z;
    de = dzd - dz;
    s = de + lambdaz*e;
    U1 = (params.m/(cos(phi)*cos(theta)))*(ddzd + params.g + lambdaz*de) + kz*sigmoid(s, params.delta);
    
    %x controller
    e = xd - x;
    de = dxd - dx;
    s = de + lambdax*e;
    Ux = (params.m/U1)*(ddxd + lambdax*de) + kx*sigmoid(s, params.delta);

    %y Controller
    e = yd - y;
    de = dyd - dy;
    s = de + lambday*e;
    Uy = (params.m/U1)*(ddyd + lambday*de) + ky*sigmoid(s, params.delta);

    %phid, thetd estimation
    phid = asin(Ux*sin(psid) - Uy*cos(psid)); phid = max(min(phid, params.phiMax), -params.phiMax);
    thetad = asin((Ux*cos(psid) + Uy*sin(psid))/cos(phid)); thetad = max(min(thetad, params.thetaMax), -params.thetaMax);

    %Attitude Controller
    %phi
    e = wrapToPi(phid - phi);
    de = -dphi;
    s = de + lambdaphi*e;
    U2 = b2*(-dphi*dpsi*a2 + lambdaphi*de) + kphi*sigmoid(s, params.delta);

    %theta
    e = wrapToPi(thetad - theta);
    de = -dtheta;
    s = de + lambdatheta*e;
    U3 = b3*(-dphi*dpsi*a3 + lambdatheta*de) + ktheta*sigmoid(s, params.delta);

    %psi
    e = wrapToPi(psid - psi);
    de = dpsid - dpsi;
    s = de + lambdapsi*e;
    U4 = b4*(ddpsid - dtheta*dphi*a4 + lambdapsi*de) + kpsi*sigmoid(s, params.delta);

    A = [
        1, 1, 1, 1;
        0, -1, 0, 1;
        -1, 0, 1, 0;
        -1, 1, -1, 1
    ];
    invA = A^(-1);
    b = [U1/params.b, U2/params.b, U3/params.b, U4/params.d]';

    wsq = invA*b;
    wsq(wsq < 0) = 0;
    wsq(wsq > params.omegaMax^2) = params.omegaMax^2;

    desEuler = [phid, thetad, psid]';
end

function x = sigmoid(s, delta)
    x = s/(abs(s) + delta);
end