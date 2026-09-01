function desKineState = generateDesKineStateVF3D(pathParams, kineState)
    %kineState: [x, y, z, dx, dy, dz, vg, chi, chid, gamma, gammad, dchid, dchic, dgammad, dgammac, r];
    x = kineState(1); y = kineState(2); z = kineState(3);
    dx = kineState(4); dy = kineState(5); dz = kineState(6);
    vg = kineState(7); chi = kineState(8); gamma = kineState(10);
    dchi = kineState(13); dgamma = kineState(15);

    wx = pathParams.w(1); wy = pathParams.w(2); %wind velocity

    %%desKineState: [xd, dxd, ddxd, yd, dyd, ddyd, zd, dzd, ddzd, psid, dpsid, ddpsid];
    %rate of change of ground velocity
    a = wx*cos(chi) + wy*sin(chi);
    b = wx*sin(chi) - wy*cos(chi);
    dvg = -b*cos(gamma)*dchi - a*sin(gamma)*dgamma...
        - (dgamma*sin(2*gamma)*a^2 + dchi*(1 + cos(2*gamma))*a*b)/(2*(vg - a*cos(gamma)));
    ddx = dvg*dx/vg - dz*cos(chi)*dgamma - dy*dchi;
    ddy = dvg*dy/vg - dz*sin(chi)*dgamma + dx*dchi;
    ddz = dvg*dz/vg + vg*cos(gamma)*dgamma;

    desKineState = [
        x, dx, 0,...
        y, dy, 0,...
        z, dz, 0,...
        chi, dchi, 0
    ]';
end