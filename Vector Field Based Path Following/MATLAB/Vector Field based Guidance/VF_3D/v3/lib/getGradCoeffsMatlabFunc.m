function [B, C] = getGradCoeffsMatlabFunc(f, g)
    %%coefficient functions vector of time derivative of desired course and climb angle
    %%chi_d_dot = B_x*x_dot + B_y*y_dot + B_Z*z_dot
    %%gamma_d_dot = C_x*x_dot + C_y*y_dot + C_z*z_dot
    syms x y z s vg kappa;

    fx = diff(f, x);
    fy = diff(f, y);
    fz = diff(f, z);

    gx = diff(g, x);
    gy = diff(g, y);
    gz = diff(g, z);

    vc = [f*fx + g*gx, f*fy + g*gy, f*fz + g*gz]';
    vc = vc./sqrt(vc'*vc);

    vs = [fy*gz - fz*gy; fz*gx - fx*gz; fx*gy - fy*gx;];
    vs = vs./sqrt(vs'*vs);

    vd = -vg*tanh(kappa*sqrt(f^2 + g^2))*vc + s*vg*sech(kappa*sqrt(f^2 + g^2))*vs;
    chid = atan(vd(2)/vd(1));
    gammad = atan(vd(3)/sqrt(vd(1)^2 + vd(2)^2));

    B = {
        matlabFunction(diff(chid, x), "Vars", [x, y, z, s, vg, kappa]),...
        matlabFunction(diff(chid, y), "Vars", [x, y, z, s, vg, kappa]),...
        matlabFunction(diff(chid, z), "Vars", [x, y, z, s, vg, kappa])
    };
    C = {
        matlabFunction(diff(gammad, x), "Vars", [x, y, z, s, vg, kappa]),...
        matlabFunction(diff(gammad, y), "Vars", [x, y, z, s, vg, kappa]),...
        matlabFunction(diff(gammad, z), "Vars", [x, y, z, s, vg, kappa])
    };
end
