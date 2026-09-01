function [Rbi, Rib, C, invC] = rotationMatrices(phi, theta, psi)
    R3 = [
        cos(psi), -sin(psi), 0;
        sin(psi), cos(psi), 0;
        0, 0, 1;
    ];
    R2 = [
        cos(theta), 0, sin(theta);
        0, 1, 0;
        -sin(theta), 0, cos(theta);
    ];
    R1 = [
        1, 0, 0;
        0, cos(phi), -sin(phi);
        0, sin(phi), cos(phi);
    ];

    Rbi = R3 * R2 * R1; %transformation matrix to transform body to inertial frame
    Rib = Rbi';

    C = [
        1, 0, -sin(theta);
        0, cos(phi), sin(phi)*cos(theta);
        0, -sin(phi), cos(phi)*cos(theta);
    ];                 %transformation matrix for euler rates to body omega

    invC = [
        1, sin(phi)*tan(theta), cos(phi)*tan(theta);
        0, cos(phi), -sin(phi);
        0, sin(phi)*sec(theta), cos(phi)*sec(theta);
    ];                  %transformation matrix for body omega to euler rates
end