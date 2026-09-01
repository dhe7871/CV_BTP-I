clc; clear; close all;

va = 10;
w = [3, 4, 5];
chi = 0;
gamma = -pi/3;
unitVg = [cos(gamma)*cos(chi), cos(gamma)*sin(chi), sin(gamma)];

figure;
hold on; grid on; grid minor; axis equal;
%axes
plot3([0, 2*va], [0, 0], [0, 0], "LineStyle", "--", "LineWidth", 1.5);
plot3([0, 0], [0, 2*va], [0, 0], "LineStyle", "--", "LineWidth", 1.5);
plot3([0, 0], [0, 0], [0, 2*va], "LineStyle", "--", "LineWidth", 1.5);


% wind vector
plot3([0, w(1)], [0, w(2)], [0, w(3)], "LineWidth", 1.5);
scatter3(w(1), w(2), w(3), "Marker", "o", "LineWidth", 1.5, "MarkerFaceColor","auto");


[X,Y,Z] = sphere(50);   % higher resolution: (n+1)-by-(n+1)
X = X * va + w(1);
Y = Y * va + w(2);
Z = Z * va + w(3);
surf(X,Y,Z, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
camlight; lighting phong


% ground velocity line
plot3([0, unitVg(1)*2*va], [0, unitVg(2)*2*va], [0, unitVg(3)*2*va], "LineWidth", 1.5, "LineStyle", ":");


A = unitVg(1)^2 + unitVg(2)^2 + unitVg(3)^2;
B = -2*(unitVg(1)*w(1) + unitVg(2)*w(2) + unitVg(3)*w(3));
C = w(1)^2 + w(2)^2 + w(3)^2 - va^2;
root = roots([A, B, C]);
magVg = root(root >= 0);
magVg = magVg(1);

Vg = unitVg * magVg;

scatter3(Vg(1), Vg(2), Vg(3), "Marker", "o", "LineWidth", 1.5, "MarkerFaceColor","auto");

plot3([w(1), Vg(1)], [w(2), Vg(2)], [w(3), Vg(3)], "LineWidth", 1.5);
xlabel("x");
ylabel("y");
zlabel("z");
hold off;
