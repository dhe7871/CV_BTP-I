%F & G are MATLAB functions, x0: approximate starting point, steps: No. of
%steps to generate the length of the curve
function curve = generate3DCurvePoints(F, G, x0, steps)
    %wrapping into vector form for fsolve
    fun = @(X) [F(X(1), X(2), X(3)); G(X(1), X(2), X(3))];
    
    opts = optimoptions("fsolve", "Display", "none","Algorithm", "levenberg-marquardt");
    P0 = fsolve(fun, x0, opts);

    step = 0.05;
    Nsteps = steps;

    curve = zeros(3, Nsteps);
    curve(:, 1) = P0;
    
    tolClose = 1e-2;
    for k = 2:Nsteps
        gradF = grad(F, curve(:, k-1)); %normal to F
        gradG = grad(G, curve(:, k-1)); %normal to G

        %tangent direction
        tdir = cross(gradF, gradG);
        tdir = tdir/norm(tdir);

        Ppred = curve(:, k-1) + step*tdir;

        Pnew = fsolve(fun, Ppred, opts);

        curve(:, k) = Pnew;
        if(norm(Pnew - P0) < tolClose && k > Nsteps/10)
            curve = curve(:, 1:k);
            break;
        end
    end
end

%gradients via finite difference
function g = grad(fun, P)
    h = 1e-6;
    x = P(1); y = P(2); z = P(3);
    g = [(fun(x+h, y, z) - fun(x-h, y, z))/(2*h);
         (fun(x, y+h, z) - fun(x, y-h, z))/(2*h);
         (fun(x, y, z+h) - fun(x, y, z-h))/(2*h)];
end