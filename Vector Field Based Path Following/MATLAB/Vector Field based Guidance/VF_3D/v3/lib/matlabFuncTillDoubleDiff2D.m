function fnm = matlabFuncTillDoubleDiff2D(fn)
    syms x y;
    fx = diff(fn, x);
    fy = diff(fn , y);
    fxx = diff(fx, x);
    fyy = diff(fy, y);
    fxy = diff(fx, y);

    fm = matlabFunction(fn, "Vars",[x, y]);
    fxm = matlabFunction(fx, "Vars",[x, y]);
    fym = matlabFunction(fy, "Vars",[x, y]);
    fxxm = matlabFunction(fxx, "Vars",[x, y]);
    fyym = matlabFunction(fyy, "Vars",[x, y]);
    fxym = matlabFunction(fxy, "Vars",[x, y]);

    fnm = {fm, fxm, fym, fxxm, fyym, fxym};
end
