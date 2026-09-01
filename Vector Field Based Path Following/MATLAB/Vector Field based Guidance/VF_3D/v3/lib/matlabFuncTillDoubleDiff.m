function fnm =  matlabFuncTillDoubleDiff(fn)
    syms x y z;

    fx = diff(fn, x);
    fy = diff(fn, y);
    fz = diff(fn, z);
    
    fxx = diff(fx, x);
    fyy = diff(fy, y);
    fzz = diff(fz, z);
    fxy = diff(fx, y);
    fyz = diff(fy, z);
    fzx = diff(fz, x);

    fm = matlabFunction(fn, "Vars",[x,y,z]);

    fxm = matlabFunction(fx, "Vars",[x,y,z]);
    fym = matlabFunction(fy, "Vars",[x,y,z]);
    fzm = matlabFunction(fz, "Vars",[x,y,z]);
    fxxm = matlabFunction(fxx, "Vars",[x, y, z]);
    fyym = matlabFunction(fyy, "Vars",[x, y, z]);
    fzzm = matlabFunction(fzz, "Vars",[x, y, z]);
    fxym = matlabFunction(fxy, "Vars",[x, y, z]);
    fyzm = matlabFunction(fyz, "Vars",[x, y, z]);
    fzxm = matlabFunction(fzx, "Vars",[x, y, z]);

    fnm = {fm, fxm , fym, fzm, fxxm, fyym, fzzm, fxym, fyzm, fzxm};
end

