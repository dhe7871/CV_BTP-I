function df = calculateDerivatives2D(fm, p)
    x = p(1);
    y = p(2);

    f = fm{1};
    fx = fm{2};
    fy = fm{3};
    fxx = fm{4};
    fyy = fm{5};
    fxy = fm{6};

    df = [f(x, y), fx(x, y), fy(x, y), fxx(x, y), fyy(x, y), fxy(x, y)];
end