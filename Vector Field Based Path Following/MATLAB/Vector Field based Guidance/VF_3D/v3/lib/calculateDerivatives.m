function df = calculateDerivatives(f, p)
    x = p(1); y = p(2); z = p(3);

    fm = f{1};
    fxm = f{2}; fym = f{3}; fzm = f{4};
    fxxm = f{5}; fyym = f{6}; fzzm = f{7};
    fxym = f{8}; fyzm = f{9}; fzxm = f{10};
    
    df = [
        fm(x, y, z),...
        fxm(x, y, z), fym(x, y, z), fzm(x, y, z),...
        fxxm(x, y, z), fyym(x, y, z), fzzm(x, y, z),...
        fxym(x, y, z), fyzm(x, y, z), fzxm(x, y, z)
    ]';
end