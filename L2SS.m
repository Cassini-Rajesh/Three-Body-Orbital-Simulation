function dxdt = L2SS(t, y, MU1)
    x1 = y(1);
    x2 = y(2);
    x3 = y(3);
    x4 = y(4);

    a = (x1+MU1)*(1-MU1);
    b = (x1+MU1-1)*(MU1);
    c = x2*(1-MU1);
    d = x2*MU1;

    p1 = sqrt( (x1 + MU1)^2 + x2^2 );
    p2 = sqrt( (x1 + MU1 - 1)^2 + x2^2 );

    Ux = x1 - ( (a) / (p1^3) ) - ( (b) / (p2^3) );
    Uy = x2 - ( (c) / (p1^3) ) - ( (d) / (p2^3) );

    dxdt = zeros(4,1);
    dxdt(1) = x3;
    dxdt(2) = x4;
    dxdt(3) = Ux + 2*x4;
    dxdt(4) = Uy - 2*x3;
end