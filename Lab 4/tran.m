function [y, vout] = tran(b, a, x, vin)
    if (nargin < 4)
        vin = [0;0];  
    end

    v1(1) = vin(1);
    v2(1) = vin(2);

    for i = 1:length(x)
        y(i) = b(1)*x(i) + v1(i);
        v1(i+1) = b(2)*x(i) - a(2)*y(i) + v2(i);
        v2(i+1) = b(3)*x(i) - a(3)*y(i);
    end

    vout = [v1(length(x)+1); v2(length(x)+1)];
end