function y = direct(b, a, x)
    y(1) = b(1)*x(1);
    y(2) = b(1)*x(2) + b(2)*x(1) - a(2)*y(1);
    
    for i = 3:length(x)
        y(i) = b(1)*x(i) + b(2)*x(i-1) + b(3)*x(i-2)-a(2)*y(i-1)-a(3)*y(i-2);
    end
end
