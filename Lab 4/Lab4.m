%{
%% Problem 1.2

x = [1, 1, 2, 1, 2, 2, 1, 1]';
b = [4, 2.4, -1.6];
a = [1, -0.5, 0.6];
vin = [0;0];

ydirect = direct(b, a, x);
ytran = tran(b, a, x);

[yfilter, voutfilter] = filter(b, a, x);
yfilter = yfilter';

display(ydirect);
display(ytran);
display(yfilter);

[yfilter, voutfilter] = filter(b, a, x);
[ytran, vouttran] = tran(b, a, x);
yfilter = yfilter';

display(vouttran);
display(voutfilter);

%% Problem 1.3

fprintf(' n    x       y       v1        v2\n') 
fprintf('--------------------------------------\n')

for i = 1:length(x) 
    [ytran, vouttran] = tran(b, a, x(i), vin); 
    fprintf('%2i %4i %9.4f %9.4f %9.4f \n', i-1, x(i), ytran(1), vin(1), vin(2)); 
    vin = vouttran;
end

fprintf('%2i %4s %7s   %9.4f %9.4f \n', 8, '-', '-', vin(1), vin(2));
%}
%% Problem 2.1 && 2.2
f = linspace(0,10,1001);
f0 = 4;
fs = 200;
w0 = (2*pi*f0)/fs;
w = (2*pi*f)/fs;
R1 = 0.980;
R2 = 0.995;

G1 = (1 - 2*R1*cos(w0) + R1^2)/(1 - 2*cos(w0) + 1); 
G2 = (1 - 2*R2*cos(w0) + R2^2)/(1 - 2*cos(w0) + 1);

b1 = G1*[1, -2*cos(w0), 1];
a1 = [1, -2*R1*cos(w0), R1^2];
b2 = G2*[1, -2*cos(w0), 1];
a2 = [1, -2*R2*cos(w0), R2^2];

fresp = @(b,a,w) polyval(flip(b),exp(-j*w))./polyval(flip(a),exp(-j*w));
mag = @(b,a,w) abs(fresp(b,a,w));

fprintf('filter 1\n');
fprintf('-----------------------------------\n');
fprintf('b = [ %1.6f %4.6f %4.6f] \n', b1(1), b1(2), b1(3));
fprintf('a = [ %1.6f %4.6f %4.6f] \n\n', a1(1), a1(2), a1(3)); 

fprintf('filter 2\n'); 
fprintf('-----------------------------------\n'); 
fprintf('b = [ %1.6f %4.6f %4.6f] \n', b2(1), b2(2), b2(3)); 
fprintf('a = [ %1.6f %4.6f %4.6f] \n', a2(1), a2(2), a2(3));

%plot of filter 1
f1 = 2;
f2 = 6;
fa = fs/pi*(1-R1);
fLa = f0 - 0.5*fa; 
fRa = f0 + 0.5*fa; 
fL =[fLa+1 fLa-1]; 
fR =[fRa+1 fRa-1];
F = @(f) mag(b1,a1,2*pi*f/fs) - 1/sqrt(2);
fLexact = fzero(F, fL);
fRexact = fzero(F, fR); 
wleft = fLexact*2*pi/fs; 
wright = fRexact*2*pi/fs; 

figure; 
plot(f, G1*mag(b1, a1, w));
hold on;
plot(fLexact, mag(b1, a1, wleft), 'g.', fRexact, mag(b1, a1, wright), 'g.', f1,mag(b1, a1, f1*2*pi/fs), 'r.', f2, mag(b1, a1, f2*2*pi/fs),'r.', f0, mag(b1, a1, f0*2*pi/fs), 'r.'); 
plot([fLexact, fRexact], [mag(b1, a1, wleft), mag(b1, a1, wright)]);
title('magnitude response R1');
xlim([0 10]);
ylim([0 1.1]);
hold off;

error = 100*(abs(fLexact - fLa) + abs(fRexact - fRa))/(fLexact + fRexact);
fprintf('\n\nF1: exact approx\n');
fprintf('------------------\n');
fprintf('fL = %.4f %.4f \n', fLexact, fLa);
fprintf('fR = %.4f %.4f \n', fRexact, fRa);
fprintf('------------------\n')
fprintf('percent error = %.4f%%\n', error);

fa = (fs/pi)*(1-R2);
fLa = f0 - 0.5*fa;
fRa = f0 + 0.5*fa;
fL = [fLa-0.1 fLa+0.1];
fR = [fRa-0.1 fRa+0.1];
F = @(f) mag(b2,a2,2*pi*f/fs) - 1/sqrt(2);
fLexact = fzero(F,fL);
fRexact = fzero(F,fR);
wleft = fLexact*2*pi/fs;
wright = fRexact*2*pi/fs;

figure;
plot(f, mag(b2,a2,w));
hold on;
plot(fLexact, G2*mag(b2, a2, wleft), 'g.', fRexact, G2*mag(b2, a2, wright), 'g.', f1, mag(b2, a2, f1*2*pi/fs), 'r.', f2, mag(b2, a2, f2*2*pi/fs), 'r.', f0, mag(b2, a2, f0*2*pi/fs), 'r.');
plot([fLexact, fRexact], [mag(b2,a2,wleft),mag(b2,a2,wright)], 'g-')
title('magnitude response R2'); 
xlim([0 10]);
ylim([0 1.1]);
hold off;

error = 100*(abs(fLexact - fLa) + abs(fRexact - fRa))/(fLexact + fRexact);
fprintf('\n\nF2: exact approx\n');
fprintf('------------------\n');
fprintf('fL = %.4f %.4f \n', fLexact, fLa);
fprintf('fR = %.4f %.4f \n', fRexact, fRa);
fprintf('------------------\n');
fprintf('percent error = %.4f%%\n', error);

%% Problem 2.3

f0 = 4; 
f1 = 2; 
f2 = 6;
fs = 200;
Ts = 1/fs;

x=@(t) cos(2*pi*f1.*t).*(t >= 0 & t < 4) + cos(2*pi*f0.*t).*(t >= 4 & t < 8) + cos(2*pi*f2.*t).*(t >= 8 & t < 12);

t = linspace(0,12,1001);

figure;
plot(t, x(t));
hold on;
title('input signal');
xlim([0 12]);
ylim([-2 2]);
hold off;

G1 = (1 - 2*R1*cos(w0) + R1^2)/(1 - 2*cos(w0) + 1); 
G2 = (1 - 2*R2*cos(w0) + R2^2)/(1 - 2*cos(w0) + 1);

b1 = G1*[1, -2*cos(w0), 1];
a1 = [1, -2*R1*cos(w0), R1^2];
b2 = G2*[1, -2*cos(w0), 1];
a2 = [1, -2*R2*cos(w0), R2^2];

t = 0:Ts:12;

figure;
y = tran(G1*b1, a1, x(t));
plot(t, y);
hold on;
title('Notch Filter Output, R = 0.980');
xlim([0 12]);
ylim([-2 2]);
hold off;

figure;
y = tran(G2*b2, a2, x(t));
plot(t, y);
hold on;
title('Notch Filter Output, R = 0.995');
xlim([0 12]);
ylim([-2 2]);
hold off;

t_eff1 = Ts * log(0.01)/log(R1);
t_eff2 = Ts * log(0.01)/log(R2);

fprintf(' R       t_eff\n');
fprintf('--------------------\n');
fprintf('%.4f  %2.4f (sec)\n', R1, t_eff1);
fprintf('%.4f  %2.4f (sec)\n\n', R1, t_eff2);