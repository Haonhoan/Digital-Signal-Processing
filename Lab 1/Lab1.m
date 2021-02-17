%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1

%Problem 1.1
%x conditions
L = 200;
K = 50;

%h impulse response
syms t
h = (1/15).* ones(1,15);

%x impulse response
n = 0:L-1;
x = double(rem(n,K) < K/2);

%convolution
y = conv(h, x);

%plot of graph
figure
plot(n,y(1:200))
hold on;
plot(n, x, '--');
title("Problem 1.1");
xlim([0,200]);
ylim([-0.5,1.5]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1.2
%x conditions
L = 200;
K = 50;

%h impulse response
t = 0:1:14;
h = (0.25*(0.75.^t));

%x impulse response
n = 0:L-1;
x = double(rem(n,K) < K/2);

%convolution
y = conv(h, x);

%plot of graph
figure
plot(n,y(1:200))
hold on;
plot(n, x, '--');
title("Problem 1.2");
xlim([0,200]);
ylim([-0.5,1.5]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1.3
%h impulse response
t = 0:1:24;
h = (0.95.^t);

%x impulse response
n = 0:1:120;
d = @(n) double(n==0);
x = d(n) + 2*d(n-40) + 2*d(n-70) + d(n-80);

%convolution
y = conv(h, x);

%plot of graph
figure
plot(n,y(1:121));
hold on;
plot(n, x, '--');
title("Problem 1.3");
xlim([0,120]);
ylim([-0.1,3]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 2

%Problem 2.1
n = 0;
d = @(n) double(n==0);
x = d(n);

y = syst(x, 179008726);

disp(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 2.2
n = 0:1:2;
d = @(n) double(n==0);
x = d(n-2);

y = syst(x, 179008726);

disp(y);

x = 3*d(n) + 2*d(n-2);

y = syst(x, 179008726);

disp(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 2.3

id = 179008726;
n = 0:1:2;
h = [1 -1 2 3 4];
d = @(n) double(n==0);
x = d(n);

x = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

y = syst(x,id);
disp(y);
y = conv(x,h);
disp(y);

x = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1];

y = syst(x,id);
disp(y);
y = conv(x,h);
disp(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3

%Problem 3.1
p = @(n) double(n>=0 & n<=L-1);
L = 40;
n = -5:1:45;

%plot of graph
figure;
stem(n, p(n));
hold on;
title("Problem 3.1 - Part 1");
xlim([-5, 45]);
ylim([0, 1.5]);
hold off;

x = linspace(-pi,pi,1001);
P = @(x) L*exp(-1i.*x*(L-1)/2).*(sinc(x.*L./(2*pi))./sinc(x./(2*pi)));
w0 = 0;
f = abs(P(x)/P(w0));
freq = freqz(p(n),1,x);

%plot of graph
figure;
plot(x/pi,f);
hold on;
title("Problem 3.1 - Part 2");
xlim([-1,1]);
ylim([0,1]);
hold off;

error1 = norm(freq);
error2 = norm(P(x));
error = error1 - error2;
disp(error1);
disp(error2);
disp(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3.2

L = 40;
w0 = 0.2*pi;
n = -5:45;
x = linspace(-pi,pi,1001);

s = @(n) sin(w0.*n).*double(n >= 0 & n <= L - 1);

figure;
stem(n,s(n));
hold on;
xlim([-5,45]);
ylim([-1.5,1.5]);
title("Problem 3.2 - Part 1");
hold off;

P = @(x) L*exp(-1i.*x*(L-1)/2).*(sinc(x.*L./(2*pi))./sinc(x./(2*pi)));
p = @(n) double(n >= 0 & n <= L - 1);

S = @(x) (1/(2*1i))*(P(x - w0) - P(x + w0));
f = abs(S(x)/S(w0));
freq = freqz(s(n),1,x);

figure;
plot(x/pi,f);
hold on;
xlim([-1,1]);
ylim([0,1]);
title("Problem 3.2 - Part 2");
hold off;

error1 = norm(freq);
error2 = norm(S(x));
error = error1 - error2;
disp(error1);
disp(error2);
disp(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3.3

w1 = .2*pi;
w2 = .4*pi;
n = -5:1:45;

p = @(n) double(n>=0 & n<=L-1);
s = @(n) (sin(w1*n) + .8.*(sin(w2*n))).*p(n);

figure;
stem(n, s(n));
hold on;
xlim([-5,45]);
ylim([-1.5,1.5]);
title('Problem 3.3 - Part 1');
hold off;

L = 40;
x = linspace(-pi,pi,1001);


P = @(x) L*exp( (-1i*(L-1)/2 ).*(x)).*sinc(x.*L./(2*pi))./(sinc(x./(2*pi)) );
S1 = @(x) (1/2*1i) .*((P(x-w1)) - (P(x+w1)));
S2 = @(x) (1/2*1i) .*((P(x-w2)) - (P(x+w2)));
S = @(x) S1(x) + .8*S2(x);
f = abs(S(x)./S(w1));
freq = freqz(s(n),1,x);

figure;
plot(x/pi, f);
hold on;
xlim([-1,1]);
ylim([0,1]);
title('Problem 3.3 - Part 2');

error1 = norm(freq);
error2 = norm(S(x));
error = error1 - error2;
disp(error1);
disp(error2);
disp(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3.4

n = -5:1:45;
w1 = .2*pi;
w2 = .4*pi;
L = 40;
x = linspace(-pi,pi,1001);

p = @(n) double(n>=0 & n<=L-1);
P = @(x) L*exp( (-1i*(L-1)/2 ).*(x)).*sinc(x.*L./(2*pi))./(sinc(x./(2*pi)) );
s = @(n) (sin(w1*n) + .8.*(sin(w2*n))).*p(n);
S1 = @(x) (1/2*1i) .*((P(x-w1)) - (P(x+w1)));
S2 = @(x) (1/2*1i) .*((P(x-w2)) - (P(x+w2)));
S = @(x) S1(x) + .8*S2(x);
f = abs(S(x)./S(w1));

f2 = @(x) -abs(S(x));
peak = fminbnd(f2,0.1*pi,0.3*pi);
disp(peak);

peak = fminbnd(f2,0.3*pi,0.5*pi);
disp(peak);

L = 40;
w0 = 0.2*pi;
P = @(x) L*exp(-1i.*x*(L-1)/2).*(sinc(x.*L./(2*pi))./sinc(x./(2*pi)));
S = @(x) (1/(2*1i))*(P(x - w0) - P(x + w0));
f1 = @(x) -abs(S(x));
peak = fminbnd(f1,0.1*pi,0.3*pi);
disp(peak);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3.5

L = 80;
w1 = 0.1987*pi;
w2 = 0.4008*pi;
w0 = 0.1996*pi;
x = linspace(-pi,pi,1001);

P = @(x) L*exp(-1i.*x*(L-1)/2).*(sinc(x.*L./(2*pi))./sinc(x./(2*pi)));
S_P = @(x) (1/(2*1i))*(P(x - w0) - P(x + w0));
s = abs(S_P(x)/S_P(w0));
freq = linspace(-1,1,length(s));

figure;
plot(freq,s);
hold on;
xlim([-1,1]);
ylim([0,1]);
title("Single Sinusoid - 80");
hold off;

S_P1 = @(x) (1/(2*1i))*(P(x - w1) - P(x + w1));
S_P2 = @(x) (1/(2*1i))*(P(x - w2) - P(x + w2));
S = @(x) S_P1(x) + 0.8*S_P2(x);
sp = abs(S(x)/S(w1));

figure;
plot(freq,sp);
hold on;
xlim([-1,1]);
ylim([0,1]);
title('Double sinusoid - 80');
hold off;

L = 160;
w1 = 0.1997*pi;
w2 = 0.4002*pi;
w0 = 0.1999*pi;
x = linspace(-pi,pi,1001);

P = @(x) L*exp(-1i.*x*(L-1)/2).*(sinc(x.*L./(2*pi))./sinc(x./(2*pi)));
S_P = @(x) (1/(2*1i))*(P(x - w0) - P(x + w0));
sp = abs(S_P(x)/S_P(w0));
freq = linspace(-1,1,length(sp));

figure;
plot(freq,sp);
hold on;
xlim([-1,1]);
ylim([0,1]);
title('Single sinusoid - 160');
hold off;

S_P1 = @(x) (1/(2*1i))*(P(x - w1) - P(x + w1));
S_P2 = @(x) (1/(2*1i))*(P(x - w2) - P(x + w2));
S = @(x) S_P1(x) + 0.8*S_P2(x);
sp = abs(S(x)/S(w1));

figure;
plot(freq,sp);
hold on;
xlim([-1,1]);
ylim([0,1]);
title('Double sinusoid - 160');
hold off;
%}