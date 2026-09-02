function [LEs,meanEnergy,boundedFlag] = InitialConditionLLE(initXYZ)
p.a = 0.3;
p.b = 1.0;
p.c = 0.1;
p.apha = 0.1;
p.beta = 0.77;
p.A = 0.828;
p.f = 0.19;
p.k1 = 1.0;
p.k2 = 1.8;
p.u0 = 0;
h = 0.005;
m = 800000;
transient = 100000;
bound = 1e6;
n = size(initXYZ,1);
LEs = NaN(n,3);
meanEnergy = NaN(n,1);
boundedFlag = false(n,1);
parfor k = 1:n
    [ly,ok,H] = calcLE(initXYZ(k,:).',p,h,m,transient,bound);
    LEs(k,:) = ly(:).';
    meanEnergy(k) = H;
    boundedFlag(k) = ok;
end
end

function [ly,ok,Hmean] = calcLE(x,p,h,m,transient,bound)
Q = eye(3);
lySum = zeros(3,1);
Hsum = 0;
count = 0;
ok = true;
for i = 1:m
    t = (i-1)*h;
    [x,Phi] = rk4(x,t,h,p);
    if any(~isfinite(x)) || max(abs(x))>bound
        ly = NaN(3,1);
        Hmean = NaN;
        ok = false;
        return
    end
    if i > transient
        [Q,R] = qr(Phi*Q);
        d = abs(diag(R));
        d(d<realmin) = realmin;
        lySum = lySum+log(d);
        Hsum = Hsum+energy(x,t,p);
        count = count+1;
    end
end
ly = lySum/(count*h);
Hmean = Hsum/count;
end

function [xNext,Phi] = rk4(x,t,h,p)
I = eye(3);
k1 = rhs(x,t,p);
J1 = jac(x,t,p);
x2 = x+0.5*h*k1;
k2 = rhs(x2,t+0.5*h,p);
J2 = jac(x2,t+0.5*h,p);
x3 = x+0.5*h*k2;
k3 = rhs(x3,t+0.5*h,p);
J3 = jac(x3,t+0.5*h,p);
x4 = x+h*k3;
k4 = rhs(x4,t+h,p);
J4 = jac(x4,t+h,p);
xNext = x+h*(k1+2*k2+2*k3+k4)/6;
Phi = I+h*(J1+2*J2*(I+0.5*h*J1)+2*J3*(I+0.5*h*J2*(I+0.5*h*J1))+J4*(I+h*J3*(I+0.5*h*J2*(I+0.5*h*J1))))/6;
end

function H = energy(x,t,p)
u = p.u0+p.A*cos(2*pi*p.f*t);
H = 0.5*x(1)^2+0.5/p.c*x(2)^2+0.5*(p.apha*x(3)+p.beta*x(3)^3)*(u-x(1));
end

function dx = rhs(x,t,p)
u = p.u0+p.A*cos(2*pi*p.f*t);
dx = [-(p.apha+p.beta*x(3)^2)*(x(1)-u)-x(2)+x(1)-x(1)^3/3;p.c*(x(1)+p.a-p.b*x(2));-p.k1*(x(1)-u)-p.k2*x(3)];
end

function J = jac(x,t,p)
u = p.u0+p.A*cos(2*pi*p.f*t);
J = [1-x(1)^2-p.apha-p.beta*x(3)^2,-1,-2*p.beta*x(3)*(x(1)-u);p.c,-p.b*p.c,0;-p.k1,0,-p.k2];
end
