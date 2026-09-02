function [xc,yc] = PhaseCenter(params,Amp,freq,numStep,transientRatio)
h = params.h;
omega = 2*pi*freq;
x = 0.2;
y = 0.1;
z = 0.05;
transientStep = round(transientRatio*numStep);
sumX = 0;
sumY = 0;
count = 0;
for step = 1:numStep
    t = (step-1)*h;
    I1 = Amp*cos(omega*t);
    dx1 = (params.apha+params.beta*z^2)*(I1-x)-y+x-x^3/3;
    dy1 = params.c*(x+params.a-params.b*y);
    dz1 = params.k1*(I1-x)-params.k2*z;
    x2 = x+0.5*h*dx1;
    y2 = y+0.5*h*dy1;
    z2 = z+0.5*h*dz1;
    I2 = Amp*cos(omega*(t+0.5*h));
    dx2 = (params.apha+params.beta*z2^2)*(I2-x2)-y2+x2-x2^3/3;
    dy2 = params.c*(x2+params.a-params.b*y2);
    dz2 = params.k1*(I2-x2)-params.k2*z2;
    x3 = x+0.5*h*dx2;
    y3 = y+0.5*h*dy2;
    z3 = z+0.5*h*dz2;
    dx3 = (params.apha+params.beta*z3^2)*(I2-x3)-y3+x3-x3^3/3;
    dy3 = params.c*(x3+params.a-params.b*y3);
    dz3 = params.k1*(I2-x3)-params.k2*z3;
    x4 = x+h*dx3;
    y4 = y+h*dy3;
    z4 = z+h*dz3;
    I4 = Amp*cos(omega*(t+h));
    dx4 = (params.apha+params.beta*z4^2)*(I4-x4)-y4+x4-x4^3/3;
    dy4 = params.c*(x4+params.a-params.b*y4);
    dz4 = params.k1*(I4-x4)-params.k2*z4;
    x = x+h*(dx1+2*dx2+2*dx3+dx4)/6;
    y = y+h*(dy1+2*dy2+2*dy3+dy4)/6;
    z = z+h*(dz1+2*dz2+2*dz3+dz4)/6;
    if step > transientStep
        sumX = sumX+x;
        sumY = sumY+y;
        count = count+1;
    end
end
if count <= 0
    error('PhaseCenter: no valid data.')
end
xc = sumX/count;
yc = sumY/count;
end
