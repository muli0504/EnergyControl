function R = FHNdetailErrors_RMS(x0,params)
h = params.h;
T = params.T;
N = params.N;
M = params.M_edge;
Amp = params.Amp(:);
omega = 2*pi*params.freq(:);
recordStep = max(1,round(getOpt(params,'detailRecordStep',10)));
q0 = getOpt(params,'q0',ones(M,1));
Y = [x0(:);q0(:)];
t = 0;
nRec = floor(T/recordStep)+2;
tRec = NaN(nRec,1);
EH = NaN(nRec,1);
Ex = NaN(nRec,1);
Ey = NaN(nRec,1);
Ephi = NaN(nRec,1);
rec = 1;
[EH(rec),Ex(rec),Ey(rec),Ephi(rec)] = errors(reshape(Y(1:3*N),N,3),Amp.*cos(omega*t),params);
tRec(rec) = t;
rec = rec+1;
for step = 1:T
    I1 = Amp.*cos(omega*t);
    I2 = Amp.*cos(omega*(t+0.5*h));
    I3 = I2;
    I4 = Amp.*cos(omega*(t+h));
    K1 = FHNfunction(Y,I1,params);
    K2 = FHNfunction(Y+0.5*h*K1,I2,params);
    K3 = FHNfunction(Y+0.5*h*K2,I3,params);
    K4 = FHNfunction(Y+h*K3,I4,params);
    Y = Y+h*(K1+2*K2+2*K3+K4)/6;
    q = max(Y(3*N+1:3*N+M),getOpt(params,'qMin',1e-12));
    if isfield(params,'qMax') && ~isempty(params.qMax)
        q = min(q,params.qMax);
    end
    if getOpt(params,'rescaleQMean',false)
        q = q/mean(q);
    end
    Y(3*N+1:3*N+M) = q;
    t = t+h;
    if mod(step,recordStep)==0
        X = reshape(Y(1:3*N),N,3);
        [EH(rec),Ex(rec),Ey(rec),Ephi(rec)] = errors(X,I4,params);
        tRec(rec) = t;
        rec = rec+1;
    end
end
last = rec-1;
R.t = tRec(1:last);
R.EH = EH(1:last);
R.Ex = Ex(1:last);
R.Ey = Ey(1:last);
R.Ephi = Ephi(1:last);
end

function [EH,Ex,Ey,Ephi] = errors(X,Iext,params)
H = NodeEnergy(X,Iext,params);
x = X(:,1);
y = X(:,2);
z = X(:,3);
EH = sqrt(mean((H-mean(H)).^2));
Ex = sqrt(mean((x-mean(x)).^2));
Ey = sqrt(mean((y-mean(y)).^2));
Ephi = sqrt(mean((z-mean(z)).^2));
end

function v = getOpt(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
