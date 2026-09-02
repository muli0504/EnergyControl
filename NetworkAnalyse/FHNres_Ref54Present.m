function R = FHNres_Ref54Present(x0,params)
h = params.h;
T = params.T;
N = params.N;
M = params.M_edge;
Amp = params.Amp(:);
omega = 2*pi*params.freq(:);
transient = params.transient;
K = T-transient;
mode = lower(params.mode);
if strcmp(mode,'present')
    Y = [x0(:);params.q0(:)];
elseif strcmp(mode,'ref54')
    Y = [x0(:);params.k0(:)];
else
    error('FHNres_Ref54Present: invalid mode.')
end
S = initStat(N);
HJsum = 0;
gsum = 0;
valid = true;
t = 0;
for step = 1:T
    I1 = Amp.*cos(omega*t);
    I2 = Amp.*cos(omega*(t+0.5*h));
    I3 = I2;
    I4 = Amp.*cos(omega*(t+h));
    K1 = rhs(Y,I1,params);
    K2 = rhs(Y+0.5*h*K1,I2,params);
    K3 = rhs(Y+0.5*h*K2,I3,params);
    K4 = rhs(Y+h*K3,I4,params);
    Y = Y+h*(K1+2*K2+2*K3+K4)/6;
    if any(~isfinite(Y))
        valid = false;
        break
    end
    X = reshape(Y(1:3*N),N,3);
    if max(abs(X),[],'all') > params.maxAbsState
        valid = false;
        break
    end
    if strcmp(mode,'present')
        q = max(Y(3*N+1:3*N+M),params.qMin);
        if any(~isfinite(q)) || max(q)>params.maxAdaptive
            valid = false;
            break
        end
        Y(3*N+1:3*N+M) = q;
    else
        k = max(Y(3*N+1:3*N+N),params.kMin);
        if any(~isfinite(k)) || max(k)>params.maxAdaptive
            valid = false;
            break
        end
        Y(3*N+1:3*N+N) = k;
    end
    t = t+h;
    if step <= transient
        continue
    end
    H = NodeEnergy(X,I4,params);
    S = updateStat(S,H,K);
    x = X(:,1);
    if strcmp(mode,'present')
        q = Y(3*N+1:3*N+M);
        g = params.gc*q/mean(q);
    else
        k = Y(3*N+1:3*N+N);
        g = 0.5*(k(params.edge_i)+k(params.edge_j));
    end
    dx = x(params.edge_i)-x(params.edge_j);
    HJsum = HJsum+mean(g.*dx.^2);
    gsum = gsum+mean(g);
end
R.valid = valid;
if valid
    R.hSF = syncFactor(S);
    R.HJ = HJsum/K;
    R.meanCoupling = gsum/K;
else
    R.hSF = NaN;
    R.HJ = NaN;
    R.meanCoupling = NaN;
end
end

function dY = rhs(Y,Iext,params)
N = params.N;
M = params.M_edge;
X = reshape(Y(1:3*N),N,3);
x = X(:,1);
y = X(:,2);
z = X(:,3);
H = NodeEnergy(X,Iext,params);
if strcmpi(params.mode,'present')
    q = Y(3*N+1:3*N+M);
    active = double(abs(H(params.edge_i)-H(params.edge_j))>=params.epsH);
    dq = params.sigmaPresent*q.*active;
    g = params.gc*q/mean(q);
    dxEdge = params.B*x;
    coup = -params.B.'*(g.*dxEdge);
    dAdaptive = dq;
else
    k = Y(3*N+1:3*N+N);
    dH = abs(H(params.edge_i)-H(params.edge_j));
    mismatch = accumarray(params.edge_i,dH,[N,1],@sum,0)+accumarray(params.edge_j,dH,[N,1],@sum,0);
    dk = params.sigmaRef54*k.*double(mismatch>=params.epsH);
    coup = k.*(params.Adj*x-params.deg.*x);
    dAdaptive = dk;
end
dx = (params.apha+params.beta*z.^2).*(Iext-x)-y+x-x.^3/3+coup;
dy = params.c*(x+params.a-params.b*y);
dz = params.k1*(Iext-x)-params.k2*z;
dX = [dx,dy,dz];
dY = [dX(:);dAdaptive(:)];
end

function S = initStat(N)
S.s1 = 0;
S.s2 = 0;
S.s3 = zeros(N,1);
S.s4 = zeros(N,1);
end

function S = updateStat(S,x,K)
F = mean(x);
S.s1 = S.s1+F^2/K;
S.s2 = S.s2+F/K;
S.s3 = S.s3+x.^2/K;
S.s4 = S.s4+x/K;
end

function sf = syncFactor(S)
den = mean(S.s3-S.s4.^2);
if den > 0
    sf = (S.s1-S.s2^2)/den;
else
    sf = NaN;
end
end
