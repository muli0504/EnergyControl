function dY = FHNfunction(Y,Iext,params)
N = params.N;
M = params.M_edge;
X = reshape(Y(1:3*N),N,3);
q = Y(3*N+1:3*N+M);
if isscalar(Iext)
    Iext = Iext*ones(N,1);
else
    Iext = Iext(:);
end
x = X(:,1);
y = X(:,2);
z = X(:,3);
xcouple = Icouple(X,q,params);
dx = (params.apha + params.beta*z.^2).*(Iext-x) - y + x - x.^3/3 + xcouple;
dy = params.c*(x + params.a - params.b*y);
dz = params.k1*(Iext-x) - params.k2*z;
dq = Qdot_energy(X,Iext,q,params);
dX = [dx,dy,dz];
dY = [dX(:);dq(:)];
end
