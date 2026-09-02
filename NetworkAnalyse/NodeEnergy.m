function energy = NodeEnergy(X,Iext,params)
N = size(X,1);
if isscalar(Iext)
    Iext = Iext*ones(N,1);
else
    Iext = Iext(:);
end
x = X(:,1);
y = X(:,2);
z = X(:,3);
energy = 0.5*x.^2 + 0.5/params.c*y.^2 + 0.5*(params.apha*z + params.beta*z.^3).*(Iext-x);
end
