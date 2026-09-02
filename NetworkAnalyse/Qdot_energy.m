function dq = Qdot_energy(X,Iext,q,params)
q = q(:);
i = params.edge_i(:);
j = params.edge_j(:);
M = numel(i);
if numel(q) ~= M
    error('Qdot_energy: q length mismatch.')
end
type = params.couplingType;
if strcmpi(type,'constant')
    dq = zeros(M,1);
    return
end
q = max(q,getOpt(params,'qMin',1e-12));
if strcmpi(type,'energyExp')
    H = NodeEnergy(X,Iext,params);
    dq = params.rho_q*q.*double(abs(H(i)-H(j)) > params.epsE);
    return
end
if strcmpi(type,'hebbianPRX')
    N = params.N;
    qBase = getOpt(params,'qBase',1);
    pHebb = getOpt(params,'pHebb',2.5);
    kHebb = getOpt(params,'kHebb',params.rho_q);
    phi = tanh(X(:,1));
    dq = -(q-qBase)/pHebb + kHebb/(pHebb*N)*(phi(i).*phi(j));
    return
end
if strcmpi(type,'pddp')
    rho = getOpt(params,'rho_pddp',params.rho_q);
    amp = getOpt(params,'pddpAmp',0.9);
    theta = atan2(X(:,2)-params.phaseCenterY,X(:,1)-params.phaseCenterX);
    qTarget = 1 + amp*cos(theta(i)-theta(j));
    dq = rho*(qTarget-q);
    return
end
error('Qdot_energy: unknown coupling type %s.',type)
end

function v = getOpt(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
