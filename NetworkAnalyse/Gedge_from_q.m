function g = Gedge_from_q(q,params)
q = q(:);
M = numel(q);
if strcmpi(params.couplingType,'constant')
    g = params.gc*ones(M,1);
    return
end
q = max(q,getOpt(params,'qMin',1e-12));
if isfield(params,'qMax') && ~isempty(params.qMax)
    q = min(q,params.qMax);
end
qMean = mean(q);
if ~isfinite(qMean) || qMean <= 0
    error('Gedge_from_q: invalid q mean.')
end
g = params.gc*q/qMean;
end

function v = getOpt(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end
