function coup = Icouple(X,q,params)
x = X(:,1);
g = Gedge_from_q(q,params);
if isfield(params,'B') && ~isempty(params.B)
    dx = params.B*x;
    coup = -params.B.'*(g.*dx);
else
    i = params.edge_i(:);
    j = params.edge_j(:);
    dx = x(j)-x(i);
    ci = g.*dx;
    N = size(X,1);
    coup = accumarray(i,ci,[N,1],@sum,0) + accumarray(j,-ci,[N,1],@sum,0);
end
end
