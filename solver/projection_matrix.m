function map = projection_matrix(low_basis)
%PROJECTION_MATRIX Stable map from LF fields to paired basis coefficients.
gram=low_basis'*low_basis;r=size(gram,1);
assert(r>0,'The selected basis is empty.');
if rcond(gram)<1e-12,gram=gram+1e-8*max(trace(gram)/r,eps)*eye(r);end
map=gram\low_basis';
end
