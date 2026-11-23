function result = evaluate_basis(low_basis,high_basis,low_eval,high_eval)
%EVALUATE_BASIS Mean scaled spatial error used in the presentation.
% E = mean_i(norm(H_i - approximation_i,2)) / (Nx*Ny).
assert(size(low_basis,2)==size(high_basis,2),'Paired basis mismatch.');
assert(isequal(size(low_eval),size(high_eval)),'Evaluation matrix mismatch.');
assert(all(isfinite([low_basis(:);high_basis(:);low_eval(:);high_eval(:)])),'Nonfinite input.');
M=size(high_eval,1); rmax=size(low_basis,2);
result.mean_scaled_error=zeros(rmax,1);
result.rmse=zeros(rmax,1);result.relative_l2=zeros(rmax,1);
for r=1:rmax
    reconstruction=high_basis(:,1:r)*(projection_matrix(low_basis(:,1:r))*low_eval);
    difference=high_eval-reconstruction;
    result.mean_scaled_error(r)=mean(sqrt(sum(difference.^2,1)))/M;
    result.rmse(r)=sqrt(mean(difference(:).^2));
    result.relative_l2(r)=norm(difference,'fro')/max(norm(high_eval,'fro'),eps);
end
result.n_basis=rmax;result.n_eval=size(high_eval,2);
end
