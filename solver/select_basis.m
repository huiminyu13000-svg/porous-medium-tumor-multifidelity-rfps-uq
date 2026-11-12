% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Original author solver; executable formulas retained.
% Xu'an Dou attribution is retained where present in the original source.
function[p,L] = select_basis(V, z, S)
% pchol -- Full-pivoting Cholesky factorization
%
% [p,L] = cholesky_node_selection_ip(V, z, S)
%     V: is the low fidelity data
%     Computes truncated Cholesky factorization: Given a function handle ipfun
%     and a (P x N) matrix z, this function computes a greedy-S cholesky
%     factorization of the Gramian G whose entries are 
%
%       G_{m,n} = ipfun(z(:,m), z(:,n))
%
%     The truncation is taken by performing a pivoted Cholesky decomposition on
%     G, and then stopping after S steps. The output is an S x S Cholesky
%     factor L and a permutation vector p such that L*L' = G(p,p).
%
%     This function is, in principle, faster than a direct Cholesky
%     factorization because it only computes the necessary inner products
%     rather than forming the full Gramian G.
rank_tol = 1e-8;
N = size(V, 2);
if isempty(z),z=zeros(0,N);end
assert(size(z,2)==N,'Snapshot and parameter columns must match');
S=min(S,N);
L = zeros([N S]);
Lsum = zeros([N 1]);
global_indices = (1:N).';
w = zeros([N 1]);
w = sum(V.*V);
w = w(:);
for n = 1:S
  % Find the next element
  [garbage, p] = max(w(n:end));
  p = p + n - 1;
  % In the case of singularity, break.
  if garbage < rank_tol
    warning(['Gramian numerically singular...stopping at ' num2str(n-1) ' instead of ' num2str(S) ' nodes']);
    n = n-1;
    break
 end
  % Exchange next element to the pivot position
  w([n p]) = w([p n]);
  L([n p],:) = L([p n],:);
  z(:, [n p]) = z(:, [p n]);
  % Change col of the data
  V(:, [n p]) = V(:, [p n]);
  % Update permutation vector
  global_indices([n p]) = global_indices([p n]);
  % Form inner product and eliminate from previous steps
  r = V(:,n+1:end)'*V(:,n); %ipfun(z(:,n+1:end), z(:,n));
  r = r(:);
  r = r - L(n+1:end,1:n-1)*L(n,1:n-1)';
  % Store information from this stage
  L(n,n) = sqrt(w(n));
  L(n+1:end,n) = r./L(n,n);
  w(n+1:end) = w(n+1:end) - L(n+1:end,n).^2;
end
% Truncate to appropriate size
L = L(1:n,1:n);
p = global_indices(1:n);
