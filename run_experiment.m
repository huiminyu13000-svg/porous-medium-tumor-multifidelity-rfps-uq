% 完整实验入口：读取原设置、调用求解器、重构并保存结果。
function results=run_experiment(mode,m)
%RUN_EXPERIMENT Recompute original BF/RFPS/TRF algorithms, retaining PDE solvers.
% cached: redo selection, enrichment, projection and errors using saved snapshots.
% pde_check: run one complete-time PME/LS/coarse-PME sample at saved settings.
% fresh: regenerate candidate LS and selected HF; reuse expensive HF reference.
% fresh_reference: additionally regenerate the complete held-out HF reference.
if nargin<1,mode='fresh';end
if nargin<2,m=8;end
m=double(m);
assert(ismember(m,[8 64]),'m must be 8 or 64');
assert(any(strcmp(mode,{'cached','short','pde_check','fresh','fresh_reference'})),'Unknown mode');
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'solver'));out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
S=load(fullfile(root,'data',sprintf('bf_compute_m%02d.mat',m)));
% MAT metadata may use integer storage; time arithmetic must use floating point.
S.T_final=double(S.T_final);S.dx=double(S.dx);S.dt_PME=double(S.dt_PME);S.dt_HS=double(S.dt_HS);
if strcmp(mode,'short')
 % 短时测试使用原网格和原步长：LS 1 步，PME/粗 PME 各 25 步。
 Tshort=.05;j=1;p=[S.lambda_t(j);S.cB_t(j);S.G0_t(j);S.R0_t(j)];
 initial=.95*double(S.X.^2+S.Y.^2<p(4)^2);
 L=solve_heleshaw(j,Tshort,S.dx,S.dx,S.dt_HS,p(1),p(2),p(3),initial);
 H=solve_pme(j,Tshort,m,S.dx,S.dx,S.dt_PME,p(1),p(2),p(3),initial);
 [Xc,Yc]=meshgrid(-2.5:.1:2.5);initial_c=.95*double(Xc.^2+Yc.^2<p(4)^2);
 C=solve_cmpme(j,Tshort,8,.1,.002,p(1),p(2),p(3),initial_c);
 assert(isequal(size(L),[101 101])&&isequal(size(H),[101 101])&&isequal(size(C),[51 51]));
 assert(all(isfinite([L(:);H(:);C(:)])),'Nonfinite short-run output');
 results=struct('T',Tshort,'m',m,'L',L,'H',H,'C',C,'passed',true);
 save(fullfile(out,sprintf('short_m%02d.mat',m)),'-struct','results','-v7');
 fprintf('PASS: short PDE m=%d, T=%g\n',m,Tshort);return;
end
if strcmp(mode,'pde_check')
 j=1;p=[S.lambda_t(j);S.cB_t(j);S.G0_t(j);S.R0_t(j)];initial=.95*double(S.X.^2+S.Y.^2<p(4)^2);
 L=solve_heleshaw(j,S.T_final,S.dx,S.dx,S.dt_HS,p(1),p(2),p(3),initial);
 H=solve_pme(j,S.T_final,m,S.dx,S.dx,S.dt_PME,p(1),p(2),p(3),initial);
 [Xc,Yc]=meshgrid(-2.5:.1:2.5);initial_c=.95*double(Xc.^2+Yc.^2<p(4)^2);
 C=solve_cmpme(j,S.T_final,8,.1,.002,p(1),p(2),p(3),initial_c);
 results=struct('L',L,'H',H,'C',C,'m',m,'sample',j,'T',S.T_final,'PME_steps',round(S.T_final/S.dt_PME), ...
 'LS_saved_max_difference',max(abs(L(:)-double(S.rhoL_test(:,j)))), ...
 'HF_saved_max_difference',max(abs(H(:)-S.rhoH_test(:,j))));
 Cfine=interp2(Xc,Yc,C,S.X,S.Y,'linear',0);
 results.coarse_saved_max_difference=max(abs(Cfine(:)-S.tri_rhoM_(:,j)));
 save(fullfile(out,sprintf('pde_check_m%02d.mat',m)),'-struct','results','-v7');
 fprintf('PDE m=%d T=%g: saved LS/HF max difference %.4g / %.4g\n',m,S.T_final,results.LS_saved_max_difference,results.HF_saved_max_difference);return;
end
candidate=double(S.rho_L);L=double(S.rhoL_test);H=double(S.rhoH_test);
if any(strcmp(mode,{'fresh','fresh_reference'}))
 for j=1:size(S.z,2)
  p=[S.lambda_samples(j);S.cB_samples(j);S.G0_samples(j);S.R0_samples(j)];
  initial=.95*double(S.X.^2+S.Y.^2<p(4)^2);
  a=solve_heleshaw(j,S.T_final,S.dx,S.dx,S.dt_HS,p(1),p(2),p(3),initial);candidate(:,j)=a(:);
  fprintf('Candidate LS %d/%d\n',j,size(candidate,2));
 end
 for j=1:size(S.z_test,2)
  p=[S.lambda_t(j);S.cB_t(j);S.G0_t(j);S.R0_t(j)];initial=.95*double(S.X.^2+S.Y.^2<p(4)^2);
  a=solve_heleshaw(j,S.T_final,S.dx,S.dx,S.dt_HS,p(1),p(2),p(3),initial);L(:,j)=a(:);
  if strcmp(mode,'fresh_reference'),a=solve_pme(j,S.T_final,m,S.dx,S.dx,S.dt_PME,p(1),p(2),p(3),initial);H(:,j)=a(:);end
 end
end
[ids,~]=select_basis(candidate,S.z,15);ids=ids(:)';
PL=candidate(:,ids);
if strcmp(mode,'cached')
 [found,loc]=ismember(ids,S.u_inp(:)');
 assert(all(found),'Recomputed selection contains HF nodes absent from cache; use fresh mode');
 PH=S.phi_H(:,loc);
else
 PH=zeros(numel(S.X),numel(ids));
 for j=1:numel(ids)
  k=ids(j);p=[S.lambda_samples(k);S.cB_samples(k);S.G0_samples(k);S.R0_samples(k)];
  initial=.95*double(S.X.^2+S.Y.^2<p(4)^2);
  a=solve_pme(k,S.T_final,m,S.dx,S.dx,S.dt_PME,p(1),p(2),p(3),initial);PH(:,j)=a(:);
 end
end
[PE,QE,added]=enrich_basis(PL,PH,L,H,S.tri_K_add,S.tri_alpha,true,S.tri_wvec,S.tri_lambda,[]);
[base,base_reconstruction]=projection_errors(PL,PH,L,H);
[enriched,enriched_reconstruction]=projection_errors(PE,QE,L,H);
% Rebuild the archived TRF seed from the original coarse snapshot pool.
C=double(S.tri_rhoM_);TL=double(S.tri_rhoL_);TH=double(S.tri_rhoH_);
if any(strcmp(mode,{'fresh','fresh_reference'}))
 TL=L;TH=H;[Xc,Yc]=meshgrid(-2.5:.1:2.5);
 for j=1:size(S.z_test,2)
  p=[S.lambda_t(j);S.cB_t(j);S.G0_t(j);S.R0_t(j)];initial_c=.95*double(Xc.^2+Yc.^2<p(4)^2);
  a=solve_cmpme(j,S.T_final,8,.1,.002,p(1),p(2),p(3),initial_c);a=interp2(Xc,Yc,a,S.X,S.Y,'linear',0);C(:,j)=a(:);
  fprintf('TRF coarse PME %d/%d\n',j,size(C,2));
 end
end
[coarse_ids,~]=select_basis(C,S.z_test,15);coarse_ids=coarse_ids(:)';
[TP,TQ,trf_added]=enrich_basis(TL(:,coarse_ids),TH(:,coarse_ids),TL,TH,S.tri_K_add,S.tri_alpha,true,S.tri_wvec,S.tri_lambda,[]);
[trf_error,trf_reconstruction]=projection_errors(TP,TQ,TL,TH);
results=struct('mode',mode,'m',m,'initial_indices',ids,'enrichment_indices',added, ...
 'trf_initial_indices',coarse_ids,'trf_enrichment_indices',trf_added, ...
 'baseline_error',base,'enriched_error',enriched,'trf_error',trf_error, ...
 'baseline_reconstruction',base_reconstruction,'enriched_reconstruction',enriched_reconstruction, ...
 'trf_reconstruction',trf_reconstruction,'computed_phiL',PE,'computed_phiH',QE, ...
 'saved_initial_order_matches',isequal(ids,S.u_inp(:)'), ...
 'saved_enriched_basis_max_difference',array_difference(PE,double(S.tri_phiL_ls)), ...
 'saved_enriched_error_max_difference',array_difference(enriched(:),S.tri_err_LS(:)), ...
 'saved_trf_error_max_difference',array_difference(trf_error(:),S.tri_err_TRI(:)));
assert(all(isfinite(enriched))&&all(isfinite(trf_error)),'Nonfinite reconstruction error');
save(fullfile(out,sprintf('computed_%s_m%02d.mat',mode,m)),'-struct','results','-v7');
f=figure('visible','off','Color','w','paperpositionmode','auto');
semilogy(1:numel(base),base,'ko-',1:numel(enriched),enriched,'bs--','LineWidth',1.5);
xlabel('HF basis count');ylabel('Original scaled spatial error');legend({'Computed BF','Computed BF + RFPS'});grid on;
print(f,fullfile(out,sprintf('computed_%s_m%02d.png',mode,m)),'-dpng','-r150');close(f);
fprintf('Computed m=%d: saved enriched basis/error difference %.4g / %.4g; TRF error difference %.4g\n',m,results.saved_enriched_basis_max_difference,results.saved_enriched_error_max_difference,results.saved_trf_error_max_difference);
end
function [errors,reconstruction]=projection_errors(PL,PH,L,H)
% Original conditional regularization and original sample/error normalization.
errors=zeros(size(PL,2),1);
for r=1:size(PL,2)
 G=PL(:,1:r)'*PL(:,1:r);rhs=PL(:,1:r)'*L;
 if r>1&&cond(G)>1e12,G=G+1e-8*trace(G)/r*eye(r);end
 reconstruction=PH(:,1:r)*(G\rhs);D=H-reconstruction;
 errors(r)=sqrt(mean(sum(D.^2,1)/size(H,1)))/size(H,2);
end
end

function value=array_difference(A,B)
% Different accepted enrichment counts are reported, not silently padded.
if isequal(size(A),size(B)),value=max(abs(A(:)-B(:)));else,value=NaN;end
end
