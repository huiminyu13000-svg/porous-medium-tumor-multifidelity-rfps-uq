function result = run_experiment(mode,m,options)
%RUN_EXPERIMENT Generate inputs, solve PDEs, reconstruct, and plot errors.
% The primary error uses 500 independent parameters. run_reference uses 1200.
% All plotting routines are local functions at the end of this file.
if nargin<1,mode='fresh';end
if nargin<2,m=8;end
if nargin<3,options=struct();end
assert(any(strcmp(mode,{'fresh','short','reference'})),'Use fresh, short, or run_reference.');
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'solver'));
if strcmp(mode,'short') && ~isfield(options,'T'),options.T=0.03;end
cfg=experiment_config(m,options);
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
grid=cfg.domain(1):cfg.dx:cfg.domain(2);[X,Y]=meshgrid(grid);
training=sample_parameters(cfg.n_train,cfg.seed);
pool=sample_parameters(cfg.n_pool,cfg.seed+1);
evaluation=sample_parameters(cfg.n_eval,cfg.seed+2);
statistics=sample_parameters(cfg.n_stats,cfg.seed+3);
assert_disjoint(training,pool,evaluation,statistics);
signature=rmfield(cfg,{'n_eval','n_stats','make_plots','resume'});
basis_path=fullfile(out,sprintf('basis_m%02d.mat',m));
model=[];
if strcmp(mode,'reference') && exist(basis_path,'file')
    saved=load(basis_path);
    if isfield(saved,'model') && isequal(saved.model.signature,signature),model=saved.model;end
end
if isempty(model)
    model=build_model(cfg,training,pool,X,Y);
    model.signature=signature;
    atomic_save(basis_path,struct('model',model));
end
if strcmp(mode,'reference')
    result=compute_statistics(cfg,model,statistics,X,Y,out);
    return;
end
timer=tic;low=zeros(numel(X),cfg.n_eval);high=low;
for j=1:cfg.n_eval
    low(:,j)=solve_sample(evaluation(:,j),'ls',cfg);
    high(:,j)=solve_sample(evaluation(:,j),'hf',cfg);
    if mod(j,10)==0 || j==cfg.n_eval,fprintf('Error evaluation %d/%d\n',j,cfg.n_eval);end
end
evaluation_seconds=toc(timer);timer=tic;
result=struct('settings',cfg,'X',X,'Y',Y,'evaluation_parameters',evaluation, ...
    'basis_parameters',model.parameters,'environment',environment_info());
names={'bf','rfps','trf'};
for j=1:3
    basis=model.basis.(names{j});
    result.(names{j})=evaluate_basis(basis.low,basis.high,low,high);
end
result.common_budget=min([result.bf.n_basis,result.rfps.n_basis,result.trf.n_basis]);
r=result.common_budget;
result.equal_budget_error=[result.bf.mean_scaled_error(r),result.rfps.mean_scaled_error(r),result.trf.mean_scaled_error(r)];
result.construction_seconds=model.construction_seconds;
result.evaluation_seconds=evaluation_seconds;result.projection_seconds=toc(timer);
result.solver_calls=model.solver_calls;
result.solver_calls.evaluation_ls=cfg.n_eval;result.solver_calls.evaluation_hf=cfg.n_eval;
assert(all(isfinite(result.equal_budget_error)),'Nonfinite reconstruction errors.');
atomic_save(fullfile(out,sprintf('errors_m%02d.mat',m)),result);
if cfg.make_plots,plot_errors(result,out);end
fprintf('m=%d, T=%g, N_eval=%d, common basis=%d\n',m,cfg.T,cfg.n_eval,r);
fprintf('Mean scaled spatial error (BF/RFPS/TRF): %.7g %.7g %.7g\n',result.equal_budget_error);
end

function model=build_model(cfg,training,pool,X,Y)
timer=tic;M=numel(X);candidate_low=zeros(M,cfg.n_train);
pool_low=zeros(M,cfg.n_pool);coarse=pool_low;
candidate_high=nan(M,cfg.n_train);pool_high=nan(M,cfg.n_pool);
counts=struct('training_ls',0,'pool_ls',0,'coarse_pme',0,'training_hf',0,'pool_hf',0);
for j=1:cfg.n_train
    candidate_low(:,j)=solve_sample(training(:,j),'ls',cfg);counts.training_ls=counts.training_ls+1;
    if mod(j,25)==0 || j==cfg.n_train,fprintf('Training LS %d/%d\n',j,cfg.n_train);end
end
for j=1:cfg.n_pool
    pool_low(:,j)=solve_sample(pool(:,j),'ls',cfg);counts.pool_ls=counts.pool_ls+1;
    coarse(:,j)=solve_sample(pool(:,j),'coarse',cfg);counts.coarse_pme=counts.coarse_pme+1;
    if mod(j,20)==0 || j==cfg.n_pool,fprintf('Supplementary pool %d/%d\n',j,cfg.n_pool);end
end
[base_ids,~]=select_basis(candidate_low,training,cfg.initial_basis+cfg.max_add);
assert(~isempty(base_ids),'Training fields have zero numerical rank.');
seed_ids=base_ids(1:min(cfg.initial_basis,numel(base_ids)));
seed_low=candidate_low(:,seed_ids);seed_high=zeros(size(seed_low));
for j=1:numel(seed_ids),seed_high(:,j)=candidate_hf(seed_ids(j));end
% Analytic form of the original spatial weights, generated without a MAT file.
weights=exp(-(X(:).^2+Y(:).^2)/(2*cfg.rfps_spatial_sigma^2));
[rfps_low,rfps_high,rfps_added]=enrich_basis(seed_low,seed_high,pool_low,@pool_hf, ...
    cfg.max_add,cfg.rfps_alpha,true,weights,cfg.rfps_weight,[]);
[tri_ids,~]=select_basis(coarse,pool,cfg.initial_basis);
assert(~isempty(tri_ids),'Coarse training fields have zero numerical rank.');
tri_high=zeros(M,numel(tri_ids));
for j=1:numel(tri_ids),tri_high(:,j)=pool_hf(tri_ids(j));end
[trf_low,trf_high,trf_added]=enrich_basis(pool_low(:,tri_ids),tri_high,pool_low,@pool_hf, ...
    cfg.max_add,cfg.rfps_alpha,true,weights,cfg.rfps_weight,[]);
base_low=candidate_low(:,base_ids);base_high=zeros(size(base_low));
for j=1:numel(base_ids),base_high(:,j)=candidate_hf(base_ids(j));end
model.basis.bf=struct('low',base_low,'high',base_high,'map',projection_matrix(base_low));
model.basis.rfps=struct('low',rfps_low,'high',rfps_high,'map',projection_matrix(rfps_low));
model.basis.trf=struct('low',trf_low,'high',trf_high,'map',projection_matrix(trf_low));
model.parameters=struct('bf',training(:,base_ids),'rfps',[training(:,seed_ids),pool(:,rfps_added)], ...
    'trf',[pool(:,tri_ids),pool(:,trf_added)]);
model.indices=struct('bf',base_ids,'rfps_seed',seed_ids,'rfps_added',rfps_added, ...
    'trf_seed',tri_ids,'trf_added',trf_added);
model.training_parameters=training;model.pool_parameters=pool;
model.construction_seconds=toc(timer);model.solver_calls=counts;

    function value=candidate_hf(index)
        if any(isnan(candidate_high(:,index)))
            candidate_high(:,index)=solve_sample(training(:,index),'hf',cfg);
            counts.training_hf=counts.training_hf+1;
        end
        value=candidate_high(:,index);
    end
    function value=pool_hf(index)
        if any(isnan(pool_high(:,index)))
            pool_high(:,index)=solve_sample(pool(:,index),'hf',cfg);
            counts.pool_hf=counts.pool_hf+1;
        end
        value=pool_high(:,index);
    end
end

function result=compute_statistics(cfg,model,z,X,Y,out)
% Online mean/M2 updates avoid storing 1200 full HF fields in memory.
names={'bf','rfps','trf','hf'};M=numel(X);
key=struct('settings',cfg,'basis_parameters',model.parameters,'parameters',z);
checkpoint=fullfile(out,sprintf('reference_checkpoint_m%02d.mat',cfg.m));
state=struct('completed',0,'mean',zeros(M,4),'M2',zeros(M,4),'seconds',0);
if cfg.resume && exist(checkpoint,'file')
    saved=load(checkpoint);
    if isfield(saved,'key') && isequal(saved.key,key),state=saved.state;end
end
start=state.completed;timer=tic;
for j=start+1:cfg.n_stats
    low=solve_sample(z(:,j),'ls',cfg);high=solve_sample(z(:,j),'hf',cfg);
    values=zeros(M,4);values(:,4)=high;
    for k=1:3,b=model.basis.(names{k});values(:,k)=b.high*(b.map*low);end
    delta=values-state.mean;state.mean=state.mean+delta/j;
    state.M2=state.M2+delta.*(values-state.mean);state.completed=j;
    if mod(j,10)==0 || j==cfg.n_stats
        state.seconds=state.seconds+toc(timer);
        atomic_save(checkpoint,struct('key',key,'state',state));timer=tic;
        fprintf('Independent statistics %d/%d\n',j,cfg.n_stats);
    end
end
result=struct('settings',cfg,'X',X,'Y',Y,'statistics_parameters',z, ...
    'sample_count',state.completed,'new_hf_calls',state.completed-start, ...
    'new_ls_calls',state.completed-start,'statistics_seconds',state.seconds, ...
    'basis_parameters',model.parameters,'environment',environment_info());
for k=1:4
    result.(names{k})=struct('mean',state.mean(:,k),'std',sqrt(max(0,state.M2(:,k)/cfg.n_stats)));
    if k<=3,result.(names{k}).n_basis=size(model.basis.(names{k}).low,2);end
end
assert(all(isfinite([state.mean(:);state.M2(:)])),'Nonfinite reference statistics.');
atomic_save(fullfile(out,sprintf('reference_m%02d.mat',cfg.m)),result);
if cfg.make_plots,plot_statistics(result,out);end
fprintf('m=%d, T=%g: mean and standard deviation from %d independent samples.\n',cfg.m,cfg.T,cfg.n_stats);
end

function assert_disjoint(varargin)
for a=1:nargin
    for b=a+1:nargin
        assert(~any(ismember(varargin{a}',varargin{b}','rows')),'Parameter sets overlap.');
    end
end
end

function info=environment_info()
info=struct('version',version,'computer',computer);
end

function atomic_save(path,content)
temporary=[tempname(fileparts(path)) '.mat'];save(temporary,'-struct','content','-v7');
movefile(temporary,path,'f');
end

function plot_errors(R,out)
names={'bf','rfps','trf'};labels={'BF','BF + RFPS','TRF'};styles={'ko-','bs--','rd-.'};
f=figure('Color','w','Position',[100 100 780 500]);hold on;
for j=1:3
    values=R.(names{j}).mean_scaled_error;
    semilogy(1:numel(values),max(values,realmin),styles{j},'LineWidth',1.5,'MarkerSize',4);
end
set(gca,'YScale','log');grid on;xlabel('Number of HF basis functions');ylabel('Mean scaled spatial error');
title(sprintf('m=%d, T=%g, %d independent evaluation samples',R.settings.m,R.settings.T,R.settings.n_eval));
legend(labels,'Location','best');print(f,fullfile(out,sprintf('errors_m%02d.png',R.settings.m)),'-dpng','-r160');
end

function plot_statistics(R,out)
names={'bf','rfps','trf','hf'};labels={'BF','BF + RFPS','TRF','HF reference'};
[~,row]=min(abs(R.Y(:,1)));
for item={'mean','std'}
    stat=item{1};reference=reshape(R.hf.(stat),size(R.X));
    fields=cell(1,4);for j=1:4,fields{j}=reshape(R.(names{j}).(stat),size(R.X));end
    limits=[min(cellfun(@(v) min(v(:)),fields)),max(cellfun(@(v) max(v(:)),fields))];
    if limits(2)<=limits(1),limits(2)=limits(1)+1;end
    f=figure('Color','w','Position',[100 100 1400 660]);
    for j=1:4
        subplot(2,4,j);imagesc(R.X(1,:),R.Y(:,1),fields{j});axis image;set(gca,'YDir','normal');
        caxis(limits);colorbar;xlabel('x');ylabel('y');title(sprintf('%s: %s',labels{j},stat));
        subplot(2,4,j+4);plot(R.X(1,:),fields{j}(row,:),'k-',R.X(1,:),reference(row,:),'r--','LineWidth',1.4);
        xlabel('x');ylabel(stat);grid on;title(sprintf('y=0, N=%d',R.sample_count));
    end
    print(f,fullfile(out,sprintf('%s_m%02d.png',stat,R.settings.m)),'-dpng','-r160');
end
end
