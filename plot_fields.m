function result=plot_fields(make_plots)
%RUN_PAPER_FIGURES Replay numeric fields from the supplied original TRF workspace.
% No new sample selection, PDE solve, smoothing or statistical recomputation.
% The numeric cache avoids loading MATLAB graphics objects in the raw workspace.
if nargin<1,make_plots=true;end
root=fileparts(mfilename('fullpath'));
result=load(fullfile(root,'data','trifidelity_paper_fields.mat'));
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
names={'CH','LS','TRI','H'};labels={'CMPME (BF)','LS (BF)','TRF','HF'};
tags={'t1a','t64'};ms=[8,64];
for stat={'mean','std'}
 statistic=stat{1};
 for i=1:2
  for j=1:4
   field=['rho' names{j} '_' statistic '_' tags{i}];
   assert(isfield(result,field),'Missing original field: %s',field);
   values=result.(field);
   assert(isequal(size(values),size(result.X))&&all(isfinite(values(:))),'Invalid original field: %s',field);
  end
 end
 if ~make_plots,continue;end
 f=figure('Color','w','visible','off','position',[100 100 1600 750],'paperpositionmode','auto');
 for i=1:2
  for j=1:3
   ids=[1 3 4];id=ids(j);field=['rho' names{id} '_' statistic '_' tags{i}];
   subplot(2,4,(i-1)*4+j);
   imagesc(result.X(1,:),result.Y(:,1),result.(field));axis image;colorbar;
   title(sprintf('%s %s (m=%d)',statistic,labels{id},ms(i)));
   xlabel('x');ylabel('y');
  end
  subplot(2,4,i*4);hold on;styles={'k-','g-','b-','r--'};
  mid=round(size(result.X,1)/2);h=zeros(1,4);
  for j=1:4
   field=['rho' names{j} '_' statistic '_' tags{i}];values=result.(field);
   h(j)=plot(result.X(1,:),values(mid,:),styles{j},'LineWidth',1.5);
  end
  title(sprintf('%s slice (m=%d)',statistic,ms(i)));xlabel('x');grid on;
  legend(h,labels,'Location','northwest');axis tight;
 end
 print(f,fullfile(out,['paper_trifidelity_' statistic '.png']),'-dpng','-r150');close(f);
end
save(fullfile(out,'paper_trifidelity_replayed_fields.mat'),'-struct','result','-v7');
fprintf('Original saved fields replayed without numerical changes: %s\n',out);
end
