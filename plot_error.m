function result=plot_error
%PLOT_BF_ENRICHMENT Plot original saved BF and RFPS enrichment errors.
% Keep the original sample set, basis order and saved error normalization.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'results');
if ~exist(out,'dir'),mkdir(out);end
result=load(fullfile(root,'data','bf_enrichment_paper.mat'));
f=figure('visible','off','Color','w','position',[100 100 1200 450],'paperpositionmode','auto');
tags={'t1a','t64'};ms=[8 64];
for k=1:2
 base=result.(['err_bi_' tags{k}]);extra=result.(['err_LS_' tags{k}]);
 assert(numel(base)==15&&numel(extra)==18,'Unexpected saved basis counts');
 assert(max(abs(base(:)-extra(1:numel(base))))<1e-15,'Saved basis prefixes differ');
 subplot(1,2,k);
 h1=semilogy(1:numel(base),base,'ko-','LineWidth',1.5);hold on;
 h2=semilogy(1:numel(extra),extra,'bs--','LineWidth',1.5);
 assert(isequal(get(h1,'YData')(:),base(:))&&isequal(get(h2,'YData')(:),extra(:)));
 xlabel('Number of HF basis snapshots');ylabel('Saved scaled spatial error');
 title(sprintf('BF and RFPS enrichment (m=%d)',ms(k)));
 legend([h1 h2],{'Original BF (15 points)','BF + RFPS (18 points)'},'Location','northeast');grid on;
 fid=fopen(fullfile(out,sprintf('bf_enrichment_m%02d.csv',ms(k))),'w');
 assert(fid~=-1,'Cannot write CSV');fprintf(fid,'basis_count,original_BF,saved_RFPS_enrichment\n');
 for r=1:numel(extra)
  v=NaN;if r<=numel(base),v=base(r);end
  fprintf(fid,'%d,%.17g,%.17g\n',r,v,extra(r));
 end
 fclose(fid);
end
print(f,fullfile(out,'bf_enrichment_error.png'),'-dpng','-r180');close(f);
save(fullfile(out,'bf_enrichment_replayed.mat'),'-struct','result','-v7');
fprintf('Saved BF and enrichment errors replayed without recomputation.\n');
end
