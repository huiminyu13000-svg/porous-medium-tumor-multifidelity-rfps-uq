function [rho,p,info] = solve_heleshaw(k,T,dx,dy,dt,lambda,cB,G0_k,rho0)
%SOLVE_HELESHAW Pressure solve, fast-marching extension and Godunov update.
% The input multiplier is converted once to the physical coefficient.
G0=0.5*G0_k; nominal_dt=dt;
assert(T>0&&dt>0,'Time and step must be positive.');
phi=0.5*min(dx,dy)*(1-2*double(rho0>0));
[phi,~]=fast_marching(phi,zeros(size(phi)),dx,dy);
Nt=ceil(T/nominal_dt-1e-12);p=zeros(size(phi));
info=struct('steps',Nt,'final_time',T,'last_dt',0,'max_courant',0);
for step=1:Nt
    dt=min(nominal_dt,T-(step-1)*nominal_dt);info.last_dt=dt;
    rho=double(phi<0);
    if ~any(rho(:)),p=zeros(size(rho));break;end
    c=nutrient(rho,lambda,cB,dx,dy);
    p=pressure(rho,G0*c,dx,dy);
    [px,py]=gradient(p,dx,dy);[fx,fy]=gradient(phi,dx,dy);
    magnitude=sqrt(fx.^2+fy.^2+eps);
    normal_speed=-(px.*fx+py.*fy)./magnitude;
    [~,speed]=fast_marching(phi,normal_speed,dx,dy);
    [xp,xm,yp,ym]=differences(phi,dx,dy);
    plus=sqrt(max(max(xm,0).^2,min(xp,0).^2)+max(max(ym,0).^2,min(yp,0).^2));
    minus=sqrt(max(min(xm,0).^2,max(xp,0).^2)+max(min(ym,0).^2,max(yp,0).^2));
    phi=phi-dt*(max(speed,0).*plus+min(speed,0).*minus);
    info.max_courant=max(info.max_courant,dt*max(abs(speed(:)))*(1/dx+1/dy));
    if mod(step,20)==0,[phi,~]=fast_marching(phi,zeros(size(phi)),dx,dy);end
    assert(all(isfinite(phi(:))),'Nonfinite level-set solution.');
end
rho=double(phi<0);
end

function [xp,xm,yp,ym]=differences(phi,dx,dy)
xp=([phi(:,2:end),phi(:,end)]-phi)/dx;
xm=(phi-[phi(:,1),phi(:,1:end-1)])/dx;
yp=([phi(2:end,:);phi(end,:)]-phi)/dy;
ym=(phi-[phi(1,:);phi(1:end-1,:)])/dy;
end

function c=nutrient(rho,lambda,cB,dx,dy)
[ny,nx]=size(rho);N=numel(rho);indices=reshape(1:N,ny,nx);
center=indices(2:end-1,2:end-1);rows=center(:);
up=indices(1:end-2,2:end-1);down=indices(3:end,2:end-1);
left=indices(2:end-1,1:end-2);right=indices(2:end-1,3:end);
inside=rho(rows)>0;diagonal=ones(N,1);
diagonal(rows)=2/dx^2+2/dy^2+1;
diagonal(rows(inside))=2/dx^2+2/dy^2+lambda;
A=sparse([(1:N)';rows;rows;rows;rows],[(1:N)';up(:);down(:);left(:);right(:)], ...
    [diagonal;repmat(-1/dy^2,2*numel(rows),1);repmat(-1/dx^2,2*numel(rows),1)],N,N);
b=cB*ones(N,1);b(rows(inside))=0;c=reshape(A\b,ny,nx);
end

function p=pressure(rho,source,dx,dy)
[ny,nx]=size(rho);N=numel(rho);indices=reshape(1:N,ny,nx);
mask=rho>0;mask([1 end],:)=false;mask(:,[1 end])=false;
rows=indices(mask);diagonal=ones(N,1);diagonal(rows)=2/dx^2+2/dy^2;
A=sparse([(1:N)';rows;rows;rows;rows], ...
    [(1:N)';rows-1;rows+1;rows-ny;rows+ny], ...
    [diagonal;repmat(-1/dy^2,2*numel(rows),1);repmat(-1/dx^2,2*numel(rows),1)],N,N);
b=zeros(N,1);b(rows)=source(rows);p=reshape(A\b,ny,nx);
end
