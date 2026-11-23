function [distance, extended] = fast_marching(phi, values, dx, dy)
%FAST_MARCHING First-order Eikonal distance and normal velocity extension.
% Subcell interface crossings initialize a min-heap ordered upwind march.
% Each sign of phi is marched separately, using the same interface seeds.
[ny,nx]=size(phi); N=numel(phi); inside=phi<0;
sx=inf(size(phi));sy=sx;
cross=inside(:,1:end-1)~=inside(:,2:end);
left=abs(phi(:,1:end-1));right=abs(phi(:,2:end));den=left+right;
a=dx*left./max(den,eps);b=dx*right./max(den,eps);a(~cross)=inf;b(~cross)=inf;
sx(:,1:end-1)=min(sx(:,1:end-1),a);sx(:,2:end)=min(sx(:,2:end),b);
cross=inside(1:end-1,:)~=inside(2:end,:);
low=abs(phi(1:end-1,:));high=abs(phi(2:end,:));den=low+high;
a=dy*low./max(den,eps);b=dy*high./max(den,eps);a(~cross)=inf;b(~cross)=inf;
sy(1:end-1,:)=min(sy(1:end-1,:),a);sy(2:end,:)=min(sy(2:end,:),b);
d=1./sqrt(1./sx.^2+1./sy.^2);d(phi==0)=0;
accepted=isfinite(d); seeds=find(accepted); extended=values;
if isempty(seeds),distance=phi;return;end
heap=zeros(N,1);position=zeros(N,1);heap_size=0;
for q=1:numel(seeds), update_neighbors(seeds(q));end
while heap_size>0
    index=heap(1);last=heap(heap_size);heap_size=heap_size-1;position(index)=0;
    if heap_size>0
        heap(1)=last;position(last)=1;k=1;
        while 2*k<=heap_size
            child=2*k;
            if child<heap_size && d(heap(child+1))<d(heap(child)),child=child+1;end
            if d(heap(k))<=d(heap(child)),break;end
            swap(k,child);k=child;
        end
    end
    accepted(index)=true;update_neighbors(index);
end
distance=d;distance(inside)=-d(inside);

    function swap(a,b)
        tmp=heap(a);heap(a)=heap(b);heap(b)=tmp;
        position(heap(a))=a;position(heap(b))=b;
    end
    function update_neighbors(index)
        row=mod(index-1,ny)+1;col=floor((index-1)/ny)+1;
        if row>1,relax(index-1);end
        if row<ny,relax(index+1);end
        if col>1,relax(index-ny);end
        if col<nx,relax(index+ny);end
    end
    function relax(index)
        if accepted(index),return;end
        row=mod(index-1,ny)+1;col=floor((index-1)/ny)+1;
        ax=inf;ay=inf;vx=0;vy=0;phase=inside(index);
        if col>1
            j=index-ny;if accepted(j)&&inside(j)==phase,ax=d(j);vx=extended(j);end
        end
        if col<nx
            j=index+ny;if accepted(j)&&inside(j)==phase&&d(j)<ax,ax=d(j);vx=extended(j);end
        end
        if row>1
            j=index-1;if accepted(j)&&inside(j)==phase,ay=d(j);vy=extended(j);end
        end
        if row<ny
            j=index+1;if accepted(j)&&inside(j)==phase&&d(j)<ay,ay=d(j);vy=extended(j);end
        end
        if ~isfinite(ax)&&~isfinite(ay),return;end
        if ~isfinite(ay),next=ax+dx;value=vx;
        elseif ~isfinite(ax),next=ay+dy;value=vy;
        else
            next=min(ax+dx,ay+dy);
            if next>max(ax,ay)
                A=1/dx^2+1/dy^2;B=ax/dx^2+ay/dy^2;
                C=ax^2/dx^2+ay^2/dy^2-1;
                next=(B+sqrt(max(0,B^2-A*C)))/A;
            end
            wx=max(next-ax,0)/dx^2;wy=max(next-ay,0)/dy^2;
            value=(wx*vx+wy*vy)/max(wx+wy,eps);
        end
        if next>=d(index),return;end
        d(index)=next;extended(index)=value;
        k=position(index);
        if k==0,heap_size=heap_size+1;k=heap_size;heap(k)=index;position(index)=k;end
        while k>1
            parent=floor(k/2);if d(heap(parent))<=d(heap(k)),break;end
            swap(k,parent);k=parent;
        end
    end
end
