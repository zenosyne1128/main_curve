function [p, G, H, SI] = pgh_3d(x,data,N,dim,kernel_sigma)
[dimx,Nx]=size(x);
for i=1:Nx
    Kvector = kernel_matrix(x(:,i),data',dim,1,N,kernel_sigma);
    Xvector = (repmat(x(:,i),1,N)-data');
    GX = (Xvector.*repmat(Kvector,dim,1));
    p= sum(Kvector,2)/N;
    G = -1/(N*(kernel_sigma^2))*sum(GX,2);
    H = 1/(N*(kernel_sigma^4)) * (GX*Xvector') - eye(dim)*p/(kernel_sigma^2);
    SI= (1/p)*H-(1/p)^2*G*G';
end
function K=kernel_matrix(data1,data2,n,N1,N2,kernel_sigma)
for j=1:N1
    %C=1./((2*pi)^(n/2)*kernel_sigma.^n);
    dx=(repmat(data1(:,j),1,N2)-data2)./repmat(kernel_sigma,n,N2);
    K(j,:)=exp(-0.5*sum(dx.^2,1));
end