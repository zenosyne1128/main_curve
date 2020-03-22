function pc_projection=pc_project_multidim_3d(X,Xinit,kernel_sigma,targetdim,maxiter)

point=Xinit;
N=size(X,2);
dim=size(X,1);
threshold=1e-3;
rangex = [min(X(1,:)),max(X(1,:))];
rangey = [min(X(2,:)),max(X(2,:))];
rangez = [min(X(3,:)),max(X(3,:))];
for ind=1:size(point,2)
    flag=0;
    [p,gra,H,SI] = pgh_3d(point(:,ind),X',N,dim,kernel_sigma);
    point_pc1=point(:,ind);
    [V,D]= eig(SI);
    [~,indexsorted]=sort(diag(abs(D)),'descend');
    %降序排列特征值并标号
    ConstrainedSpace = V(:,indexsorted(1:end-targetdim));
    for tim=1:size(ConstrainedSpace,2)
        direction = ConstrainedSpace(:,tim);
        direction = direction * sign(direction'*gra);
        ConstrainedSpace(:,tim)=direction;
        %调整特征向量符号
    end
    if abs(gra'*H*gra/(norm(gra'*H)*norm(gra'*H)))<0.01
        %此时g在H下的投影与g正交，即g为H某一特征向量，xind在主曲面上
        flag=1;
        pc(ind,:)= point_pc1;
    end
        if ~flag
            for a=1:maxiter  %单点最大迭代次数
                G = kernel_matrix(point_pc1,X',dim,1,N,kernel_sigma);
                num1 = sum(repmat(G,dim,1).*X,2);
                den1 = sum(G,2);
                if(direction'*gra<0),keyboard;end %异常输出
               %[p1,gra1,H,SI] = pgh_3d(point_pc1,X',N,dim,kernel_sigma);
                point_pc1_old=point_pc1;
                for c=1:size(ConstrainedSpace,2)
                    direction=ConstrainedSpace(:,c);
                    point_pc1 = point_pc1 + direction * (direction'*(num1/den1-point_pc1));
                end

                if abs(point_pc1_old-point_pc1)<threshold,break;end
            end
            pc(ind,:)= point_pc1;
        end
       
        pc_projection=pc;
      

end

%Xdeflated=X-repmat(direction,1,N).*repmat((direction'*X),dim,1);

function K=kernel_matrix(data1,data2,n,N1,N2,kernel_sigma)
for j=1:N1,
    C=1./((2*pi)^(n/2)*kernel_sigma.^n);
    dx=(repmat(data1(:,j),1,N2)-data2')./repmat(kernel_sigma,n,N2);
    K(j,:)=C.*exp(-0.5*sum(dx.^2,1));
end