%% demo for 3d 
t = linspace(0,10*pi,500);
x = [5*sin(t),5*sin(t)+randn(1,500)];
y = [5*cos(t),5*cos(t)+0.8*randn(1,500)];
z = [t,t+0.3*randn(1,500)];
dataset=[x;y;z];
fig1 = plot3(x,y,z,'.')
print
% create a dataset of a helix line added noise in 3 dimension 
pc_projection=pc_project_multidim_3d(dataset,dataset,kernel_sigma,targetdim,20);
pc_projection = pc_projection';
figure(2),plot3(dataset(1,:),dataset(2,:),dataset(3,:),'.b');hold on,
title('PC projection with KDE pdf')
plot3(pc_projection(1,:),pc_projection(2,:),pc_projection(3,:),'.r');axis equal;
print('-f2','main curve of helix curve','-dpng')