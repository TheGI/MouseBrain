%% Get image files
clear all;
close all;

imgDir = uigetdir;
imgDir = [imgDir '/'];
imgList=[dir([imgDir '*.jpg']);dir([imgDir '*.png']);dir([imgDir '*.tif'])];
n_img=length(imgList);
n_group = 2;
n_indiv = 17;

%% Analyze Image Colors
% created sorted color map
% number of bins for each R,G,B
n_bin=10;
% total number of colors
n_color = n_bin^3;
% define edges for color histogram
edg=(0:n_bin)/n_bin;
% define middle point of edges
edgM = edg(2:end) - diff(edg)/2;

% create rgb bins
i_rgb = 1;
for b = 1:n_bin
    for g = 1:n_bin
        for r = 1:n_bin
            rgb(i_rgb,:) = [edgM(r), edgM(g), edgM(b)];
            i_rgb = i_rgb + 1;
        end
    end
end

% sort rgb bins into more human perceptible rainbow scale
hsv = rgb2hsv(rgb);
hsv = sortrows(hsv,[-1 -3 2]);
rgb = hsv2rgb(hsv);
cmap = colormap(rgb);

% this is to inclue the boundary value
edg(end) = edg(end) + 1e-5;

% create histogram matrices
hTotal=zeros(n_img,n_color);
hTotalraw=zeros(n_img,n_color);

% go through each image and bins colors & create histograms
for i=1:n_img
    img=imread([imgDir imgList(i).name]);  
    
    % mask to remove background color
    maskR = (img(:,:,1) > 245);
    maskG = (img(:,:,2) > 245);
    maskB = (img(:,:,3) > 245);
    maskRGB = maskR & maskG & maskB;
     
    img=im2double(img);
    temp = rgb2ind(img,cmap);
    temp = temp(~maskRGB);
    [U,~,ic] = unique(temp);
    UC = accumarray(ic(:),1);
    hTotal(i,U) = UC / sum(UC);
    hTotalraw(i,U) = UC;
    
    disp(['Progress: ' num2str(i) '/' num2str(n_img)]);
end
%% Overview of Group Color Histogram
i_hMean = 1;
hMean = [];
for i = 1:n_indiv:n_group*n_indiv
    hMean(i_hMean,:) = mean(hTotal(i:i+n_indiv-1,:));
    i_hMean = i_hMean+1;
end

i_c = find(sum(hMean,1) ~= 0);
i_c1 = find(hMean(1,:) ~= 0);
i_c2 = find(hMean(2,:) ~= 0);
%i_c3 = find(hMean(3,:) ~= 0);
%i_c4 = find(hMean(4,:) ~= 0);

hold on;
for i = 1:length(i_c)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,1,1],cmap(i_c(i),:),'FaceAlpha', 1);
end

p1 = plot(hMean(1,i_c),'w^-','MarkerSize',8);
p2 = plot(hMean(2,i_c),'wo-','MarkerSize',8);
%p3 = plot(hMean(3,i_c),'ws-','MarkerSize',8);
%p4 = plot(hMean(4,i_c),'w*-','MarkerSize',8);
hold off;
lgd = legend([p1,p2],...
    'C','F',...
    'Location','Best');
% lgd = legend([p1,p2,p3,p4],...
%     'T0','T1','T2','T3',...
%     'Location','Best');
lgd.Color = [0 0 0];
lgd.TextColor = [1 1 1];
xlabel('Colors'), ylabel('Count');
title('Histogram of Mouse Brain');
%% Which color passes ANOVA?
hTotal3D = reshape(hTotal,n_indiv,n_group,n_color);

i_cSig = [];
for j = 1:n_color
    p = anova1(hTotal3D(:,:,j),[],'off');
    if p < 0.01 && max(mean(hTotal3D(:,:,j))) > 0.0005
        i_cSig(end+1) = j;
    end
end

hold on;
for i = 1:length(i_cSig)
    patch([i-.5,i+.5,i+.5,i-.5],[0,0,1,1],cmap(i_cSig(i),:),'FaceAlpha', 1);
end

p1Sig = plot(hMean(1,i_cSig),'w^-','MarkerSize',10);
p2Sig = plot(hMean(2,i_cSig),'wo-','MarkerSize',10);
% p3Sig = plot(hMean(3,i_cSig),'ws-','MarkerSize',10);
% p4Sig = plot(hMean(4,i_cSig),'w*-','MarkerSize',10);
hold off;
% lgd = legend([p1Sig,p2Sig,p3Sig,p4Sig],...
%     'T0','T1','T2','T3',...
%     'Location','Best');
lgd = legend([p1Sig,p2Sig],...
    'C','F',...
    'Location','Best');
lgd.Color = [0 0 0];
lgd.TextColor = [1 1 1];
xlabel('Colors'), ylabel('Count');
title('Histogram of Significantly Different Colors (p < 0.01)');

%% Chi-Square Histogram Distance
dist_func=@chi_square_statistics_fast;
D=pdist2(hTotal,hTotal,dist_func);
hm = HeatMap(D);

%% Eigenspace Analysis
cSig = cmap(i_cSig,:);
c1 = cmap(i_c1,:);
c2 = cmap(i_c2,:);
% c3 = cmap(i_c3,:);
% c4 = cmap(i_c4,:);
cSigM = mean(cSig);
cSigV = cov(cSig);
[V,D] = eig(cSigV);
S = zeros(3,3);
S(:,1) = V(:,3)*D(3,3);
S(:,2) = V(:,2)*D(2,2);
S(:,3) = V(:,1)*D(1,1);

S = S/norm(S);

scatter3(cSig(:,1),cSig(:,2),cSig(:,3),200,cSig,'filled');
xlabel('red'), ylabel('green'), zlabel('blue');
title('Most Significantly Different Colors with Eigenvectors');
axis equal;
view(19,31);
hold on;
quiver3(repmat(cSigM(1),1,3),repmat(cSigM(2),1,3),repmat(cSigM(3),1,3),...
    S(1,:),S(2,:),S(3,:));
hold off;
grid on;

% plot mean color point clouds for T0, T1, T2, T3
figure;
hold on;
scatter3(c1(:,1) + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    c1(:,2) + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    c1(:,3) + + (rand(length(c1(:,1)),1)-0.5)*0.05,...
    100,'r','filled','o');
scatter3(c2(:,1) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    c2(:,2) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    c2(:,3) + (rand(length(c2(:,1)),1)-0.5)*0.05,...
    100,'g','filled','s');
% scatter3(c3(:,1) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
%     c3(:,2) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
%     c3(:,3) + (rand(length(c3(:,1)),1)-0.5)*0.05,...
%     100,'b','filled','^');
% scatter3(c4(:,1) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
%     c4(:,2) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
%     c4(:,3) + (rand(length(c4(:,1)),1)-0.5)*0.05,...
%     100,c4,'filled','^');
% legend('T0','T1','T2','T3','Location','Best');
legend('C','F','Location','Best');
quiver3(repmat(cSigM(1),1,3),repmat(cSigM(2),1,3),repmat(cSigM(3),1,3),...
    S(1,:),S(2,:),S(3,:));
hold off;
xlabel('red'), ylabel('green'), zlabel('blue');
title('Mean Per Group Colors with Eigenvectors');
axis equal;
view(17,18);
grid on;

%% Principal Component Color Complements
for i = 1:n_group
    ind = dsearchn(cSig,S(:,i)'*10);
    eigcolors1(i,:)=cSig(ind,:);
    ind = dsearchn(cSig,-S(:,i)'*10);
    eigcolors2(i,:)=cSig(ind,:);
end

for i = 1:n_group
%     patch([i-1 i i i-1], [0 0 0.5 0.5], (repmat(0.5,3,1)+S(:,i)/norm(S(:,i))/2)')
%     patch([i-1 i i i-1], [0.5 0.5 1 1], (repmat(0.5,3,1)-S(:,i)/norm(S(:,i))/2)')
    patch([i-1 i i i-1], [0 0 0.5 0.5], eigcolors1(i,:))
    patch([i-1 i i i-1], [0.5 0.5 1 1], eigcolors2(i,:))
end