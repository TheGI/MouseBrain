%% Get image files
clear all;
close all;

imgDir = uigetdir;
imgDir = [imgDir '/'];
imgList=[dir([imgDir '*.jpg']);dir([imgDir '*.png']);dir([imgDir '*.tif'])];
n_img=length(imgList);
n_group = 3;
n_indiv = 9;

rgb = [1,0,0;0,1,0;0,0,1;1,1,0;1,1,1];
lab = rgb2lab(rgb);

%%
he = imread([imgDir imgList(1).name]);
imshow(he), title('H&E image');

%%
cform = makecform('srgb2lab');
lab_he = applycform(he,cform);

ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

%%
nColors = 5;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3,'Start',repmat(lab(:,2:3),1,1,3));
%%                                  
pixel_labels = reshape(cluster_idx,nrows,ncols);
imshow(pixel_labels,[]), title('image labeled by cluster index');

%%
segmented_images = cell(1,5);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = he;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end
subplot(1,5,1);
imshow(segmented_images{1}), title('objects in cluster 1');
subplot(1,5,2);
imshow(segmented_images{2}), title('objects in cluster 2');
subplot(1,5,3);
imshow(segmented_images{3}), title('objects in cluster 3');
subplot(1,5,4);
imshow(segmented_images{4}), title('objects in cluster 4');
subplot(1,5,5);
imshow(segmented_images{5}), title('objects in cluster 5');

%%
mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
blue_cluster_num = idx(1);

L = lab_he(:,:,1);
blue_idx = find(pixel_labels == 5);
L_blue = L(blue_idx);
is_light_blue = imbinarize(L_blue);

nuclei_labels = repmat(uint8(0),[nrows ncols]);
nuclei_labels(blue_idx(is_light_blue==false)) = 1;
nuclei_labels = repmat(nuclei_labels,[1 1 3]);
blue_nuclei = he;
blue_nuclei(nuclei_labels ~= 1) = 0;
imshow(blue_nuclei), title('blue nuclei');