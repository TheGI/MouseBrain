%% Get image files
clear all;
close all;

imgDir = uigetdir;
imgDir = [imgDir '/'];
imgList=[dir([imgDir '*.jpg']);dir([imgDir '*.png']);dir([imgDir '*.tif'])];
n_img=length(imgList);

%%
img = imread([imgDir imgList(23).name]);
figure(1), imshow(img), title('img');

%%
load regioncoordinates;

nColors = 5; %'background','red','green','blue','yellow'
sample_regions = false([size(img,1) size(img,2) nColors]);

for count = 1:nColors
    sample_regions(:,:,count) = roipoly(img);
end


%%
lab_fabric = rgb2lab(img);
a = lab_fabric(:,:,2);
b = lab_fabric(:,:,3);
color_markers = zeros([nColors, 2]);

for count = 1:nColors
    color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
    color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
end

% 1 = background, 2 = red, 3 = green, 4 = blue, 5 = yellow
color_labels = 1:nColors;
a = double(a);
b = double(b);
distance = zeros([size(a), nColors]);

for count = 1:nColors
    distance(:,:,count) = ( (a - color_markers(count,1)).^2 + ...
        (b - color_markers(count,2)).^2 ).^0.5;
end

[~, label] = min(distance,[],3);
label = color_labels(label);
clear distance;

%%
rgb_label = repmat(label,[1 1 3]);
segmented_images = zeros([size(img), nColors],'uint8');

for count = 1:nColors
    color = img;
    color(rgb_label ~= color_labels(count)) = 0;
    segmented_images(:,:,:,count) = color;
end
subplot(2,3,1);
imshow(img),title('total image');
subplot(2,3,2);
imshow(segmented_images(:,:,:,1)), title('background');
subplot(2,3,3);
imshow(segmented_images(:,:,:,2)), title('red');
subplot(2,3,4);
imshow(segmented_images(:,:,:,3)), title('green');
subplot(2,3,5);
imshow(segmented_images(:,:,:,4)), title('blue');
subplot(2,3,6);
imshow(segmented_images(:,:,:,5)), title('yellow');

%%
plot_labels = {'k', 'r', 'g', 'b', 'y'};

figure;
for count = 1:nColors
    plot(a(label==count),b(label==count),'.','MarkerEdgeColor', ...
        plot_labels{count}, 'MarkerFaceColor', plot_labels{count});
    hold on;
end

title('Scatterplot of the segmented pixels in ''a*b*'' space');
xlabel('''a*'' values');
ylabel('''b*'' values');

%%
totalcount = numel(find(label == 2 | label == 3 | label == 4 | label == 5));
colorprop(1) = numel(find(label == 2)) / totalcount;
colorprop(2) = numel(find(label == 3)) / totalcount;
colorprop(3) = numel(find(label == 4)) / totalcount;
colorprop(4) = numel(find(label == 5)) / totalcount;

bar(colorprop);

%%
segmented_images = zeros([n_img, size(img), nColors],'uint8');
    for count = 1:nColors
        color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
        color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
    end
    
for i=1:n_img
    img=imread([imgDir imgList(i).name]);
    
    lab_fabric = rgb2lab(img);
    a = lab_fabric(:,:,2);
    b = lab_fabric(:,:,3);
    color_markers = zeros([nColors, 2]);
    

    
    % 1 = background, 2 = red, 3 = green, 4 = blue, 5 = yellow
    color_labels = 1:nColors;
    a = double(a);
    b = double(b);
    distance = zeros([size(a), nColors]);
    
    for count = 1:nColors
        distance(:,:,count) = ( (a - color_markers(count,1)).^2 + ...
            (b - color_markers(count,2)).^2 ).^0.5;
    end
    
    [~, label] = min(distance,[],3);
    label = color_labels(label);
    clear distance;
    
    totalcount = numel(find(label == 2 | label == 3 | label == 4 | label == 5));
    for j = 1:4
        colorprop(i,j) = numel(find(label == j+1)) / totalcount;
    end
    
    rgb_label = repmat(label,[1 1 3]);
    
    
    for count = 1:nColors
        color = img;
        color(rgb_label ~= color_labels(count)) = 0;
        segmented_images(i,:,:,:,count) = color;
    end
    
    disp(['Progress: ' num2str(i) '/' num2str(n_img)]);
end

%%
i = 1;
img=imread([imgDir imgList(i).name]);
subplot(2,3,1);
imshow(img),title('total image');
subplot(2,3,2);
imshow(reshape(segmented_images(i,:,:,:,1),1440,1920,3)), title('background');
subplot(2,3,3);
imshow(reshape(segmented_images(i,:,:,:,2),1440,1920,3)), title('red');
subplot(2,3,4);
imshow(reshape(segmented_images(i,:,:,:,3),1440,1920,3)), title('green');
subplot(2,3,5);
imshow(reshape(segmented_images(i,:,:,:,4),1440,1920,3)), title('blue');
subplot(2,3,6);
imshow(reshape(segmented_images(i,:,:,:,5),1440,1920,3)), title('yellow');

%%
C = mean(colorprop(1:17,:),1);
F = mean(colorprop(18:34,:),1);

R(:,1) = colorprop(1:17,1);
R(:,2) = colorprop(18:34,1);
G(:,1) = colorprop(1:17,2);
G(:,2) = colorprop(18:34,2);
B(:,1) = colorprop(1:17,3);
B(:,2) = colorprop(18:34,3);
Y(:,1) = colorprop(1:17,4);
Y(:,2) = colorprop(18:34,4);

[hR, pR, cR, sR] = ttest(R(:,1),R(:,2));
[hG, pG, cG, sG] = ttest(G(:,1),G(:,2));
[hB, pB, cB, sB] = ttest(B(:,1),B(:,2));
[hY, pY, cY, sY] = ttest(Y(:,1),Y(:,2)); 
