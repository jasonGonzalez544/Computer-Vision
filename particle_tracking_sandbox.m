close all;
clear;

%load images
folder = dir('narrowday3exp1_out');
nfiles = length(folder);

%load background
bg0 = imread(folder(3).name);

%loop through
% for ii=8000:4000:50000
for ii = 40000
%for ii=53000
    %read images
    im0 = imread(folder(ii).name);
    % figure
    % imshow(im0)
    % saveas(gca , sprintf('Experimental data\\Particle tracking\\%d_original.tiff' , ii))

    %show original image
    %figure;
    %imagesc(im0);
    
    %substract image from background
    im1 = bg0 - im0;
    
    %binarize image
    threshold = 40;
    im1(im1 <= threshold) = 0;
    im1(im1 > threshold) = 1;
    
    %erode image
    radius = 4;
    rel = strel("disk" , radius);
    im2 = imerode(im1 , rel);

    
%       %Method 1: only erosion
%     count objects in image
%     cc = bwconncomp(im2 , 8);
%     
%     overlay counted particles to image
%     figure;
%     A = labeloverlay(im0 , im2*500 , 'Colormap',"hot");
%     imshow(A)
%     title(sprintf('%d objects', cc.NumObjects));
%     saveas(gca , sprintf('Experimental data\\Particle tracking\\%d.tiff' , ii))
    
    %   Method 2: Distance transform of eroded image
    dist2 = bwdist(~im2);
    L = watershed(-dist2);
    L(~im2) = 0;
    figure;
    imagesc(dist2)
    axis equal
    saveas(gca , sprintf('Experimental data\\Particle tracking\\Secondpass\\thresh_%d.tiff' , ii))
    %A2 = labeloverlay(im0 , L);
    %imshow(A2)
    
    %   Method 3: foreground vs background distance transform
    diff = im1-im2;
    dist3 = bwdist(diff);
    dist3(~im2) = 0;
    dist3 = dist3;
    L3 = watershed(-dist3);
    L3(~im2) = 0;
%     figure;
%     imagesc(dist3);
%     saveas(gca , sprintf('Experimental data\\Particle tracking\\Secondpass\\diff_dist_%d.tiff' , ii))
%     figure;
%     imagesc(L3);
%     figure;
%     imagesc(diff)
%     saveas(gca , sprintf('Experimental data\\Particle tracking\\Secondpass\\thresh_diff_%d.tiff' , ii))
    
    figure;
    A3 = labeloverlay(im0 , L3);
    imshow(A3)
%     saveas(gca , sprintf('Experimental data\\Particle tracking\\Secondpass\\secondpass_%d.tiff' , ii))
    %figure;
    %imshow(im0)
%     
end