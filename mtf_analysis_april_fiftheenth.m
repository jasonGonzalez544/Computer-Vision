close all;
clear;

%load images
folder = dir('Experimental data\\MTF_heights\\2025-04-15_MTFs');
nfiles = length(folder);
%remove non-tifs, first file from folder for ease of counting
for i=1:3
    folder(1) = [];
end

%get distance per pixel
pixel_space = 871; %pipxel count
measured_space = 2700;  %measured distance in um
dpp = floor(measured_space / pixel_space);

%simulate imge of particles with calculated size
circ_bg = uint32(ones(1000,1000)*72);
cx = 1:1000;
cy = 1:1000;
[imcols, imrows] = meshgrid(cx , cy);
centerX = [100 300 500 700 500 255 352];
centerY = [300 100 800 500 500 773 100];
radius = [33 16 33 16 33 16 33];
circlePixels = any((imrows(:) - centerY).^2 ...
    + (imcols(:) - centerX).^2 <= radius.^2, 2);
circlePixels = reshape(circlePixels, 1000, 1000);
circ_bg(circlePixels) = uint32(12);
% figure
% imagesc(circ_bg)
% title('Simulated particles')
% saveas(gca,'Experimental data\\MTF_heights\\figures\\simulated_circles.png')

% iterate over data
items = [1:3 ; 477:479 ; 830:832 ; 1281:1283; 1734:1736];
%items = [1:3];
figure;
hax1 = axes;
hold on
for i=1:length(items)
    %label for plotting
    lab = ((i-1)/2)*100;
%for i=1
%   load in raw images
    im_raw_1 = uint32(imread(folder(items(i,1)).name));
    im_raw_2 = uint32(imread(folder(items(i,2)).name));
    im_raw_3 = uint32(imread(folder(items(i,3)).name));
%   average for noise reduction
    im = (im_raw_1 + im_raw_2 + im_raw_3) / 3;
%   flip for analysis
    im = flip(im,2);
    % figure
    % imagesc(im)
    % cb = colorbar();
    % ylabel(cb, 'intensity','Rotation',270)
    % axis equal
    % xlabel('Pixels')
    % ylabel('Pixels')
    % set(gca,"FontName","times")
    % saveas(gca, sprintf('Experimental data\\MTF_heights\\figures\\edge_%d.png' , i))

%   1st derivative of line (LSF) for ESF
    ROI = double(im(300,350:450));
    % figure
    % plot(ROI / max(ROI))
    % xlabel('pixels')
    % ylabel('normalized intensity')
    % title(sprintf('ESF for %.1d um displacement' , lab))
    % xlim([1,length(ROI)])
    % set(gca,"FontName","times")
    % saveas(gca,sprintf('Experimental data\\MTF_heights\\figures\\ESF_%d.png' , i))

%   Line sprad function
    LSF = diff(double(ROI));
    LSF(LSF<0) = 0;
    % figure
    % plot(ESF)
    % xlabel('pixels')
    % ylabel('signal difference')
    % title(sprintf('LSF for %.1d um displacement' , lab))
    % saveas(gca , sprintf('Experimental data\\MTF_heights\\figures\\LSF_%d.png' , i))
    %close
    
%   fit to gaussian (PSF)
    f = fit((1:length(LSF)).' , LSF.' , 'gauss1');
    % % figure
    % plot(ESF)
    % hold on
    % plot(f)
    % % close
    % hold off
    % xlabel('pixels')
    % ylabel('signal difference')
    % title(sprintf('LSF fit for %.1d um displacement' , lab))
    % saveas(gca , sprintf('Experimental data\\MTF_heights\\figures\\fit_%d.png' , i))

    %use gaussian parameters to find convolution kernel (PSF)
    std = f.c1;
    amp = f.a1;
    x = -ceil(3*std):ceil(3*std); 
    G1 = exp(-(x/std).^2);  % Create 1D gaussian 
    %plot(G1)
    x_plot = -25:25;
    G1_plot = exp(-(x_plot/std).^2);
    %plot(G1_plot , 'DisplayName',sprintf('height = +%i' , i-1))
    G2 = G1'*G1;  % Compute the 2D gaussian out of 1D gaussian.
    % figure
    % imagesc(G2)
    % axis equal
    % xlabel('pixels')
    % ylabel('pixels')
    % title(sprintf('PSF for %.1d um displacement\nStandard Deviation: %.2f pixels' , lab,std))
    % saveas(gca , sprintf('Experimental data\\MTF_heights\\figures\\PSF_%d.png' , i))

    %convolve simulated particles with PSF
    sim = conv2(circ_bg , G2, "same");
    sim = (sim/max(sim,[],"all"));
    noise_sim = imnoise(sim, "gaussian",0,0.001);
    % figure
    % imagesc(noise_sim)
    % xlabel('pixels')
    % ylabel('pixels')
    % title(sprintf('Simulated particles for %.1d um displacement' , lab))
    % saveas(gca , sprintf('Experimental data\\MTF_heights\\figures\\noisy_circs_%d.png' , i))

    %deconvolve
    signal = mean(im(:));
    noise = 0.0001;
    % nsr = 10000000 / std2(double(im(:)));
    nsr = 3000;

    im_deconv = deconvlucy(noise_sim,G2);
    % figure
    % imagesc(im_deconv)
    % xlabel('pixels')
    % ylabel('pixels')
    % title(sprintf('deconvolved particles for %.1d um displacement' , lab))
    % saveas(gca , sprintf('Experimental data\\MTF_heights\\figures\\recreated_circs_lucy_%d.png' , i))
    
    %transform PSF to get MTF
    MTF = fft(ifftshift(G1_plot));
    MTF = MTF/max(MTF);
    Fs = 0:(1/(length(G1_plot)-1)):1;
    Fs = Fs*(1/dpp);

    plot(hax1,Fs , abs(MTF),'DisplayName',sprintf('%.1d um' , lab),'LineWidth',1)
    xlim(hax1,[0, 0.5*(1/dpp)])
    xlabel(hax1,'Cycles per um')
    ylabel(hax1,'Modulation Transfer')
    title(hax1,'MTF from focal plane')
    legend(hax1)


end
legend
set(hax1,"FontName","times")
% saveas(hax1 ,'Experimental data\\MTF_heights\\figures\\MTFs.png')