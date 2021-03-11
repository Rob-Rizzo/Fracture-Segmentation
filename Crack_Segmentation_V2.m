% CRACK SEGMENTATION Version 2.0

% The code is built to automatically extract linear features (i.e., fractures)
% from two-dimensional images.

% New version of Crack Segmentation more specific for large scale images (photos)
% Main differences in this version are: 
% 1 - Applies a Gaussian Filter on the Gray-scale image to reduce noise
% 2 - Applies a Frangi vesselness filter during the pre-processing of the
% image to find the lineaments (i.e., fractures) in the analysed image.

% PARAMETERS TO CHANGE: 
% Gaussian Filter; 
% Median Filter Aperture; 
% Frangi Filter size; 
% Binarised bw_threshold; 
% Morphological Parameter: gap_length1 and conversion

% Author: Roberto Rizzo @ Uni of Edinburgh / Uni of Aberdeen
% Date: March 2021
% Based on Version 1.0 of the code, Luke Griffiths (crack-extract Python Code) and Franz Abes code (FRAg)


clear;
close all;

workspace;
format short g;
format compact;

%% ============= SET PARAMETERS =================
% The choice for the best parameters its a matter of trial and error. 

% FILTER and BINARIZATION PARAMETERS

% Gaussian Filter
gaus_filt = 2;

% Median Filter
aperture1 = [64 64];
aperture2 = [3 3];
% Frangi filter
frangi_size1 = 7;
frangi_size2 = 2;

%binarization threshold (1-100)
bw_threshold = 30;

% neighbourhoodsize (x,y) for wiener filter
wiener_area1 = [10 10];    
wiener_area2 = [2 2]; 

% min-max range to rescale grayscale values
rescale_range  = [0 40];    

% MORPHOLOGICAL PARAMETERS
% minimum fracture length in pixels (given by the major axes length, regionprops function)
frac_length = 20;    % default 8
% largest object size (nr of pixels) to be removed
object_limit  = round((2/3).*frac_length);  % default 2/3 of fracture length (if known)
% gap length closed by morphological dilation-erosion
gap_length1 = round((1/3).*frac_length);   % default 1/3 of fracture length (if known)
% Strel parameters used to perform the morphological close 
str_size = 5;
str_bound = 4;

% TARGETED CLOSING
% eccentricity threshold (0 is circular, 1 is line)
ecc_threshold = 0.95;
% variation in alignment between major axis of component and end-point of a second component (-angle_var to + angle_var, in degrees)
angle_var = 15;
% gap length between line segments closed by targeted connections
gap_length2 = 10; %default 10  

% PRUNING
% length of branches to be pruned
pruning_length  = frac_length; %default none
% last size-based filter
frac_length2 = 1.*frac_length;     %default 1


%% ============================ READ DATA - DO NOT MODIFY================================
[FILENAME, PATHNAME] = uigetfile({'*.*'},'Select an image');
im_original = imread([PATHNAME FILENAME]);

figure; 
imshow(im_original,[]);
axis on;
axis image;
axis tight;
box on;
set(gca, 'FontSize',14);
caption = sprintf('Input Image');
title(caption, 'Interpreter', 'None');
xlabel('X [px]');
ylabel('Y [px]');
%save input image
print('-djpeg', '-r300', 'Alexis_Input.jpeg');

% Get the dimensions of the image.
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[rows, columns, numberOfColorChannels] = size(im_original);
if numberOfColorChannels > 1
  % It's not really gray scale like we expected - it's color.
  % Use weighted sum of ALL channels to create a gray scale image.
  im_gr = rgb2gray(im_original);
else
   im_gr = im_original;
end


%% ========================== FILTERING ===============================
%FILTER1 - Filter out small components
tic;
% Gaussian Filter
im_gaus = imgaussfilt(im_gr, gaus_filt);
%median filter
im_filt1 = medfilt2(im_gaus, aperture1,'symmetric');
% difference image
im_filt1 = im_filt1 - im_gaus;
% Frangi filter
im_frangi1 = fibermetric(im_filt1, frangi_size1, 'ObjectPolarity','bright');
% remove noise
im_wien1 = wiener2(im_frangi1, wiener_area1);      

%FILTER2 - Filter out large components
% Gaussian Filter
im_gaus2 = imgaussfilt(im_gr, gaus_filt);
%median filter
im_filt2 = medfilt2(im_gaus, aperture2,'symmetric');
% difference image
im_filt2 = im_filt2 - im_gaus2;
% Frangi filter
im_frangi2 = fibermetric(im_filt2, frangi_size2, 'ObjectPolarity','bright');
% remove noise
im_wien2 = wiener2(im_frangi2, wiener_area2);

% Final filtred image

im_final = im_wien1 + im_wien2;

toc;

figure; 
subplot(1,2,1);
imshow(im_gr, []);
set(gca, 'FontSize',14);
title('Original Gray-scale Image');
axis on;
axis image;
axis tight;
box on;
xlabel('X [px]');
ylabel('Y [px]');

subplot(1,2,2);
imshow(im_final, []);
set(gca, 'FontSize',14);
title('Binarized Image');
axis on;
axis image;
axis tight;
box on;
xlabel('X [px]');
ylabel('Y [px]');

%save the filtered images
%print('-djpeg', '-r300', 'FileName_Filtered.jpeg');
%% ========================== BINARIZATION ===============================

im_bin = imbinarize(im_final);

im_bin = bwareaopen(im_bin, object_limit);

tic;
% angle of the line object
% for alpha = 1:181  
%     % set orientation of line object
%     se = strel('line', gap_length1, alpha);     
%     % close gaps between line segments by dilation and erosion
%     im_close = imclose(im_bin, se);                                          
% end

se = strel('disk',str_size,str_bound);
im_close = imclose(im_bin, se);

 % skeletonize binary image
im_thin = bwmorph(im_close,'thin', Inf);                                 
toc;

figure; 
subplot(1,2,1);
imshow(im_close, []);
set(gca, 'FontSize',14);
title('Small Object Removed');
axis on;
axis image;
axis tight;
box on;
xlabel('X [px]');
ylabel('Y [px]');

subplot(1,2,2);
imshow(im_thin, []);
set(gca, 'FontSize',14);
title('Skeletonised Image');
axis on;
axis image;
axis tight;
box on;
xlabel('X [px]');
ylabel('Y [px]');

%save the skelotnised and binarised images
%print('-djpeg', '-r300', 'FileName_Skeletonise.jpeg');

%% =================== CONNECTION OF LINEAR OBJECTS ======================
% ------------------------ MORPHOLOGICAL CLOSING ------------------------
% ------------------------- TARGETED CLOSING -----------------------------

% GET CONNECTED COMPONENT PROPERTIES
% get properties of connected components
STATS = regionprops(im_thin,'Eccentricity','MajorAxisLength','Orientation','Centroid');                         
STATS = struct2table(STATS);
% get indices for all connected components
CC = bwconncomp(im_thin); 
% get end-points for ALL connected components
ENDS = bwmorph(im_thin, 'endpoints');    
% x,y-coordinates of end-points
[y_ends,x_ends] = find(ENDS == 1);                                          

% filter connected components
I = find(STATS.MajorAxisLength >= frac_length & STATS.Eccentricity > ecc_threshold);                    

% set image that will be connected
im_bin_conn = im_thin;                                                  

disp(['Found ' num2str(length(I)) ' suited connected objects']);

%--- LOOP OVER FILTERED CONNECTED COMPONENTS-------
for j = 1:length(I)
    
   % disp(['analysing object ' num2str(j) ' of ' num2str(length(I))]);
    
    im_work = im_bin;
    % remove component under consideration from work image
    im_work(CC.PixelIdxList{I(j)}) = 0;  
    % get end-notes of component under consideration
    ENDS_obj  = bwmorph(im_bin - im_work, 'endpoints');
    % get x-y coordinates of end-nodes
    [y_obj,x_obj] = find(ENDS_obj == 1);  
    
%--- END-POINTS OF CONSIDERED COMPONENT
    for i = 1:length(y_obj)
        % set end-points of considered component to zero
        J = find(y_ends == y_obj(i) & x_ends == x_obj(i));  
        y_ends(J) = NaN;
        x_ends(J) = NaN;
        
    end
    
%--- FIND AND CONNECT COMPONENTS FOR EACH END-POINT
    for i = 1:length(y_obj)
        % obtain distances between considered end-point and all end-points
        dist = sqrt(abs(x_obj(i) - x_ends).^2 + abs(y_obj(i) - y_ends).^2);  
        % obtain angles between above
        angle = atand(-(y_obj(i) - y_ends)./(x_obj(i) - x_ends)) - STATS.Orientation(I(j));
        % apply filter based on angle and gap distance
        K = find(dist <= gap_length2 & abs(angle) <= angle_var);
        
        if isempty(K)
            
            continue
            
        elseif length(K) >= 2
            % if multiple K, sort on shortest distance and only use that one
            [~,Y] = sort(dist(K));
            K = K(Y(1));
            
        end
        % x-distance between considered and target end-point
        dx = x_obj(i) - x_ends(K); 
        % y-distance between considered and target end-point
        dy = y_obj(i) - y_ends(K);                                          
        
% ------- GET X,Y-COORDINATES OF POINTS TO BE FILLED
        
        % for x-distance larger than y-distance
        if abs(dx) >= abs(dy)  
            % x-indices
            x_fill = 0:1:(abs(dx));
            
            % for target end-point 'in front of' considered end-point
            if dx < 0                                                       
                y_fill = x_fill .* (dy./dx);
                x_fill = x_fill + x_obj(i);
                y_fill = y_fill + y_obj(i);
            % for target end-point 'behind' considered end-point
            else                                                            
                y_fill = x_fill .* (dy./dx);
                x_fill = x_fill + x_ends(K);
                y_fill = y_fill + y_ends(K);
            end
            
            % make integer y-indices
            y_fill = round(y_fill);
            
        % for y-distance larger than x-distance    
        else                                                                
            y_fill = 0:1:(abs(dy)); 
            
            % y-indices
            if dy < 0 
            % for target end-point 'in front of' considered end-point    
                x_fill = -(y_fill .* (dx./dy));
                y_fill = y_fill + y_obj(i);
                y_fill = flip(y_fill);
                x_fill = x_fill + x_ends(K);
            % for target end-point 'behind' considered end-point    
            else                                                            
                x_fill = -(y_fill .* (dx./dy));
                y_fill = y_fill + y_ends(K);
                y_fill = flip(y_fill);
                x_fill = x_fill + x_obj(i);
            end
            % make integer x-indices
            x_fill = round(x_fill);                                     
        end
        
        % fill connection on image
        fill = sub2ind(size(im_bin),y_fill,x_fill);            
        
        clearvars y_fill x_fill
        im_bin_conn(fill) = 1;
    end
end

%% ==================== FILTER CIRCULAR COMPONENTS =======================
im_bin_conn = bwmorph(im_bin_conn,'thin', Inf);                        % skeletonize binary image

% get new stats
STATS   = regionprops(im_bin_conn,'MajorAxisLength','Eccentricity');        
STATS   = struct2table(STATS);
% get indices new connected components
CC = bwconncomp(im_bin_conn); 
% find objects smaller than frac_length .* 2
I = find(STATS.MajorAxisLength <= frac_length.*2& STATS.Eccentricity < ecc_threshold);                            

% remove smaller object
for i = 1:length(I)
    im_bin_conn(CC.PixelIdxList{I(i)}) = 0;
end
figure; 
imshow(im_bin_conn);
set(gca, 'FontSize', 14);
title('Skeletonised Image');
axis on;
axis tight;
axis image;
box on;
xlabel('X, px');
ylabel('Y, px');
%% ============================= PRUNING =================================
%close all;

disp('pruning...');
tic;
% EROSION
BRANCH = bwmorph(im_bin_conn,'branchpoints');

% get end-points for connected components
for i = 1:pruning_length
    ENDS= bwmorph(im_bin_conn, 'endpoints');               
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    END_OLD{i} = ENDS;
    im_bin_conn = im_bin_conn - ENDS;
    im_bin_conn = logical(im_bin_conn);
end

% DILATION
nhood = ones(3,3);
% get end-points for connected components
for i = 1:pruning_length
    ENDS= bwmorph(im_bin_conn, 'endpoints');    
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    ENDS = imdilate(ENDS,nhood);
    ADD = ENDS & END_OLD{(pruning_length+1 - i)};
    im_bin_conn = im_bin_conn + ADD;
    im_bin_conn = logical(im_bin_conn);
end
toc;

x = 1;
% get new stats
STATS= regionprops(im_bin_conn,'Eccentricity','MajorAxisLength','Orientation','Centroid'); 
STATS = struct2table(STATS);
% get indices new connected components
CC = bwconncomp(im_bin_conn);  
% find objects smaller than frac_length .* x
I = find(STATS.MajorAxisLength <= frac_length2);   
% remove smaller objects
for i = 1:length(I)
    im_bin_conn(CC.PixelIdxList{I(i)}) = 0;   
end


% x,y-coordinates of end-points
ENDS_f = bwmorph(im_bin_conn, 'endpoints');
[y_final,x_final] = find(ENDS_f == 1);     

%% ============================== FIGURE =================================
[B,L]   = bwboundaries(im_bin_conn);

% figure; 
% hold on;
figure; 
% set(gcf, 'PaperPositionMode', 'manual') ; 
% set(gcf, 'PaperUnits', 'centimeters') ; 
% set(gcf, 'PaperPosition', [ 1 1 20 20 ]) ; 

% Print the original image together with the outputs produced in this code
subplot(1,2,1); imshow(im_original,[]);
axis on;
axis image;
axis tight;
box on;
set(gca, 'FontSize', 14);
title('Input image');
xlabel('X [px]');
ylabel('Y [px]');

subplot(1,2,2);
imshow(im_original);
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end


% conv = [' conversion px/unit ',num2str(conversion)];
% disp(conv);
% addMM=@(x) sprintf('%.1f',x/conversion);
% xticklabels(cellfun(addMM,num2cell(xticks'),'UniformOutput',false));
% yticklabels(cellfun(addMM,num2cell(yticks'),'UniformOutput',false));
axis on;
box on;
axis tight;
set(gca, 'FontSize',14);
xlabel('X [px]');
ylabel('Y [px]');
title({'Segmented fractures, N = ', num2str(length(B))});
hold off;

%save the segmented fracture network
%print('-djpeg', '-r300', 'FileName_Network.jpeg');

figure;

%plot histogram of fracture length distribution

length_fr = cat(1,STATS.MajorAxisLength);
subplot(1,2,1);
h1 = histogram(length_fr,'FaceColor','r','FaceAlpha',.3);
box on
grid on
axis square;
set(gca, 'FontSize', 14);
xlabel('Fracture length [px]');
ylabel('Frequency');
% caption = sprintf({'Histogram of'; 'fracture distribution'});
title({'Hist of fracture distribution, N = ', num2str(length(B))});

%Save fracture length for statistical analysis
%dlmwrite('FileName_FractureLength.txt', length_fr); % <--- Change "FileName" with current data name

%plot equal area rose diagram of fracture orientation
orient = cat(1, 90-STATS.Orientation);

orient2 = [orient; orient + 180];
for i = 1:max(size(orient2))
    if orient2(i)<0
        orient2(i) = orient2(i) +360;
    end
end
   
subplot(1,2,2) ;
roseEqualArea(orient2, 18, 0, 0, 'red') ;
set(gca, 'FontSize', 14);
title({'Fracture';'orientations'});


%save the histogram + rose diagram of the segmented fracture network 
%print('-djpeg', '-r300', 'FileName_Hist-Rose.jpeg');

figure;
for k = 1:length(B)
    boundary2 = B{k};
    plot(boundary2(:,2), boundary2(:,1),'k', 'LineWidth',1);
    hold on;
end
hold off;
axis on;
axis ij;
axis tight;
box on;
set(gca, 'FontSize', 14)
xlabel('X [px]');
ylabel('Y [px]');
title({'Segmented fractures, N = ', num2str(length(B))});
%save the segmented fracture network to file
%print('-djpeg', '-r300', 'FileName_Segmented.jpeg');

%% ======================== PRINT DATA TO FILE ===========================
% Save coordinates of segmented fracture into a txt for FracPaQ
% Reshape fracture coordinates in cell array B, into a format readable by
% FracPaQ: each line in the .txt file correspond to a segmented fracture
% with pairs of xn - yn coordinates.

for q = 1:length(B)
   
    bb2 = B{q,1}; 

% % Ramer-Douglas-Peucker algorithm for curve simplification (https://en.wikipedia.org/wiki/Ramer?Douglas?Peucker_algorithm)
   bb2 = RDPsimplify(bb2,1);
  
    if size(bb2,1) > 1
        k = 0; ss = zeros(1,2*length(bb2)); 
        for i = 1:length(bb2)
            ss(i+k) = bb2(i,2); 
            ss(i+k+1) = bb2(i,1);
            k = k+1; 
        end
       % dlmwrite('FracPaQ-Input_FileName.txt',ss,'delimiter','\t','-append');
    end
 end
