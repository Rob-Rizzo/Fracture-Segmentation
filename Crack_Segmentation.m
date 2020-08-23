% AUTO-FRACTURE TRACING 
% Fracture tracing script, currently only tested on SEM/BSE images, small size photographs and CT micrographs 
%
% Most important parameters to changes: MEDIAN FILTER aperture, BINARIZED bw_threshold, MORPHOLOGICAL PARAMETER gap_length1, conversion. 
% 
% Author: Roberto Rizzo @ HWU Edinburgh / University of Aberdeen
% Date: February 2020
% Based on Luke Griffiths (crack-extract Python Code) and Franz Abes code (FRAg)

clear;
close all;
workspace;
format short g;
format compact;

%% --------------------------- SET PARAMETERS ----------------------------
%INPUT Value of pixel per mm in the analysed image
conversion=700;

%MEDIAN FILTER
% n x m area for median filter. preferably smaller than the pore size and larger than the fracture aperature
aperture1      = [32 32];  % default [24 24]
aperture2      = [2 2]; %default [3 3]

% BINARIZATION AND FILTER PARAMETERS
% binarization threshold (1-100)
bw_threshold   = 30;  %default 30
% neighbourhoodsize (x,y) for wiener filter
wiener_area1    = [10 10];    
wiener_area2    = [3 3]; 
% min-max range to rescale grayscale values
rescale_range  = [0 40];   

% MORPHOLOGICAL PARAMETERS
% minimum fracture length in pixels (given by the major axes length,
% regionprops function)
frac_length     = 8;    % default 8
% largest object size (nr of pixels) to be removed
object_limit    = round((2/3).*frac_length);  
% gap length closed by morphological dilation-erosion
gap_length1     = round((1/3).*frac_length);   

% TARGETED CLOSING
% eccentricity threshold (0 is circular, 1 is line)
ecc_threshold   = 0.90;
% variation in alignment between major axis of component and end-point of a second component (-angle_var to + angle_var, in degrees)
angle_var       = 15;
% gap length between line segments closed by targeted connections
gap_length2     = 10;  

% PRUNING
% length of branches to be pruned
pruning_length  = frac_length; 
% last size-based filter
frac_length2 = 1.*frac_length;     

%% ============================ READ DATA ================================
[FILENAME, PATHNAME] = uigetfile({'*.*'},'Select an image');
Im_original = imread([PATHNAME FILENAME]);
%Im_original = imread("Cor_Pre.png");

figure; imshow(Im_original,[]);
axis on;
axis image;
axis tight;
box on;
caption = sprintf('Input Image');
title(caption, 'Interpreter', 'None');
xlabel('X, px');
ylabel('Y, px');


% Get the dimensions of the image.
% numberOfColorChannels should be = 1 for a gray scale image, and 3 for an RGB color image.
[rows, columns, numberOfColorChannels] = size(Im_original);
if numberOfColorChannels > 1
  % It's not really gray scale like we expected - it's color.
  % Use weighted sum of ALL channels to create a gray scale image.
  Im_original = rgb2gray(Im_original);
  % ALTERNATE METHOD: Convert it to gray scale by taking only the green channel,
  % which in a typical snapshot will be the least noisy channel.
  % grayImage = grayImage(:, :, 2); % Take green channel.
end

% SETTINGS FOR CUTTING DOWN SIZE FOR SPECIFIC IMAGE SIZES
% Im_original = Im_original(1:1696,1:end);
% Im_original = Im_original(1:849,1:end);

%%========================== BINARIZATION ===============================
%FILTER1 - Filter out small components
tic;
%median filter
Im_filt = medfilt2(Im_original, aperture1,'symmetric');
% difference image
Im_filt = Im_filt - Im_original;
% remove noise
Im_filt = wiener2(Im_filt, wiener_area1);                                    
toc;

% FILTER2 - filter out large components
tic;
% median filter
Im_filt2 = medfilt2(Im_original, aperture2,'symmetric'); 
% difference image
Im_filt2 = Im_filt2 - Im_original;
% remove noise
Im_filt2 = wiener2(Im_filt2, wiener_area2);  
toc;

% Final filtered image
Im_filt = Im_filt+Im_filt2;

% BINARIZE
tic;
% rescale grayscale values
Im_bin = mat2gray(Im_filt,rescale_range);      
% set values to zero, based on non-rescaled image
Im_bin(Im_filt < bw_threshold)= 0;   
% set values to one, based on non-rescaled image
Im_bin(Im_filt >= bw_threshold)= 1;
% binarise image (convert to logical array)
Im_bin = imbinarize(Im_bin);                            
toc;
%display figure
figure; imshowpair(Im_filt, Im_bin, 'montage'); title('Filtred vs. Binarized Images');
axis on;
axis image;
axis tight;
box on;
xlabel('X, px');
ylabel('Y, px');


% remove smallest objects
Im_bin = bwareaopen(Im_bin,object_limit);   
tic;

% angle of the line object
for alpha = 1:181  
    % set orientation of line object
    se = strel('line' ,gap_length1 ,alpha);     
    % close gaps between line segments by dilation and erosion
    Im_bin = imclose(Im_bin, se);                                          
end

 % skeletonize binary image
Im_bin= bwmorph(Im_bin,'thin', Inf);                                 
toc;

figure; imshowpair(Im_filt,Im_bin,'montage'); title('Filtred vs. Skeletonized Images');
axis on;
axis tight;
axis image;
box on;
xlabel('X, px');
ylabel('Y, px');

%% =================== CONNECTION OF LINEAR OBJECTS ======================
% ------------------------ MORPHOLOGICAL CLOSING ------------------------
% ------------------------- TARGETED CLOSING -----------------------------

% GET CONNECTED COMPONENT PROPERTIES
% get properties of connected components
STATS = regionprops(Im_bin,'Eccentricity','MajorAxisLength','Orientation','Centroid');                         
STATS = struct2table(STATS);
% get indices for all connected components
CC = bwconncomp(Im_bin); 
% get end-points for ALL connected components
ENDS = bwmorph(Im_bin, 'endpoints');    
% x,y-coordinates of end-points
[y_ends,x_ends] = find(ENDS == 1);                                          

% filter connected components
I = find(STATS.MajorAxisLength >= frac_length & STATS.Eccentricity > ecc_threshold);                    

% set image that will be connected
Im_bin_conn = Im_bin;                                                  

disp(['Found ' num2str(length(I)) ' suited connected objects']);

%--- LOOP OVER FILTERED CONNECTED COMPONENTS-------
for j = 1:length(I)
    
   % disp(['analysing object ' num2str(j) ' of ' num2str(length(I))]);
    
    Im_work = Im_bin;
    % remove component under consideration from work image
    Im_work(CC.PixelIdxList{I(j)}) = 0;  
    % get end-notes of component under consideration
    ENDS_obj  = bwmorph(Im_bin - Im_work, 'endpoints');
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
        fill = sub2ind(size(Im_bin),y_fill,x_fill);            
        
        clearvars y_fill x_fill
        Im_bin_conn(fill) = 1;
    end
end

%% ==================== FILTER CIRCULAR COMPONENTS =======================
% Im_bin_conn      = bwmorph(Im_bin_conn,'thin', Inf);                        % skeletonize binary image

% get new stats
STATS   = regionprops(Im_bin_conn,'MajorAxisLength','Eccentricity');        
STATS   = struct2table(STATS);
% get indices new connected components
CC = bwconncomp(Im_bin_conn); 
% find objects smaller than frac_length .* 2
I = find(STATS.MajorAxisLength <= frac_length.*2& STATS.Eccentricity < ecc_threshold);                            

% remove smaller object
for i = 1:length(I)
    Im_bin_conn(CC.PixelIdxList{I(i)}) = 0;
end
figure; imshow(Im_bin_conn);title('Skeletonised Image');
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
BRANCH = bwmorph(Im_bin_conn,'branchpoints');

% get end-points for connected components
for i = 1:pruning_length
    ENDS= bwmorph(Im_bin_conn, 'endpoints');               
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    END_OLD{i} = ENDS;
    Im_bin_conn = Im_bin_conn - ENDS;
    Im_bin_conn = logical(Im_bin_conn);
end

% DILATION
nhood = ones(3,3);
% get end-points for connected components
for i = 1:pruning_length
    ENDS= bwmorph(Im_bin_conn, 'endpoints');    
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    ENDS = imdilate(ENDS,nhood);
    ADD = ENDS & END_OLD{(pruning_length+1 - i)};
    Im_bin_conn = Im_bin_conn + ADD;
    Im_bin_conn = logical(Im_bin_conn);
end
toc;

x = 1;
% get new stats
STATS= regionprops(Im_bin_conn,'Eccentricity','MajorAxisLength','Orientation','Centroid'); 
STATS = struct2table(STATS);
% get indices new connected components
CC = bwconncomp(Im_bin_conn);  
% find objects smaller than frac_length .* x
I = find(STATS.MajorAxisLength <= frac_length2);   
% remove smaller objects
for i = 1:length(I)
    Im_bin_conn(CC.PixelIdxList{I(i)}) = 0;   
end


% x,y-coordinates of end-points
ENDS_f = bwmorph(Im_bin_conn, 'endpoints');
[y_final,x_final] = find(ENDS_f == 1);     

%% ============================== FIGURE =================================
[B,L]   = bwboundaries(Im_bin_conn);

% figure; 
% hold on;
figure; 
% set(gcf, 'PaperPositionMode', 'manual') ; 
% set(gcf, 'PaperUnits', 'centimeters') ; 
% set(gcf, 'PaperPosition', [ 1 1 20 20 ]) ; 

% Print the original image together with the outputs produced in this code
subplot(1,2,1); imshow(Im_original,[]);
axis on;
axis image;
axis tight;
box on;
title('Input image');
xlabel('X, px');
ylabel('Y, px');

subplot(1,2,2);
imshow(Im_original);
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end


conv = [' conversion px/unit ',num2str(conversion)];
disp(conv);
addMM=@(x) sprintf('%.1f',x/conversion);
xticklabels(cellfun(addMM,num2cell(xticks'),'UniformOutput',false));
yticklabels(cellfun(addMM,num2cell(yticks'),'UniformOutput',false));
axis on;
box on;
axis tight;
xlabel('X, mm');
ylabel('Y, mm');
title({'Segmented fracture network, N = ', num2str(length(B))});
hold off;

%save the fracture hist and rose diagram
print('-djpeg', '-r300', 'TEST_ceramics_Segmented.jpeg');

figure;

%plot histogram of fracture length distribution
length_fr = cat(1,STATS.MajorAxisLength./conversion);
subplot(1,2,1);
h1 = histogram(length_fr,'FaceColor','r','FaceAlpha',.3);
box on
grid on
axis square;
xlabel('Fracture length, mm');
ylabel('Frequency');
% caption = sprintf({'Histogram of'; 'fracture distribution'});
title({'Histogram of'; 'fracture distribution'});

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
title({'Fracture';'orientations'});


%save the original image + segmented,and histogram + rose diagram fracture network to file
print('-djpeg', '-r300', 'TEST_ceramics_Hist-Rose.jpeg');

figure;
for k = 1:length(B)
    boundary2 = B{k};
    plot(boundary2(:,2), boundary2(:,1),'k', 'LineWidth',1);
    hold on;
end
hold off;
axis on;
% axis ij;
axis tight;
box on;
xlabel('X, px');
ylabel('Y, px');
title({'Segmented fractures, N = ', num2str(length(B))});
%save the segmented fracture network to file
print('-djpeg', '-r300', 'TEST_ceramics_Network.jpeg');

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
        dlmwrite('TEST_ceramics.txt',ss,'delimiter','\t','-append');
    end
 end
