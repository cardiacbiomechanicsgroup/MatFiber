% MatFiber fiber orientation analysis code.
% Source: 
%     Cardiac Biomechanics Group
%     PI: Jeffrey W. Holmes, MD,PhD
%         Biomedical Engineering, Medicine, Robert M. Berne Cardiovascular Research Center
%         University of Virginia
%     If used, please reference Fomovsky & Holmes. "Evolution of scar structure, mechanics, and ventricular 
%     function after myocardial infarction in the rat." AJP Heart Circ Phys 298: H221-228, 2010.

% Description:
%     MatFiber is a Matlab adaptation of Fiber3, a C++ program that uses the intensity gradient method described by 
%     Karlon, Covell, McCulloch, Hunter, and Omens. The Anatomical Record 252:612-625. 1998.

clear all;

% SPECIFY IMAGE FILE (8-BIT GRAYSCALE), SUBREGION SIZE, AND THRESHOLD LEVEL

% SUBREGION SIZE (S) determines how many pixels (S x S) you wish to group into only 1 output angle.
% The code will divide the original image into N1 x N2 subregions of SxS
% pixels each, and calculate 1 orientation value for each subregion. This
% orientation will correspond to the angle perpendicular to the strongest 
% mean intensity gradient direction within that subregion.

% THRESHOLD (Thresh) determines the minimum value of intensity gradient
% within a subregion that is required for the code to output an orientation
% for that subregion. Subregions with summed intensity gradients below 'Thresh' 
% will not return an orientation angle. Warning: because the threshold is compared to intensity
% gradients summed within each subregion of SxS pixels, an appropriate value of 'Thresh' may vary depending upon the
% subregion size S.

% Sample images:
% 'SIMULATED1' ... this is an artifically generated image of straight lines, populated from a distribution with MVL = 0.3 and Mean Angle = 0.
% 'SIMULATED2' ... this is an artifically generated image of straight lines, populated from a distribution with MVL = 0.7 and Mean Angle = -60.
% 'CP15S219_SUBTRACTED' ... subtracted image of infarct scar (processing methods descibed in Fomovsky 2010)
% Also included are output images after analyzing 'CP15S219_SUBTRACTED'
% with settings we considered appropriate for this image (S = 40, Thresh = 100000),
% as well as with substantially higher/lower subregion size, or thresholds.
% Note that increasing subregion size and threshold resulted in an increased
% mean vector length (MVL) measurement, while decreasing subregion size and
% threshold resulted in decreased MVL.

filename = 'CP15S219_SUBTRACTED';
S = 40;                                      % Size of the square subregion (pixels)
Thresh = 100000;
iImage = imread([filename '.tif'], 'tif');   % Import the image
dImage = double(iImage)+1;                   % Convert into double

% Initial computations related to the image and subregions
% Because gradient mask (next step) is 9x9, convolution result will be
% erroneous in pixels near image border. Subregions are shifted 5 pixels
% away from the image borders to avoid this error.
D = size(dImage);
N1 = floor((D(1)-10)/S);           % Number of subregions vertically;
N2 = floor((D(2)-10)/S);           % Number of subregions horizontally;

% Gradient masks computation
% for i=-4:4;
%     for j=-4:4;
%         MX(i+5,j+5) = 2*i/9*exp(-(i^2+j^2)/4);
%         MY(i+5,j+5) = 2*j/9*exp(-(i^2+j^2)/4);
%     end
% end
% Hard coded here to save computation time.
MX = ...
    [-0.00029819    -0.001716	-0.0059893	-0.012679	-0.016281	-0.012679	-0.0059893	-0.001716	-0.00029819;...
    -0.001287      -0.007406	-0.025849	-0.054723	-0.070266	-0.054723	-0.025849	-0.007406	-0.001287;...
    -0.0029946     -0.017233	-0.060149	-0.12734	-0.1635     -0.12734	-0.060149	-0.017233	-0.0029946;...
    -0.0031698     -0.018241	-0.063668	-0.13478	-0.17307	-0.13478	-0.063668	-0.018241	-0.0031698;...
    0              0           0           0           0           0           0           0           0;...
    0.0031698      0.018241	0.063668	0.13478     0.17307     0.13478     0.063668	0.018241	0.0031698;...
    0.0029946      0.017233	0.060149	0.12734     0.1635      0.12734     0.060149	0.017233	0.0029946;...
    0.001287       0.007406	0.025849	0.054723	0.070266	0.054723	0.025849	0.007406	0.001287;...
    0.00029819     0.001716	0.0059893	0.012679	0.016281	0.012679	0.0059893	0.001716	0.00029819];
MY = ...
    [-0.00029819    -0.001287	-0.0029946	-0.0031698	0	0.0031698	0.0029946	0.001287	0.00029819;...
    -0.001716      -0.007406	-0.017233	-0.018241	0	0.018241	0.017233	0.007406	0.001716;...
    -0.0059893     -0.025849	-0.060149	-0.063668	0	0.063668	0.060149	0.025849	0.0059893;...
    -0.012679      -0.054723	-0.12734	-0.13478	0	0.13478     0.12734     0.054723	0.012679;...
    -0.016281      -0.070266	-0.1635     -0.17307	0	0.17307     0.1635      0.070266	0.016281;...
    -0.012679      -0.054723	-0.12734	-0.13478	0	0.13478     0.12734     0.054723	0.012679;...
    -0.0059893     -0.025849	-0.060149	-0.063668	0	0.063668	0.060149	0.025849	0.0059893;...
    -0.001716      -0.007406	-0.017233	-0.018241	0	0.018241	0.017233	0.007406	0.001716;...
    -0.00029819	-0.001287	-0.0029946	-0.0031698	0	0.0031698	0.0029946	0.001287	0.00029819];

% Convolve masks with the image
GX = conv2(dImage,MX,'same');
GY = conv2(dImage,MY,'same');
clear MX MY

% Edge image and gradient direction
E = GX.^2+GY.^2;
phi = 180/pi*atan2(GY, GX);
clear GX GY

% Determine local orientation in each subregion
bins = 0:1:179;
FiberAngle = zeros([N1*N2 1]);
FiberPosX = zeros([N1*N2 1]);
FiberPosY = zeros([N1*N2 1]);
count = 0;

for q = 1:1:N1
    for r = 1:1:N2
        lx = 5 + S*(r - 1) + 1;
        ux = 5 + S*r;
        ly = 5 + S*(q - 1) + 1;
        uy = 5 + S*q;
        AthetaW = zeros(size(bins));
        for n = 1:1:numel(bins)
            C = E(ly:uy,lx:ux).*exp(2*cos(2*pi/180*(bins(n) - phi(ly:uy,lx:ux))))/exp(2);
            AthetaW(n) = sum(sum(C)); % AthetaW is a sum of alignment values between pixel gradient directions and each direction from 'bins', weighted by pixel gradient magnitudes.
                                      % The sum includes all pixels within the subregion bounded by lx, ux, ly, uy.
        end
        if max(AthetaW) > Thresh;
            [~, Ang] = max(AthetaW); % Sets 'Ang' as the entry from 'bins' corresponding to the maximum 'AthetaW' value.
            if Ang > 90
                Ang = Ang - 180;
            end
            count = count + 1;
            FiberAngle(count) = Ang;
            FiberPosY(count) = round(5 + S*q - S/2);
            FiberPosX(count) = round(5 + S*r - S/2);
        end
        clear AthetaW Ang;
    end
end

if count < numel(FiberAngle)
    FiberAngle(count+1:end) = [];
    FiberPosY(count+1:end) = [];
    FiberPosX(count+1:end) = [];
end

FiberOrientationX = cos(FiberAngle*pi/180);
FiberOrientationY = sin(FiberAngle*pi/180);
FiberPosY = -(FiberPosY-D(1));

% Calculate MA, MVL, CSD
c2m = mean(cos(2*FiberAngle*pi/180));
s2m = mean(sin(2*FiberAngle*pi/180));
MA = 180/pi*1/2*atan2(s2m, c2m);        % mean angle (deg)
MVL = sqrt(s2m^2 + c2m^2);              % mean vector length
CSD = 180/pi*1/2*sqrt(-2*log(MVL));     % circular standard deviation (deg)
AD = 180/pi*1/2*sqrt(2*(1-MVL));        % angular deviation (deg)

% Display Results
figure;
clf;
imagesc(iImage);
colormap(gray);
hold on;
quiver(FiberPosX, -FiberPosY+D(1), FiberOrientationX, -FiberOrientationY, 0.5, 'y', 'LineWidth', 1);
hold off;
title(['MatFiber: Thresh = ' num2str(Thresh) ', Size = ' num2str(S) ', MA = ' num2str(MA), ', MVL = ' num2str(MVL) ', CSD = ' num2str(CSD)]);
set(gca, 'DataAspectRatio', [1 1 1]);
display(filename);
display([Thresh S]);
display([MA MVL CSD]);

% Save results
clear phi N1 N2 E C;
% save([filename '_' num2str(S) 'px_FiberAngles.mat'], '-mat');
