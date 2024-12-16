%TITLE:MRI simulation code and phantom data sets
%Matlab script to generate figures using Sobol's GRE simulation
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 10/28/2024
%Email: dtamada@wisc.edu


%%
%Figure 5
%Plot healthy volunteer images

%3T heanlth vol
R2map = dicomread('../phantom_raw_T2_data/Volunteer/volunteer_DICOM_1p5T/image13.dcm');
figure;
imagesc(R2map, [0 75]);
axis square;
r2s_3T = dicomread('../phantom_raw_T2_data/Volunteer/volunteer_R2s_1p5T/I0018.dcm');
figure;
imagesc(r2s_3T, [0 100]);axis square;


%1.5T heanlth vol
R2map = dicomread('../phantom_raw_T2_data/Volunteer/volunteer_DICOM_3T/image14.dcm');
figure;
imagesc(R2map, [0 75]);
axis square;
r2s_1p5T = dicomread('../phantom_raw_T2_data/Volunteer/volunteer_R2s_3T/I0018.dcm');
figure;
imagesc(r2s_1p5T, [0 100]);axis square;



%%
%Code for generating R2 map

%3T
file_3T = '../phantom_raw_T2_data/Volunteer/volunteer_raw_data_3T.mat';
load(file_3T, 'img_recon', 'header', 'params', 'kspace', 'G','filename');
beta_map = 1.0 + 0.0*double(img_recon(:,:,:,1))/100.0;
% imaging parameters
TR      = params.TR;
FA      = params.FA;
dphi    = params.dphi;
R1_guess = 2.0;
field_strength = header.rdb_hdr_exam.magstrength/1e4;
T1 = 1000/R1_guess;

%Define coefficient of relationship between R1 and R2 for in vivo.
if field_strength > 1.5
    coeff = [0.0038 0.89];
else
    coeff = [0.0065 1.23];
end

% Obtain phase
[LUTs] = build_LUTs_vdphi_T1B1(dphi*(pi/180), FA*(pi/180), TR*1000, false, coeff);

kern = ones(3,3,1)/(3*3*1);
img_recon_conv = convn(img_recon,kern,'same');
ph_conv = angle(img_recon_conv(:,:,:,2).*conj(img_recon_conv(:,:,:,1)))/2;

%Calculate T2 map slice-by-slice
T2map = T2map_phase((ph_conv(:,:,16)), LUTs);
R2map = 1./T2map;
figure; imagesc(R2map, [0 75]); title('PB R2 at 3.0T');

%%

%1.5T
file_1p5T = '../phantom_raw_T2_data/Volunteer/volunteer_raw_data_1p5T.mat';
load( file_1p5T, 'img_recon', 'header', 'params', 'kspace', 'G','filename');
beta_map = 1.0 + 0.0*double(img_recon(:,:,:,1))/100.0;

% imaging parameters
TR      = params.TR;
FA      = params.FA;
dphi    = params.dphi;
R1_guess = 2.0;
field_strength = header.rdb_hdr_exam.magstrength/1e4;

T1 = 1000/R1_guess;

%Define coefficient of relationship between R1 and R2 for in vivo.
if field_strength > 1.5
    coeff = [0.0038 0.89];
else
    coeff = [0.0065 1.23];
end

% Obtain phase
[LUTs] = build_LUTs_vdphi_T1B1(dphi*(pi/180), FA*(pi/180), TR*1000, false, coeff);

kern = ones(3,3,1)/(3*3*1);
img_recon_conv = convn(img_recon,kern,'same');
ph_conv = angle(img_recon_conv(:,:,:,2).*conj(img_recon_conv(:,:,:,1)))/2;

%Calculate T2 map slice-by-slice %1+0*
T2map = T2map_phase((ph_conv(:,:,16)), LUTs);
R2map = 1./T2map;
figure; imagesc(R2map, [0 75]); title('PB R2 at 1.5T');

