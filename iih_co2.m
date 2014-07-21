%  Perform IIH correction over MRI images using CO2 method. 
%  
%  It works with nifti toolbox available in: http://www.rotman-baycrest.on.ca/~jimmy/NIfTI/#update 
%  
%  Usage: corrected = iih_co2(original_path,linear_path,sup_path,[save_path],[interp_method])
%
%  Where:
%
%	original_path:  Path to the image file in nifti format. *.hdr or
%	*.nii. It can be directly a matrix variable or a nii structure (see
%	load_nii)
%
%	linear_path:    Path to the image file PD which is acquired with linear
%	coil. It can be directly a matrix variable or a nii structure (see
%	load_nii)
%
%   sup_path:    Path to the image file PD which is acquired with
%   superficial coil. It can be directly a matrix variable or a nii
%   structure (see load_nii)
%
%	save_path (optional):	Path where to save result as nifti format.
%
%	interp_method (optional):	'linear' as default. You can choose: 'cubic',
%	'spline' or 'nearest'. (see interp3). Is the interpolation method for
%	IIH estimation.
%
%	
%
%  e.g.:
%       nii_original = load_nii('./Images/Originales/Original.hdr', [], 1);
%       nii_lin = load_nii('./Images/Body/Body.hdr', [], 1);
%       nii_sup = load_nii('./Images/Sup/Sup.hdr', [], 1);
%       
%       corrected_img_linear=
%       iih_co2(nii_original.img,nii_lin.img,nii_sup.img,'./Images/Originales/Original_correctedIIH.hdr','linear');
%
%   NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - E. Marinetto (nenetto@gmail.com)
%



function corrected_img = iih_co2(original_path,linear_path,sup_path,save_path,interp_method)
if nargin < 5
    interp_method = 'linear';
end


%% Carga de los datos

if(isstruct(original_path))
    nii_original = original_path;
elseif(ischar(original_path))
    nii_original = load_nii(original_path, [], 1);
else
    nii_original = make_nii(original_path);
end

if(isstruct(linear_path))
    nii_lin = linear_path;
elseif(ischar(linear_path))
    nii_lin = load_nii(linear_path, [], 1);
else
    nii_lin = make_nii(linear_path);
end

if(isstruct(sup_path))
    nii_sup = sup_path;
elseif(ischar(sup_path))
    nii_sup = load_nii(sup_path, [], 1);
else
    nii_sup = make_nii(sup_path);
end


if(nii_lin.hdr.dime.dim ~= nii_sup.hdr.dime.dim)
    error('Images Linear and Sup has not the same dimmensions.'); 
elseif (nii_lin.hdr.dime.pixdim ~= nii_sup.hdr.dime.pixdim)    
    error('Images Linear and Sup has not the same pixel dimmensions.');    
end

dim = nii_lin.hdr.dime.dim(2:4);
DIM = nii_original.hdr.dime.dim(2:4);

% Preparado de los datos y reformateo
lineal = double(reshape(nii_lin.img,numel(nii_lin.img),1));
sup = double(reshape(nii_sup.img,numel(nii_sup.img),1));
clear('nii_lin','nii_sup');

% Estimación inicial del IIH
iih_est = lineal(sup~=0)./sup(sup~=0);
mean_iih = mean(iih_est(sup~=0));
iih_est(sup==0) = mean_iih;
iih_est3D = reshape(iih_est,dim);
clear('lineal','sup','mean_iih','iih_est');

%% Filtrado de mediana
iih_est3D = medfilt3(iih_est3D);

%% Redimensionado

[XI,YI,ZI] = meshgrid(interp(1:dim(1),DIM(1)/dim(1)),interp(1:dim(2),DIM(2)/dim(2)),interp(1:dim(3),DIM(3)/dim(3)));
iih_est3D = interp3(iih_est3D,XI,YI,ZI,interp_method,0);
clear('XI','YI','ZI','dim');

%% Normalización a la media del IIH

iih_est3D = reshape(iih_est3D,numel(iih_est3D),1);
mean_iih = mean(iih_est3D);
iih_est3D = iih_est3D/mean_iih;
clear('mean_iih');

%% Corrección CO2
original = double(reshape(nii_original.img,numel(nii_original.img),1));
MAX_ori = max(original);
MIN_ori = min(original);
corrected = original;
corrected(iih_est3D~=0) = original(iih_est3D~=0).*iih_est3D(iih_est3D~=0);
clear('original');

%% Normalización
corrected = corrected/mean(corrected);
corrected = corrected/max(corrected);
corrected = corrected*(MAX_ori-MIN_ori)+MIN_ori;

corrected = reshape(corrected,DIM(1),DIM(2),DIM(3));
nii_corrected = nii_original;
nii_corrected.img = corrected;
nii_corrected.hdr.dime.datatype = 64;
clear('DIM','original','MAX_ori','MIN_ori');

%% Saving

if(isstruct(save_path))
    save_path = nii_corrected;
    corrected_img = nii_corrected;
elseif(ischar(save_path))
    save_nii(nii_corrected,save_path);
    corrected_img = nii_corrected;
else
    corrected_img = nii_corrected.img;
end

end