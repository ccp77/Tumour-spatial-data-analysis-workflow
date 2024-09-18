%% This script aims at generating cell densities from cell clouds, and to output masks based on thresholded densities.
clear

file_directory = 'path/';
single_image_sub_directory = 'sample/';

vessel_csv = 'data1.csv'; 
CD8_csv = 'data2.csv'; 
tumour_csv = 'data3.csv';
pixel_scale=1;


input_file_directory = [file_directory,single_image_sub_directory,'/'];
output_file_directory = [input_file_directory,'output_masks/'];
epsilon = 1e-20;


if ~isfolder(output_file_directory)
    mkdir(output_file_directory);
    disp(['Created directory  ', output_file_directory]);
else
    disp(['Directory already exists: ', output_file_directory]);
end

point_cloud= readtable([file_directory,single_image_sub_directory,vessel_csv]); % create the mask based on the TCF7 coordinates

szy = ceil(max(point_cloud.Y) + 1000);
szx = ceil(max(point_cloud.X) + 1000);
mask = zeros(szy,szx);
gridx = 0:min(szx,szy)/1000:szx;
gridy = 0:min(szx,szy)/1000:szy;
[x,y] = meshgrid(gridx, gridy);
whos gridx
xa = x(:);
ya = y(:);
xi = [xa ya];
mask_resize = imresize(mask,size(x));

[vessel,f_vessel_2d] = output_density_shift_points_by_500(input_file_directory, vessel_csv, pixel_scale, xi,szx,szy,x,mask_resize,gridx,gridy);
[CD8,f_CD8_2d] = output_density_shift_points_by_500(input_file_directory, CD8_csv, pixel_scale, xi,szx,szy,x,mask_resize,gridx,gridy);


%% Generate mask based on density
% for the vessel: generate masks using two threshold, wither 0.3 or 0.5
f_vessel_2d_resized = imresize(f_vessel_2d,size(mask));

figure, imagesc(f_vessel_2d)
colorbar()
figure, imagesc(f_vessel_2d_resized)
colorbar()

data = f_vessel_2d_resized;
threshold = 0.5*10^-7; % choose the threshold
data(data < threshold) = 0;
data(data >=threshold) = 1;

figure, imagesc(data)
colorbar()

nb_of_row = size(data,1);
idx = find(data >= 1);
Matrix = [ceil(idx/nb_of_row) idx- (ceil(idx/nb_of_row)-1)*nb_of_row]; % first is column index, second is row index 
Matrix_minus_500 = Matrix-500;
T = array2table(Matrix_minus_500, 'VariableNames', {'X', 'Y'});

figure, scatter(T.X,T.Y,2,"red",".")
hold on;
scatter(point_cloud.X,point_cloud.Y,2,"black",".")

% prepare the data for export: downsample points in the mask, and add
% columns if necessary to fin data structure necessary in following
% pipelines (For example here add 6 columns with 0 values and with headers
% matching other csv files that will analysed jointly in a next step.

T_mat = table2array(T);
newTabCol = zeros(height(T_mat), 6);
new = [T_mat newTabCol];
new_table = array2table(new, 'VariableNames', {'X','Y','nCount_RNA', 'nFeature_RNA','RNA_snn_res.1','CCR7','CXCL9','volume'});

new_table_downsample = downsample(new_table,100);

output_filename = [output_file_directory,'mask_vessel_downsample_thresholdDensity0.5.csv'];
writetable(new_table_downsample, output_filename)

% for the CD8
f_CD8_2d_resized = imresize(f_CD8_2d,size(mask));
figure, imagesc(f_CD8_2d)
colorbar()
figure, imagesc(f_CD8_2d_resized)
colorbar()

data_CD8 = f_CD8_2d_resized;
threshold = 0.5*10^-7;
data_CD8(data_CD8 < threshold) = 0;
data_CD8(data_CD8 >=threshold) = 1;

figure, imagesc(data_CD8)
colorbar()

nb_of_row = size(data_CD8,1);
idx = find(data_CD8 >= 1);
Matrix_CD8 = [ceil(idx/nb_of_row) idx- (ceil(idx/nb_of_row)-1)*nb_of_row]; % first is column index, second is row index 
Matrix_CD8_minus_500 = Matrix_CD8-500;
T_CD8 = array2table(Matrix_CD8_minus_500, 'VariableNames', {'X', 'Y'});

newTabCol = zeros(height(T_CD8), 6);
T_CD8_mat = table2array(T_CD8);
new = [T_CD8_mat newTabCol];
new_table = array2table(new, 'VariableNames', {'X','Y','nCount_RNA', 'nFeature_RNA','RNA_snn_res.1','CCR7','CXCL9','volume'});

new_table_downsample = downsample(new_table,100);

output_filename = [output_file_directory,'mask_CD8_downsample_thresholdDensity0.5.csv'];
writetable(new_table_downsample, output_filename)

%% Plot the position of masks and the cells

figure, scatter(T_CD8.X,T_CD8.Y,"blue")
%set(gca, 'YDir','reverse')
hold on;
scatter(T.X,T.Y,"red")
%set(gca, 'YDir','reverse')
hold on;
scatter(point_cloud.X,point_cloud.Y,2,"black",".")


%% Define function

function [point_cloud,f_density_2d] = output_density_shift_points_by_500(file_directory, file_csv, pixel_scale, xi,szx,szy,x,mask,gridx,gridy)
    point_cloud= readtable([file_directory,file_csv]);
    point_cloud.X=point_cloud.X/pixel_scale;
    point_cloud.Y=point_cloud.Y/pixel_scale;
    point_cloud.X=point_cloud.X+500;
    point_cloud.Y=point_cloud.Y+500;
    [f_density,xi]=ksdensity([point_cloud.X,point_cloud.Y],xi,'Bandwidth',100,'Support',[0 0;szx szy],'BoundaryCorrection','reflection');
    f_density_2d=reshape(f_density,size(x));
    %f_density_2d=reshape(f_density,size(x)).*mask;
    A_density =trapz(gridx,trapz(gridy,f_density_2d,1));
    f_density_2d = f_density_2d./A_density;
    disp(['Density completed for ', file_csv]);
end

