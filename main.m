clc;
clear;
global x_array
global y_array
global lambda_ele
global h_element
global space
global Nsteps
global lambda
%%
%   the size of metasurface
%                     x_array
%                |←          →|
%               ---
%               ↑
%         y_array
%               ↓
%               ---
%
%   the order of metasurface units
%
%               1→2→3→
%                        ↓
%               ←←←←←
%               ↓
%                →4→5→6
%
%   the order of the boundary
%                    list2
%                 -----------
%                 |         |
%           list1 |         | list3
%                 |         |
%                 -----------
%                   list4
%%
x_array = 9;
y_array = 9;
%select current folder(the root directory of this program)
selpath = uigetdir;
% parameters set
h_element = '800[nm]';%The unit must be enclosed in square brackets!(only support m/ mm/ um/ nm)
% si substrate's thickness
h_substrate = '0.3[um]';%The unit must be enclosed in square brackets!!!!!
% air's thickness
h_air = '0.5[um]';
% element's period
lambda_ele = '800[nm]';
lambda = '1550[nm]';
% margin
space = '300[nm]';
% input beam's type: gaussian/plane/dipole(lowercase)
%the refractive index of nano structure
index_nano = 3.478;

%transfer string to number
h_element = str2data(h_element);
h_substrate = str2data(h_substrate);
lambda_ele = str2data(lambda_ele);
space = str2data(space);
h_air = str2data(h_air);
lambda = str2data(lambda); 
%--------------------------------------------------------------------------

% init elements
param_design(1:x_array * y_array) = lambda_ele*0.5;
% param_design=importdata('param.txt');
% load('param');
% param_design=radius_mask';
% param_design = rand (1, x_array * y_array) * lambda_ele * 0.8;
param_design(param_design <= 100e-9) = 100e-9;
param_design(param_design >= lambda_ele) = lambda_ele*0.95; 
%--------------------------------------------------------------------------

%flag for forward(0) and backward simulation(1)
flag = 0;

%generate adjoint field(metalens)
gene_E;

Nsteps=500;
step_size= 2.5e-2;
balance = 1;%optional
%adam parameter(default)
mopt = zeros(1,x_array * y_array);
vopt = zeros(1,x_array * y_array);

%if the folder isn't exist, creat
check_folder('electrc_field_profile_flag_0', ...
             'electrc_field_profile_flag_1', 'FOM');
if exist('FOM_list.txt', 'file')
    delete('FOM_list.txt');
end
if exist('para_data.txt', 'file')
    delete('para_data.txt');
end
%%
%%please fold this section,,,,it is not elegant
code1 = code_gene('forward_field', 'Ex', 'FOM', 'forward_Ex');
code2 = code_gene('forward_field', 'Ey', 'FOM', 'forward_Ey');
code3 = code_gene('forward_field', 'Ez', 'FOM', 'forward_Ez');
code4 = code_gene('forward_field', 'x', 'FOM', 'forward_x');
code5 = code_gene('forward_field', 'y', 'FOM', 'forward_y');
code_for = [code1, code2, code3, code4, code5];
%--------------------------------------------------------------------------
mesh_accuracy = 3;
mesh_accuracy_meta = 0.02e-6;
%%
for i = 1:Nsteps
    %%%flag forward or backward simulation
    %%%(0):forward;   (1):backward
    
    %forward simulation----------------------------------------------------
    flag = 0;
    tic;
    if(i >= Nsteps / 5 * 4)
        mesh_accuracy_meta = 0.01e-6;
    end
    parameter_name = ['all_parameters','.mat'];
    save(parameter_name, 'flag', 'h_air', 'h_element', 'h_substrate',...
    'index_nano','lambda', 'lambda_ele',  'code_for', 'selpath',...
    'Nsteps', 'param_design', 'space', 'x_array', 'y_array', 'i',...
    'mesh_accuracy', 'mesh_accuracy_meta');
    surface_fdtd = appopen('fdtd');
    code=strcat('cd("', selpath, '");');
    appevalscript(surface_fdtd, code);
    code=strcat('metasurface_gene;');
    appevalscript(surface_fdtd, code);
    appclose(surface_fdtd);
    
    %get forward field-----------------------------------------------------
    filename = 'FOM\forward_Ex.mat';
    load(filename);
    Ex_for = temp;

    filename = 'FOM\forward_Ey.mat';
    load(filename);
    Ey_for = temp;

    filename = 'FOM\forward_Ez.mat';
    load(filename);
    Ez_for = temp;
    
    filename = 'FOM\forward_x.mat';
    load(filename);
    X_for = temp;

    filename = 'FOM\forward_y.mat';
    load(filename);
    Y_for = temp;
    %----------------------------------------------------------------------
    %get backward field----------------------------------------------------
    [X_formesh,Y_formesh]=meshgrid(X_for,Y_for);
    back_Ex = interp2(X, Y, E_x, X_formesh', Y_formesh', 'nearest');
    back_Ey = interp2(X, Y, E_y, X_formesh', Y_formesh', 'nearest');
    back_Ez = zeros(length(Y_for),length(X_for));
    
    FOM = conj(back_Ex) .* Ex_for + conj(back_Ey) .* Ey_for ...
    + conj(back_Ez) .* Ez_for;
    FOM_total = abs(sum(FOM(:)));
    disp(['epoches:',num2str(i)]);
    disp(['FOM: ',num2str(FOM_total)]);
    %backward simulation---------------------------------------------------
    flag = 1;
    save(parameter_name, 'flag', 'h_air', 'h_element','h_substrate',...
    'index_nano', 'lambda', 'lambda_ele','Nsteps', 'param_design',...
    'space', 'x_array', 'y_array', 'i','selpath', 'mesh_accuracy','mesh_accuracy_meta');
    if strcmp(type_beam,'gaussian')
        save(parameter_name,'w0','distance','-append');
    end
    surface_fdtd = appopen('fdtd');
    code=strcat('metasurface_gene;');
    appevalscript(surface_fdtd, code);
    appclose(surface_fdtd);
    
 
    toc;
    gradient = 0;
    gra1 = [];
    gra2 = [];
    for ii = 1:x_array * y_array
        for jj = 1:4
            filename = ['electrc_field_profile_flag_0\','Ez_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            forward = C_x0;
            filename = ['electrc_field_profile_flag_1\','Ez_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            backward = C_x0/balance;
            gradient_ele = forward .* backward;
            gradient = gradient + sum(gradient_ele(:));
        end
        gra1(ii) = 0.5*1i*lambda*(index_nano^2-1.444^2)*gradient;
    end
    
    gradient1 = 0;
    gradient2 = 0;
    for ii = 1:x_array * y_array
        for jj = 1:2:3
            filename = ['electrc_field_profile_flag_0\','D_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            forward = C_x0;
            filename = ['electrc_field_profile_flag_1\','D_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            backward = C_x0/balance;
            gradient_ele = forward .* backward;
            gradient1 = gradient1 + sum(gradient_ele(:));
        end
        for jj = 2:2:4
            filename = ['electrc_field_profile_flag_0\','D_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            forward = C_x0;
            filename = ['electrc_field_profile_flag_1\','D_list_adjoint_data',num2str(ii),'_',num2str(jj),'.mat'];
            load(filename);
            backward = C_x0/balance;
            gradient_ele = forward .* backward;
            gradient2 = gradient2 + sum(gradient_ele(:));
        end
        gra2(ii) = 0.5*1i*lambda*(index_nano^2-1.444^2)*(gradient1 + gradient2);
    end
    F = conj(back_Ex) .* Ex_for + conj(back_Ey) .* Ey_for ...
    + conj(back_Ez) .* Ez_for;
    gradient_origin = 2 * real (conj(sum(F(:))) * (gra1 + gra2));

    [grad_adam, mopt, vopt] = step_adam(gradient_origin, mopt, vopt, i);
    param_design = param_design - step_size * grad_adam * lambda_ele/(i^0.5);
    %enforce constraints
    param_design(param_design <= 100e-9) = 100e-9;
    param_design(param_design >= lambda_ele) = lambda_ele*0.9;

    
    %write iteration data
    fid_fom = fopen('FOM_list.txt', 'a');
    fid_para = fopen('para_data.txt', 'a');
    fprintf(fid_fom, '%s,', num2str(FOM_total));
    fprintf(fid_para, '%s', num2str(param_design));
    fprintf(fid_para, '\n');
    fclose(fid_fom);
    fclose(fid_para);
    %----------------------------------------------------------------------
end
