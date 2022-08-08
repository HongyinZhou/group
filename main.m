clc;
clear;
global x_array
global y_array
global shape_ele
global lambda_ele
global mesh 
global h_element
global space
global Nsteps
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
x_array = 10;
y_array = 1;
%select current folder(the root directory of this program)
selpath = uigetdir;
% parameters set
% element's shape is a cylinder(1) or cuboid(0)
shape_ele = 1;
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
% param_design=importdata('para.txt');
% load('param');
% param_design=radius_mask';
% param_design = rand (1, x_array * y_array) * lambda_ele * 0.8;
param_design(param_design <= 100e-9) = 100e-9;
param_design(param_design >= lambda_ele) = lambda_ele*0.95; 
%--------------------------------------------------------------------------

type_beam = 'plane';
if strcmp(type_beam,'gaussian')
	%define waist radius
	w0 = 0.5e-6;
	%define distance from waist
	distance = -5e-6;
end
%flag for forward(0) and backward simulation(1)
flag = 0;

%desired electric field distribution
phi = 0;%azimuth angle(degrees)
theta = 10;%vertical angle(degrees)
tolerance = 6;%optional

% [list1, list2, list3, list4] = generator_boundary_position(param_design);

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
%import FDTD --------------get forward and backward field------------------
code1 = code_gene('source', 'ExAmp', 'FOM', 'backward_Ex_a');
code2 = code_gene('source', 'ExPhs', 'FOM', 'backward_Ex_p');
code3 = code_gene('source', 'EyAmp', 'FOM', 'backward_Ey_a');
code4 = code_gene('source', 'EyPhs', 'FOM', 'backward_Ey_p');
code5 = code_gene('source', 'EzAmp', 'FOM', 'backward_Ez_a');
code6 = code_gene('source', 'EzPhs', 'FOM', 'backward_Ez_p');
code_back = [code1, code2, code3, code4, code5, code6];
code1 = code_gene('forward_field', 'Ex', 'FOM', 'forward_Ex');
code2 = code_gene('forward_field', 'Ey', 'FOM', 'forward_Ey');
code3 = code_gene('forward_field', 'Ez', 'FOM', 'forward_Ez');
code_for = [code1, code2, code3];
%--------------------------------------------------------------------------
mesh_accuracy = 3;
mesh_accuracy_meta = 0.02e-6;
FOM_pre = 0;
param_pre = [];
gradient_pre = [];
step_size_record = [];
multip_factor = 1;
symbol = 1;
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
    save(parameter_name, 'flag', 'mesh', 'h_air', 'h_element', ...
        'h_substrate','index_nano','lambda', 'lambda_ele',  'code_for',...
        'selpath', 'Nsteps', 'param_design', 'phi', 'theta', 'shape_ele',...
        'shape_ele', 'space', 'tolerance', 'type_beam', 'x_array',...
        'y_array', 'i', 'mesh_accuracy', 'mesh_accuracy_meta');
    if strcmp(type_beam,'gaussian')
        save(parameter_name,'w0','distance','-append');
    end
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
    %----------------------------------------------------------------------
    
    %backward simulation---------------------------------------------------
    flag = 1;
    save(parameter_name, 'flag', 'mesh', 'h_air', 'h_element', ...
        'h_substrate','index_nano', 'lambda', 'lambda_ele', 'code_back',...
    'Nsteps', 'param_design', 'phi', 'theta', 'shape_ele', 'shape_ele'...
    , 'space', 'tolerance', 'type_beam', 'x_array', 'y_array', 'i',...
    'selpath', 'mesh_accuracy','mesh_accuracy_meta');
    if strcmp(type_beam,'gaussian')
        save(parameter_name,'w0','distance','-append');
    end
    surface_fdtd = appopen('fdtd');
    code=strcat('metasurface_gene;');
    appevalscript(surface_fdtd, code);
    appclose(surface_fdtd);
    
    %get backward field and calculate FOM----------------------------------
    filename = 'FOM\backward_Ex_a.mat';
    load(filename);
    ExAmp = temp;
    filename = 'FOM\backward_Ex_p.mat';
    load(filename);
    ExPhs = temp;
    Ex_back = ExAmp .* (cos(ExPhs) + 1i .* sin(ExPhs));

    filename = 'FOM\backward_Ey_a.mat';
    load(filename);
    EyAmp = temp;
    filename = 'FOM\backward_Ey_p.mat';
    load(filename);
    EyPhs = temp;
    Ey_back = EyAmp .* (cos(EyPhs) + 1i .* sin(EyPhs));

    filename = 'FOM\backward_Ez_a.mat';
    load(filename);
    EzAmp = temp;
    filename = 'FOM\backward_Ez_p.mat';
    load(filename);
    EzPhs = temp;
    Ez_back = EzAmp .* (cos(EzPhs) + 1i .* sin(EzPhs));
    
    %normalization
%     operator_back = (conj(Ex_back) .* Ex_back + conj(Ey_back) .* Ey_back ...
%     + conj(Ez_back) .* Ez_back).^0.5;
%     operator_for = (conj(Ex_for) .* Ex_for + conj(Ey_for) .* Ey_for ...
%     + conj(Ez_for) .* Ez_for).^0.5;
%     operator_back = max(operator_back(:));
%     operator_for = max(operator_for(:));
%     Ex_for = Ex_for / operator_for;
%     Ey_for = Ey_for / operator_for;
%     Ez_for = Ez_for / operator_for;
%     Ex_back = Ex_back / operator_back;
%     Ey_back = Ey_back / operator_back;
%     Ez_back = Ez_back / operator_back;
    FOM = conj(Ex_back) .* Ex_for + conj(Ey_back) .* Ey_for ...
    + conj(Ez_back) .* Ez_for;
    FOM_total = sum(FOM(:)) * conj(sum(FOM(:)));
    %get backward field and calculate FOM----------------------------------
%     diff = FOM_total - FOM_pre;
%     if diff < 0
%         if step_size > 0.1
%             symbol = 1;
%             multip_factor = 1;
%         elseif step_size < -0.1
%             symbol = -1;
%             multip_factor = 1;
%         else
%         end
%         step_size = step_size - 1e-4 * (multip_factor ^ 2) * symbol;
%         disp(['step_size: ',num2str(step_size)]);
%         param_design = param_pre + step_size * gradient_pre * lambda_ele;
%         param_design(param_design <= 100e-9) = 100e-9;
%         param_design(param_design >= lambda_ele) = lambda_ele*0.9;
%         disp('adjust step_size ...');
%         multip_factor = multip_factor + 1;
%         continue;
%     else
%         step_size_record = [step_size_record,step_size];
%         multip_factor = 1;
%     end
    disp(['epoches:',num2str(i)]);
    disp(['FOM: ',num2str(FOM_total)]);
%     FOM_pre = FOM_total;
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
%             disp('relat_elec');
%             disp([num2str(max(abs(forward(:)))),'++++++', num2str(max(abs(backward(:))))]);
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
%             disp('relat_13_D');
%             disp([num2str(max(abs(forward(:)))),'++++++', num2str(max(abs(backward(:))))]);
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
%             disp('relat_24_D');
%             disp([num2str(max(abs(forward(:)))),'++++++', num2str(max(abs(backward(:))))]);
            gradient_ele = forward .* backward;
            gradient2 = gradient2 + sum(gradient_ele(:));
        end
        gra2(ii) = 0.5*1i*lambda*(index_nano^2-1.444^2)*(gradient1 + gradient2);
    end
    F = conj(Ex_back) .* Ex_for + conj(Ey_back) .* Ey_for ...
    + conj(Ez_back) .* Ez_for;
    gradient_origin = 2 * real (conj(sum(F(:))) * (gra1 + gra2));
%     grad_max = max(abs(gradient_origin));
%     grad_norm = gradient_origin / grad_max;
    [grad_adam, mopt, vopt] = step_adam(gradient_origin, mopt, vopt, i);
%     GRADIENT = grad_adam .* grad_max;
%     param_design = param_design + step_size * grad_adam * lambda_ele * ...
%         (0.5 * exp(-i/20) + (1 / i ^ 2) * 0.5);%attenuation factor
%     param_pre = param_design;
%     gradient_pre = grad_adam;
    param_design = param_design - step_size * grad_adam * lambda_ele/(i ^ 0.5);
    %enforce constraints
    param_design(param_design <= 100e-9) = 100e-9;
    param_design(param_design >= lambda_ele) = lambda_ele*0.9;
%     para_data = [para_data;param_design];
    
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
filename='stepsize_list.txt';
save(filename,'step_size_record','-ascii');