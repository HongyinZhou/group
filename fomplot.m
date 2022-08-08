clc;
clear;
fom = importdata('FOM_list.txt');
x = [1:length(fom)];
plot(x,fom);