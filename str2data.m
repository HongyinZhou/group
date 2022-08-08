function out = str2data(in)
startn = regexpi(in, '[');
endn = regexpi(in, ']');
unit = in((startn + 1):(endn - 1));
data = str2num(in(1:startn-1));

if strcmp(unit,'m')
    out = data;
elseif strcmp(unit,'mm')
    out = data * 10^-3;
elseif strcmp(unit,'um')
    out = data * 10^-6;
elseif strcmp(unit,'nm')
    out = data * 10^-9;
else
    disp('please check the unit of input');
end