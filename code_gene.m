function code = code_gene(ele_name, data_name, fold_name, file_name)
    code1 = ['temp=getdata("', ele_name, '", "', data_name, '");'];
    code2 = ['filename = base_path + "\', fold_name, '"+ "\', file_name, '";'];
    code3 = 'matlabsave(filename, temp);';
    code =  [code1, code2, code3];
end