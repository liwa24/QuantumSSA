function out = writeIC(directory)

addpath(directory);
files = dir(directory);
files(ismember({files.name},{'.','..'})) = [];
files(endsWith({files.name},{'.mat'})) = [];
out = [];

classes_names = ["BN"; "BS"; "DN"; "DPO"; "DRO"; "DS"; "L1_A"; "L1_HN"; "L1_HS"; ...
"L1_L"; "L1_V"; "L2_A"; "L2_HN"; "L2_HS"; "L2_L"; "L2_V"; "L3_A"; "L3_HN"; "L3_HS";...
"L3_L"; "L3_V"; "L4_A"; "L4_LP"; "L4_SP"; "L4_V"; "L5_A"; "L5_LP"; "L5_SP"; "L5_V";...
"LPOE"; "LPOW"; "R11"; "R12"; "R13"; "R14"; "R21"; "R23"; "R31"; "R32"; "R34"; ...
"R41"; "R43"];

keys = (1:length(classes_names))';
dict_IC = dictionary(classes_names,keys);

for i = 1:1:length(files)
  filename = files(i).name;
  fd = fopen(filename,'r');
  start_index = find(filename == '_',1) +1 ;
  end_index = find(filename == '.') -1; 
  entry = filename(start_index:end_index);
  key = dict_IC(entry)
  assert(fd > 0,['Can''t open file ' files(i).name ' for reading.']);

  iter = readmatrix(files(i).name,'NumHeaderLines',1);
  out = [out; [ones(size(iter,1),1)*key,iter(:,[2:9, 11])]];
end
fprintf("processed %d files.\n", length(files));