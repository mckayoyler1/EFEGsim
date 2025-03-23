dirlist  = dir('/Users/mckayoyler/Documents/Code/Samu/fieldtrip-master/template/electrode/*.*');
filename = {dirlist(~[dirlist.isdir]).name}';

for i=1:length(filename)
   clear pos; clear path; clear name; clear ext;
   elec = ft_read_sens(filename{i});

   [path,name,ext] = fileparts(filename{i});
   pos = elec.elecpos;
   writematrix(pos,sprintf('%s.csv',name));
end

