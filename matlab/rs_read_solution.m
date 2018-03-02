function [rs] = rs_read_solution(folder)

file = strcat(folder,'/output/rs-solution.nc');

rs = dlg_read_netcdf(file);

end