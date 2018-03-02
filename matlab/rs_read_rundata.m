function [rsR] = rs_read_rundata(folder)

file = strcat(folder,'/output/rs-rundata.nc');

rsR = dlg_read_netcdf(file);

end