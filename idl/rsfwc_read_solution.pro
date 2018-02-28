function rsfwc_read_solution, runFolderName

    SolutionFile = expand_path(runFolderName)+'/output/rs-solution.nc'

    rs = dlg_read_netcdf(SolutionFile)

    return, rs 

end
