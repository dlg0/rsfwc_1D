function rs_read_solution, runFolderName

    runDataFile = expand_path(RunFolderName)+'/output/rs-solution.nc'

    return, dlg_read_netcdf(runDataFile)

end
