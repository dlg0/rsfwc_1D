function rs_read_rundata, runFolderName

    runDataFile = expand_path(RunFolderName)+'/output/rs-rundata.nc'

    return, dlg_read_netcdf(runDataFile)

end
