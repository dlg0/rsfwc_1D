function rs_read_rundata, runFolderName

    runDataFile = expand_path(RunFolderName)+'/output/rs-rundata.nc'

    ncdf_list, runDataFile, /quiet, /variables, vName=vars

    ncdf_get, runDataFile, vars, data

    h = hash()

    N = n_elements(vars)

    for i=0,N-1 do begin

        h0 = hash(vars[i],data[vars[i],'value'])
        h = h + h0

    endfor 

    return, h

end
