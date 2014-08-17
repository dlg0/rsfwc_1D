pro rsfwc_read_iterations, nIterations, RightName, MPE_FileName = MPE_FileName, $
		vorpal=vorpal, freq = freq

    ;RightName = '/Users/dg6/scratch/rsfwc_1d/ar2_vorpal/right_simple_'
    ;CompareDir = '~/scratch/rsfwc_1d/ar2_vorpal/simple_full'
    CompareDir = '~/scratch/aorsa2d/ar2_vorpal/simple_full'

    for n=0,nIterations-1 do begin

        right_ThisDir =  expand_path(RightName)+StrCompress(string(n),/rem)
        print, 'Reading ... ', right_ThisDir
		if keyword_set(vorpal) then begin
    		right_s = ar2_read_vorpal (runFolderName = right_ThisDir, /oneD, freq = freq )
		endif else begin
        	right_s = rsfwc_read_solution (right_ThisDir)
		endelse

        if data_jP_r eq !null then begin
            nR = n_elements(right_s.jP_r)

            data_E_r = ComplexArr(nR,nIterations)
            data_E_t = ComplexArr(nR,nIterations)
            data_E_z = ComplexArr(nR,nIterations)

            data_jP_r = ComplexArr(nR,nIterations)
            data_jP_t = ComplexArr(nR,nIterations)
            data_jP_z = ComplexArr(nR,nIterations)
 
        endif

        data_jP_r[*,n] = right_s.jP_r
        data_jP_t[*,n] = right_s.jP_t
        data_jP_z[*,n] = right_s.jP_z

        data_E_r[*,n] = right_s.E_r
        data_E_t[*,n] = right_s.E_t
        data_E_z[*,n] = right_s.E_z

    endfor

    ;compare_s = rsfwc_read_solution (CompareDir)
    compare_s = ar2_read_solution (CompareDir, 1)

    ;compare_r = total(compare_s.jP_r,3)
    ;compare_t = total(compare_s.jP_t,3)
    ;compare_z = total(compare_s.jP_z,3)

    compare_r = compare_s.E_r
    compare_t = compare_s.E_t
    compare_z = compare_s.E_z

    mpe_jP_r = kj_mpe(data_jP_r)
    mpe_jP_t = kj_mpe(data_jP_t)
    mpe_jP_z = kj_mpe(data_jP_z)

    mpe_E_r = kj_mpe(data_E_r)
    mpe_E_t = kj_mpe(data_E_t)
    mpe_E_z = kj_mpe(data_E_z)

    ;mpe_r = mpe_jP_r
    ;mpe_t = mpe_jP_t
    ;mpe_z = mpe_jP_z

    mpe_r = mpe_E_r
    mpe_t = mpe_E_t
    mpe_z = mpe_E_z

    ;data_r = data_jP_r
    ;data_t = data_jP_t
    ;data_z = data_jP_z

    data_r = data_E_r
    data_t = data_E_t
    data_z = data_E_z

    range = max(abs([mpe_r,mpe_t,mpe_z]))

    p=plot(right_s.r, data_r[*,0],layout=[2,3,1], window_title=MPE_FileName)
    for n=1,nIterations-1 do p=plot(right_s.r,data_r[*,n],/over)
    p=plot(compare_s.r, compare_r, color='g',thick=3,/over)
    p=plot(right_s.r, mpe_r,/over,color='b',thick=3,yRange=[-range,range])

    p=plot(right_s.r, data_t[*,0],layout=[2,3,3],/current)
    for n=1,nIterations-1 do p=plot(right_s.r,data_t[*,n],/over)
    p=plot(compare_s.r, compare_t, color='g',thick=2,/over)
    p=plot(right_s.r, mpe_t,/over,color='b',thick=2,yRange=[-range,range])

    p=plot(right_s.r, data_z[*,0],layout=[2,3,5],/current)
    for n=1,nIterations-1 do p=plot(right_s.r,data_z[*,n],/over)
    p=plot(compare_s.r, compare_z, color='g',thick=2,/over)
    p=plot(right_s.r, mpe_z,/over,color='b',thick=2,yRange=[-range,range])

    p=plot(right_s.r, imaginary(data_r[*,0]),layout=[2,3,2],/current,color='r')
    for n=1,nIterations-1 do p=plot(right_s.r,imaginary(data_r[*,n]),/over,color='r')
    p=plot(compare_s.r, imaginary(compare_r), color='g',thick=2,/over)
    p=plot(right_s.r, imaginary(mpe_r),/over,color='b',thick=2,yRange=[-range,range])

    p=plot(right_s.r, imaginary(data_t[*,0]),layout=[2,3,4],/current,color='r')
    for n=1,nIterations-1 do p=plot(right_s.r,imaginary(data_t[*,n]),/over,color='r')
    p=plot(compare_s.r, imaginary(compare_t), color='g',thick=2,/over)
    p=plot(right_s.r, imaginary(mpe_t),/over,color='b',thick=2,yRange=[-range,range])

    p=plot(right_s.r, imaginary(data_z[*,0]),layout=[2,3,6],/current,color='r')
    for n=1,nIterations-1 do p=plot(right_s.r,imaginary(data_z[*,n]),/over,color='r')
    p=plot(compare_s.r, imaginary(compare_z), color='g',thick=2,/over)
    p=plot(right_s.r, imaginary(mpe_z),/over,color='b',thick=2,yRange=[-range,range])


    fName = 'restart_mpe.nc'
    if keyword_set(MPE_FileName) then fName = MPE_FileName
	nc_id = nCdf_create (fName, /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(right_s.R) )
	if not keyword_set(vorpal) then $
			nrH_id = nCdf_dimDef ( nc_id, 'nR_', n_elements(right_s.R_) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	if not keyword_set(vorpal) then $
		rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', nr_id, /float )
	if not keyword_set(vorpal) then $
		zH_id = nCdf_varDef ( nc_id, 'z_', nrH_id, /float )

	mpe_E_r_re_id = nCdf_varDef ( nc_id, 'E_r_re', nr_id, /float )
	mpe_E_r_im_id = nCdf_varDef ( nc_id, 'E_r_im', nr_id, /float )
	mpe_E_t_re_id = nCdf_varDef ( nc_id, 'E_p_re', nr_id, /float )
	mpe_E_t_im_id = nCdf_varDef ( nc_id, 'E_p_im', nr_id, /float )
	mpe_E_z_re_id = nCdf_varDef ( nc_id, 'E_z_re', nr_id, /float )
	mpe_E_z_im_id = nCdf_varDef ( nc_id, 'E_z_im', nr_id, /float )

	mpe_jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	mpe_jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	mpe_jP_t_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	mpe_jP_t_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	mpe_jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	mpe_jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, r_id, right_s.r
	if not keyword_set(vorpal) then $
		nCdf_varPut, nc_id, rH_id, right_s.r_

	nCdf_varPut, nc_id, z_id, right_s.z
	if not keyword_set(vorpal) then $
		nCdf_varPut, nc_id, zH_id, right_s.z_

	nCdf_varPut, nc_id, mpe_E_r_re_id, real_part(mpe_E_r)
	nCdf_varPut, nc_id, mpe_E_r_im_id, imaginary(mpe_E_r) 
	nCdf_varPut, nc_id, mpe_E_t_re_id, real_part(mpe_E_t) 
	nCdf_varPut, nc_id, mpe_E_t_im_id, imaginary(mpe_E_t) 
	nCdf_varPut, nc_id, mpe_E_z_re_id, real_part(mpe_E_z) 
	nCdf_varPut, nc_id, mpe_E_z_im_id, imaginary(mpe_E_z) 

	nCdf_varPut, nc_id, mpe_jP_r_re_id, real_part(mpe_jP_r)
	nCdf_varPut, nc_id, mpe_jP_r_im_id, imaginary(mpe_jP_r) 
	nCdf_varPut, nc_id, mpe_jP_t_re_id, real_part(mpe_jP_t) 
	nCdf_varPut, nc_id, mpe_jP_t_im_id, imaginary(mpe_jP_t) 
	nCdf_varPut, nc_id, mpe_jP_z_re_id, real_part(mpe_jP_z) 
	nCdf_varPut, nc_id, mpe_jP_z_im_id, imaginary(mpe_jP_z) 

	nCdf_close, nc_id

end
