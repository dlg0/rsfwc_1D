pro rs_plot_solution, antLoc, dR, nR, $
	e1, e2, e3, $
	kR = kR, r_kR = r_kR, $
	r1 = r1, r2 = r2, r3 = r3

	common dlg_colors
	common plotSwitches
    common switches
   
    ;	Visualise solution

		eRange	= max(sqrt ( abs(e1)^2 + abs(e2)^2 + abs(e3)^2 ))
		eRange	= max(abs(e1))
   		p_r = plot ( r1, e1, layout=[1,3,1],$
				title='Er',ytitle='Er [V/m]',name='Re',window_title='rsfwc_1d',$
                yRange = [-eRange,eRange])
		p_i = plot ( r1, imaginary(e1), color='red',/over,name='Im')

        ;if ar2EField then begin
        ;    p = plot(r1,kj.eR*kj.replace,/over,color='blue',linestyle='--')
        ;    p = plot(r1,imaginary(kj.eR)*kj.replace,/over,color='red',linestyle='--')
        ;    ;p = plot(r1,e1-(kj.er*kj.replace),/over,color='g',thick=2)
        ;    ;p = plot(r1,imaginary(e1-(kj.er*kj.replace)),/over,color='g',thick=1)
        ;endif

		l = legend(target=[p_r,p_i],position=[0.98,0.95],/norm,font_size=10,horizontal_alignment='RIGHT')
    	eRange	= max(abs(e2))
		p_r = plot ( r2, e2, layout=[1,3,2],/current,$
				title='Et',ytitle='Et [V/m]',name='Re',yRange = [-eRange,eRange])
		p_i = plot ( r2, imaginary(e2), color='red',/over,name='Im')
        ;if ar2EField then begin
        ;    p = plot(r1,kj.et*kj.replace,/over,color='blue',linestyle='--')
        ;    p = plot(r1,imaginary(kj.et)*kj.replace,/over,color='red',linestyle='--')
        ;    ;p = plot(r1,e2-(kj.et*kj.replace),/over,color='g',thick=2)
        ;    ;p = plot(r1,imaginary(e2-(kj.et*kj.replace)),/over,color='g',thick=1)
        ;endif


		l = legend(target=[p_r,p_i],position=[0.98,0.62],/norm,font_size=10,horizontal_alignment='RIGHT')
    	eRange	= max(abs(e3))
		p_r = plot ( r3, e3, layout=[1,3,3],/current,$
				title='Ez',ytitle='Ez [V/m]',name='Re', yRange=[-eRange,eRange])
		p_i = plot ( r3, imaginary(e3), color='red',/over,name='Im')
        ;if ar2EField then begin
        ;    p = plot(r1,kj.ez*kj.replace,/over,color='blue',linestyle='--')
        ;    p = plot(r1,imaginary(kj.ez)*kj.replace,/over,color='red',linestyle='--')
        ;    ;p = plot(r1,e3-(kj.ez*kj.replace),/over,color='g',thick=2)
        ;    ;p = plot(r1,imaginary(e3-(kj.ez*kj.replace)),/over,color='g',thick=1)
        ;endif

end
