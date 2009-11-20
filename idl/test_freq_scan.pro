;   scan a range of frequencies for a cavity to 
;   look for resonances

pro test_freq_scan

    w = 30d6 * 2d0 * !pi
    nR = 128 
    
    nFreq   = 301
    eR_all  = complexArr ( nR, nFreq )
    
    for i = 0, nFreq-1 do begin
    
        print, 'Freq: ', w + w * i / 500.0, i
        rsfwc_1d, eR = eR_tmp, $
                    w = w + w * i/500.0, $
                    nR = nR, $
                    /noPlot
            eR_all[*,i] = eR_tmp / (i+1)
    
    endfor
    
    loadct, 3, /sil
    device, decomposed = 0
    !p.background = 255
    levels = fIndGen ( 21 ) * max ( abs ( eR_all ) ) / 100.0
    colors = 254 - bytScl ( levels, top = 253 ) + 1 
    contour, abs ( eR_all ), $
        c_colors = colors, $
        levels = levels, $
        /fill, $
        color = 0
stop
end
