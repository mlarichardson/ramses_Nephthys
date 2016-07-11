	; This routine requires the Starburst99 outputs generated with galaxy_extra_field.f
	; for five different metallicities. 



	; starburst99 output prefix (ex. m41, m42, m43..)
	;------------------------------------------
	output_prefix  = 'm4'
	output_suffix  = '_bh50'
	;------------------------------------------
	szgrid = ['0.0004','0.004','0.008','0.02','0.05']
	zgrid  = double(szgrid)
	nz     = n_elements(zgrid)
	nt_tar = 250L ; number of time bin we would to print out


	; output file
	output_file = 'swind_krp_pagb'+output_suffix+'.dat'
	;------------------------------------------
	openw,11,output_file,/f77_unform
	;------------------------------------------


	; If the mass of the indivual species is correct, we will be using those values
	; However, the summation of H-Fe yield from SNII is not the same as the total mass field from SNII.
	; Make sure which one to use.
	individual_species_correct = 0 



	set_plot,'ps'
	!p.font=0
	!p.charsize=1.3
	!p.thick=3	
	device,file='SW.ps',xs=14,ys=14,/color
	xr = [1d5,1.4d10]
	yr = [0,0.40]
	pos = [0.16,0.12,0.96,0.96]
	xtit = 'age [yr]'
	ytit = 'total ejecta (<t)'
	
	plot, xr,yr, xr=xr, yr=yr, /xlog, xtit=xtit, ytit=ytit, /nodata, pos=pos



	spec = ['H','He','C','N','O','Mg','Si','S','Fe']
	nspec = n_elements(spec)

	loadct,33
	icol = [20,80,120,200,250]

	for iz=0,4 do begin
		siz = strtrim(iz+1,2)
		filename = output_prefix+siz+output_suffix+'/standard.yield1'
		TKNOTE, filename
		readcol,filename,skip=7,time,H1,He1,C1,N1,O1,Mg1,Si1,S1,Fe1,$
		                             H2,He2,C2,N2,O2,Mg2,Si2,S2,Fe2,Msw,Msn,Mall,/silent

;	sanity check
		sanity_check = 0 
		if iz eq 3 and sanity_check then begin
			device,/close
			set_plot,'x'
			H1 = 10d0^H1
			He1 = 10d0^He1
			C1  = 10d0^C1
			N1  = 10d0^N1
			O1  = 10d0^O1
			Mg1  = 10d0^Mg1
			Fe1  = 10d0^Fe1
			Si1  = 10d0^Si1
			S1  = 10d0^S1
			Msw = 10d0^Msw
			Zsw = C1+ N1 + O1+Mg1+Fe1+Si1+S1
			
			;plot, time, H1/Msw, /xlog
			plot, time, Zsw, /xlog
			oplot, time, O1, linestyle=2
			oplot, time, Mg1, linestyle=3
			oplot, time, C1, linestyle=1
			stop
		endif

		;----------------------------------------
		; Starburst99 is normalised to 10^6 Msun
		;----------------------------------------
		Msw -= 6d0
		Msn -= 6d0
		Mall -= 6d0
		for ispec = 0, nspec-1 do begin
			atom   = spec[ispec]
			; for stellar winds
			command = atom+'1 -= 6d0' ; H1 -= 6d0
			dum    = execute(command)
			; for SNII
			command = atom+'2 -= 6d0' ; H2 -= 6d0
			dum    = execute(command)
			
		endfor


		;----------------------------------------
		; make log scale to non-log
		;----------------------------------------
		Msw = 10d0^Msw
		Msn = 10d0^Msn
		Mall = 10d0^Mall

		for ispec = 0, nspec-1 do begin
			atom   = spec[ispec]
			; for stellar winds
			command = atom+'1 = 10d0^'+atom+'1' ; H1 =  10d0^H1
			dum    = execute(command)
			; for SNII
			command = atom+'2 = 10d0^'+atom+'2' ; H2 =  10d0^H2
			dum    = execute(command)
		endfor



		;----------------------------------------
		; create array for cumulative numbers
		;----------------------------------------
		nn = n_elements(time)
		cMsw = dblarr(nn)
		cMsn = dblarr(nn)
		cMall = dblarr(nn)

		cX1  = dblarr(nn)	
		cY1  = dblarr(nn)
		cZ1  = dblarr(nn)

		cX2  = dblarr(nn)	
		cY2  = dblarr(nn)
		cZ2  = dblarr(nn)

		for ispec = 0, nspec-1 do begin
			atom   = spec[ispec]
			; for stellar winds
			command = 'c'+atom+'1 = dblarr(nn)' ; cH1 = dblarr(nn)
			dum    = execute(command)
			; for SNII
			command = 'c'+atom+'2 = dblarr(nn)' ; cH2 = dblarr(nn)
			dum    = execute(command)
		endfor


		;----------------------------------------
		; collect data
		;----------------------------------------
		if iz eq 0 then begin
			; for stellar winds
			cMsw_data = dblarr(nn,nz)
			cEsw_data = dblarr(nn,nz)
			cMZsw_data = dblarr(nn,nz)
			for ispec = 0, nspec-1 do begin
				atom   = spec[ispec]
				; stellar winds
				command = 'cM'+atom+'sw_data = dblarr(nn,nz)' ; cMHsw_data = dblarr(nn,nz)
				dum    = execute(command)
				; SNII
				command = 'cM'+atom+'sn_data = dblarr(nn,nz)' ; cMHsn_data = dblarr(nn,nz)
				dum    = execute(command)
			endfor

			; for SNII
			cMsn_data = dblarr(nn,nz)
			cMZsn_data = dblarr(nn,nz)
		endif

		;----------------------------------------
		; time integration
		;----------------------------------------
		tprev = 0d0
		for i=0L,nn-1L do begin
			dt    = time[i]-tprev
			j     = i-1
			if i eq 0 then j=i

			cMsw[i] = cMsw[j] + Msw[i]*dt
			cMsn[i] = cMsn[j] + Msn[i]*dt
			cMall[i] = cMall[j] + Mall[i]*dt

			for ispec = 0, nspec-1 do begin
				atom   = spec[ispec]

				; stellar winds
				command = 'c'+atom+'1[i] = c'+atom+'1[j] + '+atom+'1[i]*dt' ; cH1 [i] = cH1 [j] + H1 [i]*dt
				dum    = execute(command)
	
				; SNII
				command = 'c'+atom+'2[i] = c'+atom+'2[j] + '+atom+'2[i]*dt' ; cH2 [i] = cH2 [j] + H2 [i]*dt
				dum    = execute(command)
			
			endfor

			cX1[i] = cH1[i]
			cY1[i] = cHe1[i]
			cZ1[i] = cC1[i] + cN1[i] + cO1[i] + cMg1[i] + cSi1[i] + cS1[i] + cFe1[i]	

			cX2[i] = cH2[i]
			cY2[i] = cHe2[i]
			cZ2[i] = cC2[i] + cN2[i] + cO2[i] + cMg2[i] + cSi2[i] + cS2[i] + cFe2[i]

				
			tprev = time[i]
		endfor

		oplot, time, cMsw , col=icol[iz]
		oplot, time, cMsn , col=icol[iz], linestyle=1
		oplot, time, cMall , col=icol[iz], linestyle=2

		if individual_species_correct then begin
			oplot, time, (cX1+cY1+cZ1), linestyle=0
			oplot, time, (cX2+cY2+cZ2), linestyle=2
		endif

		cMZsw = cZ1

		if individual_species_correct then begin
			cMsn  = cX2 + cY2 + cZ2
			cMZsn = cZ2
			correction_factor = 1
		endif else begin
			n0 = n_elements(cMsn)
			correction_factor = cMsn[n0-1]/(cX2[n0-1] + cY2[n0-1] + cZ2[n0-1])
			cMZsn = cZ2*correction_factor

			fd = where(cMZsn gt cMZsn[n_elements(cMZsn)-1], nfd)
			if nfd gt 0 then cMZsn[fd] = cMZsn[n_elements(cMZsn)-1]
		endelse


		;----------------------------------------
		; collect the cumulative data
		;----------------------------------------
		cMsw_data[*,iz] = cMsw
		cMZsw_data[*,iz] = cMZsw
		cMsn_data[*,iz] = cMsn
		cMZsn_data[*,iz] = cMZsn
		for ispec = 0, nspec-1 do begin
			atom   = spec[ispec]

			; stellar winds
			command = 'cM'+atom+'sw_data[*,iz] = c'+atom+'1' ; cMHsw_data[*,iz] = cH1
			dum    = execute(command)
	
			; SN II
			command = 'cM'+atom+'sn_data[*,iz] = c'+atom+'2 * correction_factor' ; cMHsn_data[*,iz] = cH2 * correction_factor
			dum    = execute(command)

		endfor

		filename = output_prefix+siz+output_suffix+'/standard.power1'
		readcol,filename,skip=7,dum,dum,dum,dum,dum,dum,cEsw,/silent
		cEsw = 10d0^cEsw/1d6

		cEsw_data[*,iz] = cEsw

		;print, max(cMsn)/1d6


	endfor

	item = 'Z='+szgrid
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2

	item2 = ['Winds','SNe','Winds+SNe']
	legend, item2, pos=[8d4,0.25],spacing=1.6,box=0,linestyle=[0,1,2],/short,charsize=1.0
	device,/close

	set_plot,'x'


	;----------------------------------------
	; make logarithmic scale 
	;----------------------------------------
	cMsw_data = alog10(cMsw_data)
	cEsw_data  = alog10(cEsw_data)
	cMZsw_data = alog10(cMZsw_data)
	cMsn_data = alog10(cMsn_data)
	cMZsn_data = alog10(cMZsn_data)
	for ispec = 0, nspec-1 do begin
		atom   = spec[ispec]

		; stellar winds
		command = 'cM'+atom+'sw_data = alog10(cM'+atom+'sw_data)' ; cMHsw_data = alog10(cMHsw_data)
		dum    = execute(command)
	
		; SNII
		command = 'cM'+atom+'sn_data = alog10(cM'+atom+'sn_data)' ; cMHsn_data = alog10(cMHsn_data) 
		dum    = execute(command)

	endfor

	time = alog10(time)
	tmin = min(time)
	tmax = max(time)
	tbin = (tmax-tmin)/double(nt_tar)
	logt = dindgen(nt_tar)*tbin + tmin

	;----------------------------------------
	; prepare data for output
	;----------------------------------------
	; only for SW
	log_cm_sw   = dblarr(nt_tar, nz)
	log_ce_sw   = dblarr(nt_tar, nz)
	log_cmz_sw   = dblarr(nt_tar, nz)
	for ispec = 0, nspec-1 do begin
		atom   = spec[ispec]

		; stellar winds
		command = 'log_cm'+atom+'_sw = dblarr(nt_tar,nz)' ; log_cmH_sw = dblarr(nt_tar,nz)
		dum    = execute(command)
	
		; SNII
		command = 'log_cm'+atom+'_sn = dblarr(nt_tar,nz)' ; log_cmH_sn = dblarr(nt_tar,nz)
		dum    = execute(command)

	endfor

	;----------------------------------------
	; linear interpolation
	;----------------------------------------
	for iz=0,nz-1 do begin
		log_cm_sw[*,iz] = interpol(cMsw_data[*,iz],time,logt)
		log_ce_sw[*,iz] = interpol(cEsw_data[*,iz],time,logt)
		log_cmz_sw[*,iz] = interpol(cMZsw_data[*,iz],time,logt)

		for ispec = 0, nspec-1 do begin
			atom   = spec[ispec]

			; stellar winds
			; log_cmH_sw[*,iz] = interpol(cMHsw_data[*,iz],time,logt)
			command  = 'log_cm'+atom+'_sw[*,iz] = '
			command += 'interpol(cM'+atom+'sw_data[*,iz],time,logt)' 
			dum    = execute(command)
	
			; SN2
			; log_cmH_sn[*,iz] = interpol(cMHsn_data[*,iz],time,logt)
			command  = 'log_cm'+atom+'_sn[*,iz] = '
			command += 'interpol(cM'+atom+'sn_data[*,iz],time,logt)' 
			dum    = execute(command)
		endfor
	endfor

	ngrid_z = n_elements(zgrid)
	ngrid_t = n_elements(logt)
	writeu,11, long(ngrid_t),long(ngrid_z)
	writeu,11, double(logt)
	writeu,11, double(alog10(zgrid))
	for iz=0,ngrid_z-1 do writeu,11, double(log_cm_sw[*,iz])
	for iz=0,ngrid_z-1 do writeu,11, double(log_ce_sw[*,iz])
	for iz=0,ngrid_z-1 do writeu,11, double(log_cmz_sw[*,iz])

	for ispec = 0, nspec-1 do begin
		atom   = spec[ispec]

		; stellar winds
		; for iz=0,ngrid_z-1 do writeu,11, double(log_cmH_sw[*,iz])
		command  = 'for iz=0,ngrid_z-1 do writeu, 11, '
		command += ' double(log_cm'+atom+'_sw[*,iz])'
		dum    = execute(command)
	endfor
	
	close, 11


	;----------------------------------------
	; sanity check
	;---------------------------------------
	nt = 0L
	nz = 0L
	openr,11,output_file,/f77_unform
	readu,11, nt, nz
	log_tgrid = dblarr(nt)
	readu,11,log_tgrid
	log_zgrid = dblarr(nz)
	readu,11,log_zgrid
	dum1 = dblarr(nt)
	log_cm = dblarr(nt,nz)
	for iz=0,nz-1 do begin
		readu,11, dum1
		log_cm[*,iz] = dum1	
	endfor
	log_ce = dblarr(nt,nz)
	for iz=0,nz-1 do begin
		readu,11, dum1
		log_ce[*,iz] = dum1	
	endfor
	log_cmz = dblarr(nt,nz)
	for iz=0,nz-1 do begin
		readu,11, dum1
		log_cmz[*,iz] = dum1	
	endfor

	for ispec = 0, nspec-1 do begin
		atom   = spec[ispec]

		comm = 'log_cm'+atom+' = dblarr(nt,nz)'
		dum = execute(comm)

		for iz=0,nz-1 do begin
			readu,11,dum1
			comm = 'log_cm'+atom+'[*,iz] = dum1'
			dum = execute(comm)	
		endfor
	endfor

	close,11
	
	set_plot,'ps'
	device,file='SW_check.ps',xs=14,ys=14,/color
	pos = [0.16,0.12,0.96,0.96]

	xtit = 'age [yr]'
	xr = [1d5,1.4d10]
	;---------------
	; mass loss
	;---------------
	yr = [0,0.40]
	ytit = 'total ejecta (<t)'
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog
	loadct,33
	for iz=0,nz-1 do begin
		oplot, 10d0^log_tgrid, 10d0^log_cm[*,iz],col=icol[iz]	
	endfor
	item = 'Z='+szgrid
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2

	;---------------
	; energy
	;---------------
	yr = [1d46,1d49]
	ytit = 'total energy (<t)'
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog,/ylog
	loadct,33
	for iz=0,nz-1 do begin
		oplot, 10d0^log_tgrid, 10d0^log_ce[*,iz],col=icol[iz]	
	endfor
	item = 'Z='+szgrid
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2

	;---------------
	; metallicity
	;---------------
	yr = [1d-4,2d-1]
	ytit = 'mean metallicity (<t)'
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog,/ylog
	loadct,33
	for iz=0,nz-1 do begin
		yy = 10d0^(log_cmz[*,iz]-log_cm[*,iz])
		oplot, 10d0^log_tgrid, yy,col=icol[iz]	
	endfor
	item = 'Z='+szgrid
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2

	;---------------
	; fraction
	;---------------
	yr = [1d-4,2d-1]

	for ispec=0,nspec-1 do begin

		comm = 'log_cmE = log_cm'+spec[ispec]+'[*,*]'
		dum  = execute(comm)

		ytit = 'mean fraction [' + spec[ispec]+']'
		loadct,0
		yr   = minmax(log_cmE-log_cm)
		yr   = [yr[0]-0.1,yr[1]+0.1]
		plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog
		loadct,33
		for iz=0,nz-1 do begin
			
			yy = log_cmE[*,iz]-log_cm[*,iz]
			oplot, 10d0^log_tgrid, yy,col=icol[iz]	
		endfor
		item = 'Z='+szgrid
		legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2

	endfor

	device,/close		


	output_file = 'sn2_krp_pagb'+output_suffix+'.dat'
	set_plot,'ps'
	device,file='SN_check.ps',xs=14,ys=14,/color
	pos = [0.16,0.12,0.96,0.96]

	xtit = 'age [yr]'
	xr = [1d5,1.4d10]

	yr = [0,0.40]
	ytit = 'total ejecta (<t)'
	;---------------
	; mass loss
	;---------------
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog
	loadct,33
	cm_sn_tot = dblarr(nz)
	nn = n_elements(cMsn_data[*,0])
	for iz=0,nz-1 do begin
		oplot, 10d0^time,10d0^cMsn_data[*,iz],col=icol[iz]
		cm_sn_tot[iz] = 10d0^cMsn_data[nn-1,iz]
	endfor
	item = 'Z='+szgrid
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2


	text ='  log_SNII_m = (/'
	for iz=0,nz-1 do begin
		text+=  strtrim(alog10(cm_sn_tot[iz]),2)
		if iz ne nz-1 then text+=','
	endfor
	text+='/)'
	print, text
	


	openw,10,output_file,/f77_unformatted
	writeu,10,long(nz)
	writeu,10,long(nspec)
	print, nz, nspec
	writeu,10,double(alog10(cm_sn_tot))
	print, alog10(cm_sn_tot)
	;---------------
	; mass loss fit
	;---------------
	loadct,0
	plot, zgrid, cm_sn_tot, /xlog,/ylog,psym=6,yr=[0.05,0.3],$
		xtit='Z',ytit='mass loss fraction due to SNe'
	afit = linfit(alog10(zgrid),alog10(cm_sn_tot))
	xx = findgen(11)/10.*(alog10(0.1)-alog10(0.0001))+alog10(0.0001)
	yy = afit[0]+afit[1]*xx
	oplot, 10d0^xx, 10d0^yy
	text = 'log10(m_sn)='+string(afit[0],format='(f7.4)')+'+ '+$
			string(afit[1],format='(f7.4)')+' log10(Z)'
	legend, text, /bottom, box=0, charsize=1.2

	;---------------
	; metallicity
	;---------------
	yr = [1d-2,1]
	ytit = 'mean metallicity (<t)'
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog,/ylog
	loadct,33
	item = 'Z='+szgrid
	text = '  log_SNII_Z=(/'
	cmz_sn_tot = dblarr(nz)
	for iz=0,nz-1 do begin
		yy = 10d0^(cMZsn_data[*,iz]-cMsn_data[*,iz])
		oplot, 10d0^time,yy,col=icol[iz]

		ej_Z  = yy[n_elements(yy)-1]
		cmz_sn_tot[iz] = ej_Z
	
		yield = ej_Z-zgrid[iz]
		item[iz] += ' (ramses yield='+string(yield,format='(f4.2)')+')'
		text += strtrim(alog10(ej_Z),2)
		if iz ne nz-1 then text+=','
	endfor
	text += '/)'
	print, text
	legend, item, spacing=1.6, color=icol, psym=6,box=0,charsize=1.2,/bottom

	writeu,10,double(alog10(cmz_sn_tot))
	text ='  log_Zgrid = (/'
	for iz=0,nz-1 do begin
		text += strtrim(alog10(zgrid[iz]),2)
		if iz ne nz-1 then text += ','
	endfor
	text += '/)'
	print,text

	;---------------
	; element
	;---------------
	xr = [1d-4,0.1]
	yr = [1d-5,1]
	ytit = 'fraction'
	xtit = 'Z'
	loadct,0
	plot, xr, yr, xr=xr,yr=yr,/xs,/ys,xtit=xtit,ytit=ytit,pos=pos,/nodata,/xlog,/ylog
	loadct,33
	nn = n_elements(cMHsn_data[*,0])
	icol2 = (findgen(nspec)+0.5)/nspec*250 
	for ispec=0,nspec-1 do begin
		comm = 'yy = cM'+spec[ispec]+'sn_data[nn-1,*]-cMsn_data[nn-1,*]'
		dum  = execute(comm)
		yy   = 10d0^(reform(yy))
		oplot, zgrid,yy,col=icol2[ispec], psym=-6

		comm = '  log_SNII_'+spec[ispec]+'=(/'
		for iz=0,nz-1 do begin
			comm +=strtrim(alog10(yy[iz]),2)
			if iz lt nz-1 then comm+=', '
		endfor
		comm += '/)'
		print, comm

		writeu,10,alog10(yy)
	endfor
	legend, spec, spacing=1.6, color=icol2, psym=6,box=0,charsize=1.2,/bottom
	
	close,10	
	device,/close
	
	set_plot,'x'
end	
