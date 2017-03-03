	

	; inverse method for random numbers

	
	omega_m = 0.288
	omega_l = 1d0 - omega_m
	H0      = 69.33
	A_snIa  = 0.0013 ; #/Msun;	Maoz +(12)


	t_ini = 50d6
	t_fin = galage(0,1000,omega_m=omega_m,lambda=omega_l,H0=H0)
	
	; DTD = A_DTD* t^-1
	; 1 = int A_DTD / t dt
	;   = A_DTD*(log_e (t_f) - log_e (t_i))
	A_DTD = 1d0/(alog(t_fin) - alog(t_ini))
	
	; n_sn(<t) = A_DTD*(log (t) - log(t_ini)) ; 0 < n_sn < 1
	; log(t) = n_sn/A_DTD + log(t_ini)
	; t = exp(n_sn/A_DTD + log(t_ini))

	nsample = 10000
	iseed   = -3997
	xx = randomu(iseed,nsample)
	tt = exp(xx / A_DTD + alog(t_ini))

	plothist,tt/1d6,bin=10,/xlog	
end
