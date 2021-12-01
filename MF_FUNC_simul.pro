pro MF

;device,decomposed=0,retain=2
;loadct,39

M        = 3
C        =0    ;the number of pixels I use to calculate
nside    = 512.
npix     = 12.*nside*nside
delta_nu = 0.25
nmaps    = 6
;miss    = -1.63750e30
w        = dblarr(npix)
a1       = dblarr(npix)
v        = dblarr(7,25)
;sign    =intarr(npix)

mask='planck_spin0_Gauss.fits'
read_fits_map,  mask, a1
ud_grade, a1, w, nside_out=nside, order_in='ring', order_out='ring'
;print,w
;stop


for j=0.,npix-1   do begin
	if w[j] ne 0.  then begin
		;sign[j]=1
		C=C+1
	endif
endfor

;print,C

for n=1001,1000+M    do begin

    pixels    = lindgen(npix)
    theta     = dblarr(npix)
    theta_cut = dblarr(C)
    phi       = dblarr(npix)
    phi_cut   = dblarr(C)
    Mink      = dblarr(C,3)  
    mapB      = dblarr(npix,nmaps)
    der0_B    = dblarr(npix)     ;B
    der1_B_1  = dblarr(npix)     ;first derivatives  
    der1_B_2  = dblarr(npix)     ;first derivatives
    der2_B_1  = dblarr(npix)     ;second derivatives
    der2_B_2  = dblarr(npix)     ;second derivatives
    der2_B_3  = dblarr(npix)     ;second  derivatives
    der0_B_cut    = dblarr(C)     ;B
    der1_B_1_cut  = dblarr(C)     ;first derivatives  
    der1_B_2_cut  = dblarr(C)     ;first derivatives
    der2_B_1_cut  = dblarr(C)     ;second derivatives
    der2_B_2_cut  = dblarr(C)     ;second derivatives
    der2_B_3_cut  = dblarr(C)     ;second  derivatives
    v0d       = dblarr(25)       ;the first MF
    v1d       = dblarr(25)       ;the second MF
    v2d       = dblarr(25)       ;the third MF
    nu        = dblarr(25)       ;given thresholds
    u         = dblarr(4,25)
    
    ;map       = dblarr(npix)
    ;map1      = dblarr(npix)


    mu_B     = 0.0        ;the expectation value of B
    sigma_B  = 0.0        ;the variance of B
    tau_B    = 0.0  
    
    read_fits_map,strcompress('map_B_no_gauss_r0_s60_real_'+string(n)+'.fits',/remove),mapB,xdata,pdata
    ;ud_grade,map1,map,nside_out=nside, order_in='ring',order_out='ring'

    ;write_fits_map,'Bmap.fits',map,coordsys='G', ordering='ring'

    ;fits2alm,index,a1,strcompress('blm_10-'+string(n)+'.fits',/remove),HDR=hdr,XHDR=xhdr
    ;alm2fits,index,a1,'alms.fits',HDR=hdr,XHDR=xhdr
    

    
    ;spawn,"synfast syn.par"
    
    ;read_fits_map, 'B_der.fits', mapB,xdata,pdata
    ud_grade, mapB(*,0),der0_B,nside_out=nside, order_in='ring',order_out='ring'
    ud_grade, mapB(*,1),der1_B_1,nside_out=nside, order_in='ring',order_out='ring'
    ud_grade, mapB(*,2),der1_B_2,nside_out=nside, order_in='ring',order_out='ring'
    ud_grade, mapB(*,3),der2_B_1,nside_out=nside, order_in='ring',order_out='ring'
    ud_grade, mapB(*,4),der2_B_2,nside_out=nside, order_in='ring',order_out='ring'
    ud_grade, mapB(*,5),der2_B_3,nside_out=nside, order_in='ring',order_out='ring'

    ;mollview,der0_B
    ;stop
    ;write_fits_map,strcompress('b_map'+string(n)+'.fits',/remove), der0_B, coordsys=’G’, ordering=’ring’
    ;fits2alm,index,alm,'blm.fits',HDR=hdr,XHDR=xhdr
    ;alm2fits,index,alm,strcompress('blm-'+string(n)+'.fits',/remove),HDR=hdr,XHDR=xhdr

    pix2ang_ring ,nside, pixels,theta,phi



;begin to cut the sky:
    d=0
    for i=0.,npix-1.   do begin
        if w[i] ne 0.  then begin
            theta_cut[d]=theta[i]
            phi_cut[d]=phi[i]
            der0_B_cut[d]=der0_B[i]
            der1_B_1_cut[d]=der1_B_1[i]
            der1_B_2_cut[d]=der1_B_2[i]
            der2_B_1_cut[d]=der2_B_1[i]
            der2_B_2_cut[d]=der2_B_2[i]
            der2_B_3_cut[d]=der2_B_3[i]
            d=d+1
        endif
    endfor




;    for i=0.,npix-1.   do begin    ;Mink(0)=B
;       ;if w[i] eq 0.  then begin
;       ;    Mink[i,0]=-100000
;       ;endif
;	 if w[i] ne 0.  then begin
;    		Mink[i,0]=der0_B[i]*w[i]
;	 endif
;    endfor
;print,Mink(*,0)
;stop		

;    pix2ang_ring ,nside, pixels,theta,phi

    for i=0.,C-1.   do begin
      ;if sign[i] eq 1  then begin
        st = sin(theta_cut[i])
        if (st lt 0.000001)   then begin    ;make the st not equal zero during the numerical computations
             st = st+0.000001
        endif
        ct     = cos(theta_cut[i])
        ut_B   = der1_B_1_cut[i]
        up_B   = der1_B_2_cut[i]
        utt_B  = der2_B_1_cut[i]
        utp_B  = der2_B_2_cut[i]-ct/st*der1_B_2_cut[i]
        upp_B  = der2_B_3_cut[i]+ct/st*der1_B_1_cut[i]
        Mink[i,1] = sqrt(ut_B*ut_B+up_B*up_B)     ;from equation 7
        Mink[i,2] = (2.0*ut_B*up_B*utp_B-ut_B*ut_B*upp_B-up_B*up_B*utt_B)/(ut_B*ut_B+up_B*up_B)   ;from equation 7
        mu_B      = mu_B    + Mink[i,0]/C
        sigma_B   = sigma_B + Mink[i,0]*Mink[i,0]/C
        tau_B     = tau_B   + Mink[i,1]*Mink[i,1]/C
      ;endif
    endfor

    mu_B    = mu_B               ;calculating the mean of B
    sigma_B = sigma_B-mu_B*mu_B  ;calculating the variance of B
    tau_B   = tau_B*(1.0/2.0)    ;calculating the variance of the gradient of B

    for i=0,24  do begin
        nu[i]=i*0.25-3.0
    endfor


    for i=0.,C-1  do begin
       ;if sign[i] eq 1  then begin
            u_B = (Mink[i,0]-mu_B)/sqrt(sigma_B)     ;a new variable derivated by B which we use it as a smoothed scalar field
            if u_B gt -3.0 and u_B lt 2.999  then begin
                k1 = floor((u_B+3.0)/delta_nu)
                for j=0,k1   do begin
                    v0d[j] = v0d[j]+1.0/npix
                endfor
            endif
      
            if((u_B LT 3.12499) and (u_B GT -3.12499)) then begin    ;calculating v1d and v2d ---- equation 6&7
                k2=floor((u_B+3.125)/delta_nu)
                ;print,k2
                v1d[k2]=v1d[k2]+(1.0/npix)*(1.0/4.0)*(1.0/(delta_nu*sqrt(sigma_B)))*Mink[i,1]  ;*(4.*!PI)
                v2d[k2]=v2d[k2]+(1.0/npix)*(1.0/2.0/!PI)*(1.0/(delta_nu*sqrt(sigma_B)))*Mink[i,2]  ;*(4.*!PI)
	        endif
	    ;d=d+1
        ;endif  
    endfor
 ;print,d

 
    for i=0,24  do begin
        v[0,i]=nu[i]
        v[1,i]=v[1,i]+v0d[i]
        v[2,i]=v[2,i]+v1d[i]
        v[3,i]=v[3,i]+v2d[i]
        v[4,i]=v[4,i]+v0d[i]*v0d[i]
        v[5,i]=v[5,i]+v1d[i]*v1d[i]
        v[6,i]=v[6,i]+v2d[i]*v2d[i]
    endfor

	for i=0,24  do begin
		u[0,i]=nu[i]
		u[1,i]=v0d[i]
		u[2,i]=v1d[i]
		u[3,i]=v2d[i]
	endfor

    openw,1,strcompress('real_Non-gaussian_r0_s60_MF_'+string(n)+'.dat',/remove)
    printf,1,u
    close,1
endfor

for i=0,24   do begin
    v[0,i]=v[0,i]
    v[1,i]=v[1,i]/M  ;mean of the MF1
    v[2,i]=v[2,i]/M  ;mean of the MF2
    v[3,i]=v[3,i]/M  ;mean of the MF3
    v[4,i]=v[4,i]/M-v[1,i]*v[1,i]   ;variance of the MF1
    v[5,i]=v[5,i]/M-v[2,i]*v[2,i]   ;variance of the MF2
    v[6,i]=v[6,i]/M-v[3,i]*v[3,i]   ;variance of the MF3
endfor

openw,1,'real_Non-gaussian_r0_s60_MF.dat'
;openw,1,'real_gaussian_r0_s60_MF.dat'
;openw,1,'ideal_Non-gaussian_r0_s60_MF.dat'
;openw,1,'ideal_gaussian_r0_s60_MF.dat'
printf,1,v
close,1

end
