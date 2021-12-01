pro MF_FUNC_B

nside    = 512.
npix     = 12.*nside*nside
delta_nu = 0.25       ;a discretization of threshold space in bins of width delta_nu
nmaps    = 6
mu_B     = 0.0        ;the expectation value of B
sigma_B  = 0.0        ;the variance of B
tau_B    = 0.0        ;the variance of the gradient of B

pixels    = lindgen(npix)
theta     = dblarr(npix)
phi       = dblarr(npix)
Mink      = dblarr(npix,3)   ;three formulas contained by covariant differentiation of B
mapB      = dblarr(npix,nmaps)
der0_B    = dblarr(npix)     ;B
der1_B_1  = dblarr(npix)     ;first derivatives  
der1_B_2  = dblarr(npix)     ;first derivatives
der2_B_1  = dblarr(npix)     ;second derivatives
der2_B_2  = dblarr(npix)     ;second derivatives
der2_B_3  = dblarr(npix)     ;second  derivatives
v0d       = dblarr(25)       ;the first MF
v1d       = dblarr(25)       ;the second MF
v2d       = dblarr(25)       ;the third MF
nu        = dblarr(25)       ;given thresholds
v          =dblarr(4,25)

B_der='smoothed_mapTEB.fits'
read_fits_map, B_der, mapB,xdata,pdata
ud_grade, mapB(*,12),der0_B,nside_out=nside, order_in='ring',order_out='ring'
ud_grade, mapB(*,13),der1_B_1,nside_out=nside, order_in='ring',order_out='ring'
ud_grade, mapB(*,14),der1_B_2,nside_out=nside, order_in='ring',order_out='ring'
ud_grade, mapB(*,15),der2_B_1,nside_out=nside, order_in='ring',order_out='ring'
ud_grade, mapB(*,16),der2_B_2,nside_out=nside, order_in='ring',order_out='ring'
ud_grade, mapB(*,17),der2_B_3,nside_out=nside, order_in='ring',order_out='ring'
;print,xdata
;print,pdata
;stop

for i=0.,npix-1.   do begin    ;Mink(0)=B
   Mink[i,0]=der0_B[i]
endfor

pix2ang_ring ,nside, pixels,theta,phi

for i=0.,npix-1.   do begin
  st = sin(theta[i])
  if (st lt 0.000001)   then begin    ;make the st not equal zero during the numerical computations
     st = st+0.000001
  endif
  ct     = cos(theta[i])
  ut_B   = der1_B_1[i]
  up_B   = der1_B_2[i]
  utt_B  = der2_B_1[i]
  utp_B  = der2_B_2[i]-ct/st*der1_B_2[i]
  upp_B  = der2_B_3[i]+ct/st*der1_B_1[i]
  Mink[i,1] = sqrt(ut_B*ut_B+up_B*up_B)     ;from equation 7
  Mink[i,2] = (2.0*ut_B*up_B*utp_B-ut_B*ut_B*upp_B-up_B*up_B*utt_B)/(ut_B*ut_B+up_B*up_B)   ;from equation 7
  mu_B      = mu_B    + Mink[i,0]/npix
  sigma_B   = sigma_B + Mink[i,0]*Mink[i,0]/npix
  tau_B     = tau_B   + Mink[i,1]*Mink[i,1]/npix
endfor

mu_B    = mu_B               ;calculating the mean of B
sigma_B = sigma_B-mu_B*mu_B  ;calculating the variance of B
tau_B   = tau_B*(1.0/2.0)    ;calculating the variance of the gradient of B

for i=0,24  do begin
  nu[i]=i*0.25-3.0
endfor

for i=0.,npix-1  do begin
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
endfor

for i=0,24  do begin
    v[0,i]=nu[i]
    v[1,i]=v0d[i]
    v[2,i]=v1d[i]
    v[3,i]=v2d[i]
endfor


openw,1,"MF_FUNC_B.dat"
printf,1,v
close,1

end