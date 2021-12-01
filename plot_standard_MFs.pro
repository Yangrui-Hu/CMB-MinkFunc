pro plot_standard_MFs

nu       =dblarr(25)
mu       =0.0
sigma    =1.0
tau      =0.0153147
v0g      =dblarr(25)
v1g      =dblarr(25)
v2g      =dblarr(25)
delta_nu =0.25
v        =dblarr(4,25)

for i=0,24  do begin
	nu[i]=i*0.25-3.0
	v0g[i]=1.0/2.0*(1-erf((nu[i]-mu)/(sqrt(2)*sigma)))
	v1g[i]=1.0/8.0/delta_nu*sqrt(tau*!PI/2)*(erf((nu[i]-mu+delta_nu/2)/(sqrt(2)*sigma))-erf((nu[i]-mu-delta_nu/2)/(sqrt(2)*sigma)))
	v2g[i]=tau/(2*!PI*sqrt(2*!PI)*sigma*delta_nu)*(exp(-(nu[i]-mu-delta_nu/2)*(nu[i]-mu-delta_nu/2)/(2*sigma*sigma))-exp(-(nu[i]-mu+delta_nu/2)*(nu[i]-mu+delta_nu/2)/(2*sigma*sigma)))
	v[0,i]=nu[i]
             v[1,i]=v0g[i]
             v[2,i]=v1g[i]
             v[3,i]=v2g[i]
endfor

openw,1,"Standard_MFs.dat"
printf,1,v
close,1

end