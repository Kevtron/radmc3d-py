;
; Some natural constants
;
AU  = 1.49598d13     ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]
MH  = 1.67372d-24    ; 
KB  = 1.3806488d-16  ;
;
; Monte Carlo parameters
;
nphot    = 1000000
;
; Grid parameters
;
;nx       = 32L
;ny       = 32L
;nz       = 32L
;sizex    = 10*AU
;sizey    = 10*AU
;sizez    = 10*AU
;
openr,lun,'MASS001_density_g_cm3_grid.dat', /GET_LUN
header=STRARR(19)
readf,lun,header
junk='' ; string variable 
nx=0L
ny=0L
nz=0L
reads,header(18),junk,nx,junk,ny,junk,nz,format='(A8,x,I0,A8,x,I0,A8,x,I0)'
array=FLTARR(nx,ny*nz)
readf,lun,array
rhod=reform(array,[nx,ny,nz]) 
sizex    = 1.12*RS
sizey    = 1.12*RS
sizez    = 1.05*RS
close,lun
FREE_LUN,lun
;
; Model parameters
;
;radius   = 5*AU
;rho0     = 1d-16
;
; Star parameters
;
mstar    = MS
rstar    = RS
tstar    = TS
;pstar    = [0.,0.,0.]
;
; Make the coordinates
;
xi       = -sizex + 2*sizex*dindgen(nx+1)/(1.d0*nx)
yi       = -sizey + 2*sizey*dindgen(ny+1)/(1.d0*ny)
zi       = -sizez + 2*sizez*dindgen(nz+1)/(1.d0*nz)
;
; Write the wavelength_micron.inp file
;
lambda1 = 0.1d0
lambda2 = 7.0d0
lambda3 = 25.d0
lambda4 = 1.0d4
n12     = 20
n23     = 100
n34     = 30
lam12   = lambda1 * (lambda2/lambda1)^(dindgen(n12)/(1.d0*n12))
lam23   = lambda2 * (lambda3/lambda2)^(dindgen(n23)/(1.d0*n23))
lam34   = lambda3 * (lambda4/lambda3)^(dindgen(n34)/(1.d0*(n34-1.d0)))
lambda  = [lam12,lam23,lam34]
nlam    = n_elements(lambda)
;
; Temp
;
openr,lun,'MASS001_u_erg_g_grid.dat', /GET_LUN
vol=(sizex/nx)*(sizey/ny)*(sizez/nz)
header=STRARR(19)
energyarray=FLTARR(nx,ny*nz)
readf,lun,header
readf,lun,energyarray
energy=reform(energyarray,nx,ny,nz)
Temp = (2*energy*MH)/(3*KB)
close,lun
FREE_LUN, lun

openw,1,'dust_temperature.dat'
printf,1,1
printf,1,nx*ny*nz 
printf,1,1
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,Temp[ix,iy,iz]
      endfor
   endfor
endfor
printf,1,''
close,1
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,0                      ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,nx,ny,nz               ; Size of grid
for i=0,nx do printf,1,xi[i]    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i]    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhod[ix,iy,iz]
;	 printf,1,0
      endfor
   endfor
endfor
close,1
;
; ions
;
;ion=(1/MH)*rhod
;openw,1,'ion_numdens.inp'
;printf,1,1                      ; Format number
;printf,1,nx*ny*nz               ; Nr of cells
;printf,1,1                      ; Nr of dust species
;for iz=0,nz-1 do begin
;   for iy=0,ny-1 do begin
;      for ix=0,nx-1 do begin
;         printf,1,ion[ix,iy,iz]
;	  printf,1,0 
;     endfor
;   endfor
;endfor
;close,1
;print,ion[nx/2,ny/2,nz/2]
;
; electrons
;
;openw,1,'electron_numdens.inp'
;printf,1,1                      ; Format number
;printf,1,nx*ny*nz               ; Nr of cells
;printf,1,1                      ; Nr of dust species
;for iz=0,nz-1 do begin
;   for iy=0,ny-1 do begin
;      for ix=0,nx-1 do begin
;         printf,1,ion[ix,iy,iz]
;;          printf,1,0 
;     endfor
;   endfor
;endfor
;close,1
;
; Star stuff
;
;
openw,1,'stellarsrc_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhod[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = 0'
printf,1,'iranfreqmode = 1'
;printf,1,'incl_freefree = 1'
close,1
;

;
; Dust opacity control file
;
openw,1,'dustopac.inp'
printf,1,'2               Format number of this file'
printf,1,'1               Nr of dust species'
printf,1,'============================================================================'
printf,1,'1               Way in which this dust species is read'
printf,1,'0               0=Thermal grain'
printf,1,'starstuff        Extension of name of dustkappa_***.inp file'
printf,1,'----------------------------------------------------------------------------'
close,1
;
; Dust Opacity File
;

openr,lun,'../../marcs-opacity/4000_wavelength.dat', /GET_LUN
wavelengtharray=FLTARR(10,107)
readf,lun,wavelengtharray
wavelength=reform(wavelengtharray/10000,1070)
close,lun
FREE_LUN, lun
openr,lun,'../../marcs-opacity/4000_opa.dat', /GET_LUN
opacityarray=FLTARR(6,357)
readf,lun,opacityarray
opacity=reform(opacityarray,2,3*357)
close,lun
FREE_LUN,lun
totalopa=2.2727d-03
openw,1,'dustkappa_starstuff.inp'
printf,1,2
printf,1,1070
for i=0,1069 do begin
    printf,1,wavelength[i],opacity[0,i]*totalopa,opacity[1,i]*totalopa
endfor
close,1
;
; Write the wavelength file
;
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
close,1
;
; Stellar Source template.
;
openw,1,'stellarsrc_templates.inp'
printf,1,2
printf,1,1
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,-4000
printf,1,RS
printf,1,MS
printf,1,''
close,1
;
; Movie control file
;
openw,1,'movie.inp'
printf,1,1
printf,1,360
for i=1,360 do begin
   printf,1, 0.,0.,0.,2*sizex,2*sizey,0.,0.,float(i)
endfor
close,1


end


