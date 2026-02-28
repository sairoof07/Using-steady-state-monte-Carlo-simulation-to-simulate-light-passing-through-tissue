program project 
implicit none 

integer:: Nbins, N
double precision, parameter:: pi=4*atan(1.0d0),c0=2.997925d10,cos90=1.0d-6
double precision:: threshold, chance,Pstatus
 !! logicals 
 integer:: Alive=1, dead=0


double precision:: x,y,z, ux,uy,uz, uxx,uyy,uzz,rp,Rarm,rb,rmuscle,Rfat !! photon's trajectory and position 
double precision:: s,psi !! step size randomnum
double precision :: costheta, sintheta, cospsi, sinpsi !! sin fimctions

integer:: Iphoton
double precision:: ifoton, w , absorb

double precision, allocatable:: Csph(:), Ccyl(:), Cplan(:) !!photon concentration
double precision:: fSph, fcyl, fplan 
double precision:: mua, mus,g,mut,nt,albedo !! optical properties 
integer:: Nphotons, Nr
double precision:: r_size, r,dr 
integer:: ir 
 double precision:: shellV, CNT 
  !! dummy variables 
  double precision:: rnd ,u ,temp,a
  integer:: i,j 
  real:: start, finish 
  
  call cpu_time(start)
  
  !! input
  mua = 0.0
mus = 0.0
g = 0.90d0
nt = 1.33d0
Nphotons = 10000
r_size = 6.0
Nbins=500
Nr=Nbins
chance=0.1d0
threshold=0.01d0

!! arm radius 
 Rarm=4.52
 Rfat=4.384
 Rmuscle=4.126 
dr=r_size/Nbins
albedo= mus/(mus+mua)
allocate(Csph(0:Nr), Ccyl(0:Nr), Cplan(0:nr))

!! initialization 
Csph(:) = 0
Ccyl(:) = 0 
Cplan(:) = 0


r=0
 
i=1

do while(i<Nphotons) !! loop for all photons 

 
if (mod(i,Nphotons/10).eq.0) then !! progress indicator
write(0,*) i*100/Nphotons,'% done'
end if
 W = 1.0d0
pstatus = ALIVE
!! beam's initialization 
x = pi/2*sqrt(-log(rand()))
y = 0
z = 0.01d0

costheta = 2.0d0*Rand() - 1.0d0
sintheta = sqrt(1.0d0 - costheta**2.0d0)
psi = 2.0d0*PI*Rand()
ux = 0.0d0
uy = 0.0d0
uz = 1.0d0
do while(pstatus==Alive)

rp=sqrt(y**2+(z-Rarm)**2)  !! the position of the photon with respect to the center of the arm
rb=sqrt((z-1.18)**2+(y-5.65)**2) !! position of the photon with respect to the center of the bone
if (rp>Rarm) then  
albedo=0 !! outside the tissue 
else 
if (rp<Rarm) then  !! inside the skin
mus=268.01 
mua= 11.0
end if 
if (rp<rfat) then !! inside the fat tissue
 mus=157.03 
 mua=0.15
end if 
if (rp<rmuscle) then !! inside muscle tissue 
 mus=255.41 
 mua=0.2
end if 

if (z>0.7  .and. z<1.79) then !! in the blood section

if (y>0.25 .and. y<0.75) then
 mus=188.3 
 mua= 11.0
end if 
end if 

if (rb < 1.15) then !! in the bones 
 mus=193.4 
 mua=20.0
end if
albedo= mus/(mus+mua) !! calculating the albedo 
end if 




!! hop
a=rand() 
s = -log(a)/(mua + mus)
x =x+ s * ux 
y =y+ s * uy
z =z+ s * uz

absorb = W*(1.0d0 - albedo) !! drop
w= w-absorb

r = sqrt(x**2.0d0 + y**2.0d0 + z**2.0d0)
ir = ceiling(r/dr)

Csph(ir) =Csph(ir)+ absorb/mua

r = sqrt(x**2.0d0 + y**2.0d0)
ir = ceiling(r/dr)
if (ir >= NR) then
ir = NR
endif
Ccyl(ir) =Ccyl(ir)+ absorb

r = abs(z)
ir = ceiling(r/dr)
if (ir >= NR) then 
ir =Nr
endif
Cplan(ir)=Cplan(ir) + absorb



temp = (1.0d0 - g**2.0d0)/(1.0d0 - g + 2*g*rand())
costheta = (1.0d0 + g**2.0d0 - temp**2.0d0)/(2.0d0*g)

sintheta=sqrt(1.0d0-costheta**2.0d0)
psi = 2.0*PI*Rand()
cospsi = cos(psi)

if (psi < PI) then 
sinpsi = sqrt(1.0d0 - cospsi**2.0d0) 
else
sinpsi = -sqrt(1.0d0 - cospsi**2.0d0)
endif


if (1.0d0 - abs(uz) <= 1.0d-12) then !! spin 
uxx = sintheta * cospsi
uyy = sintheta * sinpsi
uzz = costheta *sign(1.0d0,2.0d0*uz-1.0d0)

else 
temp = sqrt(1.0 - uz * uz)
uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
uzz = -sintheta * cospsi * temp + uz * costheta
end if 

ux = uxx
uy = uyy
uz = uzz
!! monte carlo 
if (W < THRESHOLD) then  
a=rand()
if (a.le.CHANCE) then
w = w/CHANCE

else 
pstatus = DEAD

endif 
endif
write(19,*) x,y,z
enddo 
i=i+1
 enddo
write(0,*)  '         100% done'


!! the output
open (file='results.txt', unit=11)
do ir=1,Nr 
r = (ir + 0.5d0)*dr
shellV = 4.0d0*dr*PI*r**2.0d0
Fsph = Csph(ir)/(Nphotons*shellv)

write(11,*) r, fsph,ir
enddo
call cpu_time(finish)

 print*, finish-start

deallocate(Csph, Ccyl, Cplan)



end program 
