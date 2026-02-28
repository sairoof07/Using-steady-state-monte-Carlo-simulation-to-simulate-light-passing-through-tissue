program project 
implicit none 

integer:: Nbins, N
double precision, parameter:: pi=4*atan(1.0d0),c0=2.997925d10,cos90=1.0d-6
double precision:: threshold, chance
 !! logicals 
 integer:: Alive=1, dead=0


double precision:: x,y,z, ux,uy,uz, uxx,uyy,uzz 
double precision:: s,psi !! step size randomnum
double precision :: costheta, sintheta, cospsi, sinpsi

integer:: Iphoton
double precision:: ifoton, w , absorb
!! weight
integer:: Pstatus
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
  
  
  !! input
  mua = 0.2d0
mus = 255.41d0
g = 0.720d0
nt = 1.33d0
Nphotons = 5000
r_size = 6.0
Nbins=500
Nr=Nbins
chance=0.1d0
threshold=0.01d0

dr=r_size/Nbins
albedo= mus/(mus+mua)
allocate(Csph(0:Nr), Ccyl(0:Nr), Cplan(0:nr))


Csph(:) = 0
Ccyl(:) = 0 
Cplan(:) = 0


r=0
 
i=1

do while(i<Nphotons)

 
if (mod(i,Nphotons/10).eq.0) then 
write(0,*) i*100/Nphotons,'% done'
end if
 W = 1.0d0
pstatus = ALIVE
x = sign(rand()*1.0d0,2.0d0*rand()-1.0d0)
y = 0.0d0
z = 0.01d0

costheta = 2.0d0*Rand() - 1.0d0
sintheta = sqrt(1.0d0 - costheta**2.0d0)
psi = 2.0d0*PI*Rand()
ux = 0
uy = 0
uz = 1.0d0
do while(pstatus==Alive)


a=rand() 
s = -log(a)/(mua + mus)
x =x+ s * ux 
y =y+ s * uy
z =z+ s * uz

absorb = W*(1.0d0 - albedo)
w= w-absorb

r = sqrt(x**2.0d0 + y**2.0d0 + z**2.0d0)
ir = floor(r/dr)

Csph(ir) =Csph(ir)+ absorb

r = sqrt(x**2.0d0 + y**2.0d0)
ir = int(r/dr)
if (ir >= NR) then
ir = NR
endif
Ccyl(ir) =Ccyl(ir)+ absorb

r = abs(z)
ir = int(r/dr)
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


if (1.0d0 - abs(uz) <= 1.0d-12) then 
uxx = sintheta * cospsi
uyy = sintheta * sinpsi
uzz = costheta *sign(1.0d0,2.0d0*uz-1.0d0)
write(09,*) uzz
else 
temp = sqrt(1.0 - uz * uz)
uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
uzz = -sintheta * cospsi * temp + uz * costheta
end if 

ux = uxx
uy = uyy
uz = uzz

if (W < THRESHOLD) then  
a=rand()
if (a.le.CHANCE) then
w = w/CHANCE

else 
pstatus = DEAD

endif 
endif

enddo 
i=i+1
 enddo
write(0,*)  '         100% done'



open (file='results.txt', unit=11)
do ir=1,Nr 
r = (ir + 0.5d0)*dr
shellV = 4.0d0*dr*PI*r**2.0d0
Fsph = Csph(ir)/(Nphotons*shellv*mua)

write(11,*) r, fsph,ir
enddo


 

deallocate(Csph, Ccyl, Cplan)



end program 
!!end program 
! boundary conditions 
! if r>rf then mu=muskin (with respect to 0,0)
! if r>rm then mu=mufat
! if z>0 and z<L then 
! if y>L2 and y<L3 then mu=mublood
! end if 
! if r>rbonem then mu=mubone (with respect to the second zero)
! if r< rbonem then mu=mubonem