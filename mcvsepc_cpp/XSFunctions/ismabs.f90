! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ISMABS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! XSPEC local model for ISM absorption
! Version 1.2 May 2015
! (Without Ni and Zn in this version)
!
! Additions to version 1.2
! - The parameter names do not have mathematical operators.
! - The subroutine names in the fortran code have been changed.
! - The lmodel.dat file has been rename to lmodel_ismabs.dat
!
!
! Additions to version 1.1
! - the model now prints a message to STDOUT if it is unable to
! read from the AtomicData.fits file
! - the ISMABSROOT xset variable has been added to allow the
! location of the AtomicData.fits file to be changed from within
! X-Spec
! - the startup does not depend on the ifl parameter
!
! To-Do:
! - Add turbulence
! - Add molecular-solid cross sections
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ismabs(ear, ne, param, ifl, photar)
!
! The main routine to call all subroutines
!
implicit none
integer,parameter :: num_param = 31, out_unit=20, nion=30
integer,parameter :: nemod=650000 !Number of tabulated energies for each cross section.
integer :: ne, ifl, i
double precision,allocatable :: emod(:), coeff(:), cion(:)
double precision,allocatable,save :: tabxs(:,:), tabener(:)
double precision :: nH
double precision :: zfac
double precision :: ear(0:ne), param(num_param),photar(ne)
logical, save :: startup=.true.
character (len=40) version

! silly line to suppress a compiler warning. note that nH is reset below
nH = dfloat(ifl)

version='1.2'
if(startup)then
   print *, ' '
   print *, 'ISMabs: ISM absorption model Version',version
   print *, 'Gatuzz, Garcia, Kallman, Mendoza, & Gorczyca (2014)'
   print *, 'Note: Default column densities are given'
   print *, 'according to Grevesse, N. & Sauval (1998)'     
   print *, 'assuming N_H = 1.E21 cm^-2'  
   print *, ' '
   ! allocate memory for saved arrays and load them
   allocate(tabxs(nemod,0:nion), tabener(nemod))
   call read_cross_sections_ismabs(nemod,tabener,tabxs)
   startup=.false.  
endif
! Model parameters
allocate(cion(0:nion))

cion(0) = param(1)*1.0d22
cion(1) = 0.1*cion(0)
do i = 2, nion
   cion(i) = param(i)*1.0d16
enddo
zfac = 1/(1.d0+param(31))


! allocate memory for temporary arrays
allocate(emod(0:nemod), coeff(nemod))

! calculate the absorption (ie exp(-tau)) at energies emod and place in the
! coeff array
call absorption_ismabs(cion,nion,zfac,emod,nemod,coeff,tabxs,tabener)

! interpolate the coeff array onto photar using the routine defined in
! XSFunctions/Utilities/xsCFortran.c.
call dlinintinteg(nemod,emod,coeff,(ne+1),ear,photar,0.0d0,1.0d0)
! divide the binsize because linintinteg does an integratio over bins
do i = 1, ne
   photar(i) = photar(i) / (ear(i) - ear(i-1))
enddo

! deallocate memory for temporary arrays
deallocate(emod, coeff, cion)

return
end subroutine ismabs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_cross_sections_ismabs(nemod,tabener,tabxs)
!
! This routine reads all cross sections and puts them on a given grid
!
! The IsmabsAtomicData.fits file is in the standard directory unless
! ISMABSROOT is set.
!
implicit none
integer,parameter :: nion=30, out_unit=20
integer :: nemod, i, status
double precision :: tabener(nemod), tabxs(nemod,0:nion)
character (*), parameter :: fileloc = '/IsmabsAtomicData.fits'
character (*), parameter :: ismreadchat = 'ismabs: reading from '
character (len=255 + 29) :: filename2 ! len(fileloc)
character (len=255) :: ismabs_root = ''
character (len=len(ismreadchat)+len(filename2)) :: chatmsg = ''
integer inunit,readwrite,blocksize
integer :: hdutype,colnum
integer :: felem=1, nulld=0
logical :: anynull
character (len=255) :: fgmstr, fgmodf
external :: fgmstr, fgmodf

! Where do we look for the data?
ismabs_root = trim(fgmstr('ISMABSROOT'))
if (ismabs_root .EQ. '') then
   ismabs_root = trim(fgmodf())
endif
! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ismabs_root) // fileloc
chatmsg=ismreadchat // filename2
call xwrite(chatmsg, 20)
! Get an unused Logical Unit Number to use to open the FITS file.
call ftgiou(inunit,status)
! Open the FITS file
call ftopen(inunit,filename2,readwrite,blocksize,status)
! Move to the extension 2 (the binary table)
call ftmahd(inunit,2,hdutype,status)
!Read the energy grid (column 1)
colnum=1
call ftgcvd(inunit,colnum,1,felem,nemod,nulld,tabener,anynull,status)

!Read cross sections
colnum=2
do i=0, nion
   call ftgcvd(inunit,colnum,1,felem,nemod,nulld,tabxs(1,i),anynull,status)
   colnum=colnum+1
enddo
! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the X-Spec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used. An alternative would be to use xwrite()
! with a low chatter level.
!
! This message could be displayed only once, but it is probaly worth
! repeating each time it is used.
if (status .ne. 0) then
   write (*,*) 'ERROR: unable to read cross sections from ', filename2
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_cross_sections_ismabs


! ======================================= !
subroutine absorption_ismabs(cion,nion,zfac, emod, nemod, coeff, tabxs, tabener)
!
! This is routine that calculates the optical depth given the column densities
! Finally returns the absorption coefficient exp(-tau)
!
implicit none
integer,parameter :: out_unit=20
integer :: nemod, nion
integer :: i, j
double precision :: tabener(nemod), tabxs(nemod,0:nion), emod(nemod)
double precision :: coeff(nemod), cion(0:nion)
double precision,allocatable :: tau(:)
double precision :: zfac
real hphoto, gphoto
external hphoto, gphoto

allocate(tau(nemod))

! Calculates the optical depth and the absorption coefficient exp(-tau)
do i=1,nemod
   emod(i)=(tabener(i)*zfac)/1.d3
   ! In case you want to read xspec hydrogen column density
   ! tabxs(0,i)=hphoto(real(emod(i-1)),real(emod(i)))
   ! tau for hydrogen column density
   tau(i)=cion(0)*tabxs(i,0)
   ! Calculates the optical depth and the absorption coefficient exp(-tau)
   tau(i)= tau(i)+(tabxs(i,1)*cion(1)) ! He I column density = 0.1 Nh
enddo

do j=2,nion
   do i=1, nemod
      tau(i)=tau(i)+(cion(j)*tabxs(i,j))
   enddo
enddo
do i=1, nemod
   coeff(i)=dexp(-tau(i))
enddo
deallocate(tau)

end subroutine absorption_ismabs
