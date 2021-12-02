!######## This file contains various routines for dealing with grid data


!!!------ Define property/origin/spacing/grid number and then save to a 3D matrix, infomode=1 means silent
!! iorb is used to choose the orbital for whose wavefunction will be calculated. This can be an arbitrary value if functype/=4
subroutine savecubmat(functype,infomode,iorb)
use defvar
use util
use function
implicit real*8 (a-h,o-z)
integer :: infomode,functype,iorb !Calculate which orbital wavefunction for fmo routine
real*8 xarr(nx),yarr(ny),zarr(nz)
character c80tmp*80,c200tmp*200,c400tmp*400,filename_tmp*200
integer*2,allocatable :: corpos(:,:,:)
logical,allocatable :: boundgrd(:,:,:)

!---- Special case, use cubegen to directly evaluate ESP grid data
alive=.false.
if (cubegenpath/=" ".and.ifiletype==1.and.functype==12) then
	inquire(file=cubegenpath,exist=alive)
	if (alive==.false.) then
		write(*,"(a)") " Note: Albeit current file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been defined, &
		the cubegen cannot be found, therefore electrostatic potential will still be calculated using internal code of Multiwfn"
	end if
end if
if (alive.and.ifiletype==1.and.functype==12) then !Use cubegen to calculate total ESP
	call walltime(iwalltime1)
	write(*,"(a)") " Since the input file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been properly defined, &
	now Multiwfn directly invokes cubegen to calculate electrostatic potential"
	
	!Generate cubegen input file
	open(10,file="ESPgridtmp.cub",status="replace")
	write(10,"(' Generated by Multiwfn')")
	write(10,"(' Totally ',i12,' grid points')") nx*ny*nz
	write(10,"(i5,3f12.6)") ncenter,orgx,orgy,orgz
	write(10,"(i5,3f12.6)") nx,gridv1
	write(10,"(i5,3f12.6)") ny,gridv2
	write(10,"(i5,3f12.6)") nz,gridv3
	close(10)
	ncubegenthreads=1 !Parallel implementation prior to G16 is buggy, so test here
	if (index(cubegenpath,"G16")/=0.or.index(cubegenpath,"g16")/=0) ncubegenthreads=nthreads
	filename_tmp=filename
	if (index(filename,".chk")/=0) call chk2fch(filename_tmp)
	write(c400tmp,"(a,i5,a)") """"//trim(cubegenpath)//"""",ncubegenthreads," potential="//trim(cubegendenstype)//" "//&
	""""//trim(filename_tmp)//""""//" ESPresult.cub -1 h ESPgridtmp.cub > nouseout"
	call runcommand(c400tmp)
	if (index(filename,".chk")/=0) call delfile(filename_tmp)
    
	!Load ESP data from cubegen resulting file
	call readcube("ESPresult.cub",1,1)
	!Delete intermediate files
    call delfile("cubegenpt.txt ESPresult.cub ESPgridtmp.cub nouseout")
	call walltime(iwalltime2)
	if (infomode==0) write(*,"(' Calculation of grid data took up wall clock time',i10,' s')") iwalltime2-iwalltime1
    
	return
end if

!--- Another special case, use slow but specifically optimized code for evaluating ESP grid data (deprecated)
if (functype==12.and.iESPcode==1) then
    call cubesp
    return
end if

!--- Below are normal cases, only use Multiwfn regular internal code
iorbsel=iorb
if (infomode==0.and.functype/=12) then !Not ESP case
	if (ifPBC>0) then
		if (expcutoff_PBC<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored ')") expcutoff_PBC
	else
		if (expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored ')") expcutoff
	end if
end if

!Writing and then reading can cut the minimal noise at the end of the coordinate, otherwise the originally symmetry points may become unsymmetry
if (ifPBC==0) then
    do k=1,nz
	    write(c80tmp,"(D20.13)") orgz+(k-1)*dz
	    read(c80tmp,*) zarr(k)
    end do
    do j=1,ny
	    write(c80tmp,"(D20.13)") orgy+(j-1)*dy
	    read(c80tmp,*) yarr(j)
    end do
    do i=1,nx
	    write(c80tmp,"(D20.13)") orgx+(i-1)*dx
	    read(c80tmp,*) xarr(i)
    end do
end if

!When ESPrhoiso/=0, only calculate ESP for grid around isosurface of rho=ESPrhoiso to save time
if (functype==12.and.ESPrhoiso/=0) then
	allocate(corpos(nx,ny,nz),boundgrd(nx,ny,nz))
    write(*,"(' Note: ESP will be calculated only for the grids around isosurface of electron density of ',f10.6,' a.u.')") ESPrhoiso
    write(*,*) "Detecting the grids for calculating ESP..."
	!$OMP PARALLEL DO SHARED(corpos) PRIVATE(ix,iy,iz) schedule(dynamic) NUM_THREADS(nthreads)
    do iz=1,nz
		do iy=1,ny
			do ix=1,nx
				if (fdens(xarr(ix),yarr(iy),zarr(iz))>=ESPrhoiso) then !internal corner
					corpos(ix,iy,iz)=0
				else !external corner
					corpos(ix,iy,iz)=1
				end if
			end do
		end do
	end do
	!$OMP END PARALLEL DO
    boundgrd=.false.
    ndet=ESPrhonlay !The thickness of detecting grid
    ndetmax=(2*ndet+1)**3-1
    do iz=1+ndet,nz-ndet
		do iy=1+ndet,ny-ndet
			do ix=1+ndet,nx-ndet
				!Cycle surrounding ndetmax grids
				icubtest=0
				do k=-ndet,ndet
					do j=-ndet,ndet
						do i=-ndet,ndet
							if (i==0.and.j==0.and.k==0) cycle
							icubtest=icubtest+corpos(ix+i,iy+j,iz+k)
						end do
					end do
                end do
                !At least one grid among the detection grids has different property to others
				if (icubtest>0.and.icubtest<ndetmax) boundgrd(ix,iy,iz)=.true.
			end do
		end do
	end do
	write(*,"(' Number of grids to calculate ESP:',i12)") count(boundgrd==.true.)
end if

call walltime(iwalltime1)

!If the function to be calculated is related to ESP, initialize LIBRETA so that faster code will be used
nthreads_old=nthreads
if (ifdoESP(functype).and.iESPcode==2) then
    call doinitlibreta
    if (isys==1.and.nthreads>10) nthreads=10
end if

!Start calculation of grid data
if (infomode==0) call showprog(0,nz)
ifinish=0
cubmat=0
!$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,densval) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
            if (ifpbc==0) then
			    tmpx=xarr(i)
                tmpy=yarr(j)
    				tmpz=zarr(k)
            else !In the case of PBC, calculate x,y,z of points according grid vector
                call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            end if
			if (functype==1513) then !Only involved by funcvsfunc routine, when RDG and sign(lambda2)rho is combined
				call signlambda2rho_RDG(tmpx,tmpy,tmpz,cubmat(i,j,k),cubmattmp(i,j,k))
			else if (functype==1614) then !Only involved by funcvsfunc routine, when promolecular RDG and sign(lambda2)rho is combined
				call signlambda2rho_RDG_prodens(tmpx,tmpy,tmpz,cubmat(i,j,k),cubmattmp(i,j,k))
            else if (functype==1599) then !Only involved by funcvsfunc routine, when IRI and sign(lambda2)rho is combined
                call IRI_s2lr(tmpx,tmpy,tmpz,cubmat(i,j,k),cubmattmp(i,j,k))
			else
				if (functype==12.and.ESPrhoiso/=0) then
					if (boundgrd(i,j,k)==.false.) cycle
                end if
				cubmat(i,j,k)=calcfuncall(functype,tmpx,tmpy,tmpz)
			end if
		end do
	end do
	if (infomode==0) then
        ifinish=ifinish+1
        call showprog(ifinish,nz)
	end if
end do
!$OMP END PARALLEL DO
if (infomode==0.and.ifinish<nz) call showprog(nz,nz)
nthreads=nthreads_old

!Set ESP of the grids neighouring to boundary grids to the ESP of neighouring grid when ESPrhonlay=1, this avoids weird color of ESP on vdW surface in VMD
if (functype==12.and.ESPrhoiso/=0.and.ESPrhonlay==1) then
	write(*,*) "Setting ESP of the grids neighbouring to boundary grids..."
    do iz=2,nz-1
		do iy=2,ny-1
			do ix=2,nx-1
				if (boundgrd(ix,iy,iz)==.true.) cycle
                if (boundgrd(ix+1,iy,iz)==.true.) cubmat(ix,iy,iz)=cubmat(ix+1,iy,iz)
                if (boundgrd(ix-1,iy,iz)==.true.) cubmat(ix,iy,iz)=cubmat(ix-1,iy,iz)
                if (boundgrd(ix,iy+1,iz)==.true.) cubmat(ix,iy,iz)=cubmat(ix,iy+1,iz)
                if (boundgrd(ix,iy-1,iz)==.true.) cubmat(ix,iy,iz)=cubmat(ix,iy-1,iz)
                if (boundgrd(ix,iy,iz+1)==.true.) cubmat(ix,iy,iz)=cubmat(ix,iy,iz+1)
                if (boundgrd(ix,iy,iz-1)==.true.) cubmat(ix,iy,iz)=cubmat(ix,iy,iz-1)
			end do
		end do
	end do
end if

if (infomode==0) then
    call walltime(iwalltime2)
    write(*,"(' Calculation of grid data took up wall clock time',i10,' s')") iwalltime2-iwalltime1
end if
end subroutine



!!------------------ Set up grid
!imode=0: Do not show the option used to load external points
!imode=1: Show the option used to load external points
!imode=2: Directly enter the option 9
!  igridsel is returned variable, corresponding to the selected index; if igridsel==100, that means user didn't set up grid here &
!but choose to load a set of point coordinates from external plain text file
!Usual calling instance: call setgrid(1,inouse)
subroutine setgrid(imode,igridsel)
use defvar
use GUI
use util
implicit real*8 (a-h,o-z)
real*8 molxlen,molylen,molzlen,tmpx,tmpy,tmpz
character*200 cubefilename,pointfilename
character c80tmp*80,c2000tmp*2000
integer imode
logical filealive
integer selatm(ncenter)

if (imode==2) then
	igridsel=9
else
	ntotlow=125000
	ntotmed=512000
	ntothigh=1728000
	do while(.true.)
		write(*,*)
		write(*,*) "Please select a method to set up grid"
		write(*,"(a,f7.3,a)") " -10 Set extension distance of grid range for mode 1~4, current:",aug3D," Bohr"
		write(*,*) "1 Low quality grid,    covering whole system, about 125000 points in total"
		write(*,*) "2 Medium quality grid, covering whole system, about 512000 points in total"
		write(*,*) "3 High quality grid,   covering whole system, about 1728000 points in total"
		write(*,*) "4 Input the number of points or grid spacing in X,Y,Z, covering whole system"
		write(*,*) "5 Input original point, grid spacings, and the number of points"
		write(*,*) "6 Input center coordinate, number of points and extension distance"
		write(*,*) "7 The same as 6, but input two atoms, the midpoint will be defined as center"
		write(*,*) "8 Use grid setting of another cube file"
		if (ifPBC>0) write(*,"(a)") " 9 Use translation vectors of current cell, manually specify origin, box lengths and grid spacing"
		write(*,*) "10 Set box of grid data visually using a GUI window"
		write(*,*) "11 Select a batch of atoms, set extension distance and grid spacing"
		if (imode==1) write(*,*) "100 Load a set of points from external file"
		read(*,*) igridsel
		if (igridsel/=-10) exit
		write(*,*) "Input extension distance in Bohr, e.g. 6.5"
		read(*,*) aug3D
	end do
end if

if (igridsel==100) then !Load points rather than set up grid
	write(*,*) "Input the path of the file containing points, e.g. C:\ltwd.txt"
	write(*,*) "Note: See program manual for the format of the file"
	do while(.true.)
		read(*,"(a)") pointfilename
		inquire(file=pointfilename,exist=filealive)
		if (filealive) then
			open(10,file=pointfilename,status="old")
			read(10,*) numextpt
			write(*,"(a,i10,a)") ' There are',numextpt,' points'
			if (allocated(extpt)) deallocate(extpt)
			allocate(extpt(numextpt,4))
			do itmp=1,numextpt
				read(10,*) extpt(itmp,1:3)
			end do
			close(10)
			exit
		else
			write(*,*) "Error: File cannot be found, input again"
		end if
	end do
	write(*,*) "Please wait..."
else
	molxlen=(maxval(a%x)-minval(a%x))+2*aug3D
	molylen=(maxval(a%y)-minval(a%y))+2*aug3D
	molzlen=(maxval(a%z)-minval(a%z))+2*aug3D
	if (molxlen==0.0D0.or.molylen==0.0D0.or.molzlen==0.0D0) then !Avoid catastrophe when aug3D=0 and system is plane
		write(*,"(a,/)") " WARNING: The box size in one of Cartesian axis is zero, &
		the calculation cannot be proceeded. Therefore, the size of corresponding direction is automatically set to 3 Bohr"
		if (molxlen==0D0) then
			molxlen=3D0
		else if (molylen==0D0) then
			molylen=3D0
		else if (molzlen==0D0) then
			molzlen=3D0
		end if
	end if
	if (igridsel==1.or.igridsel==2.or.igridsel==3) then
		if (igridsel==1) dx=(molxlen*molylen*molzlen/dfloat(ntotlow))**(1.0D0/3.0D0)
		if (igridsel==2) dx=(molxlen*molylen*molzlen/dfloat(ntotmed))**(1.0D0/3.0D0)
		if (igridsel==3) dx=(molxlen*molylen*molzlen/dfloat(ntothigh))**(1.0D0/3.0D0)
		dy=dx
		dz=dx
		nx=nint(molxlen/dx)+1
		ny=nint(molylen/dy)+1
		nz=nint(molzlen/dz)+1
		orgx=minval(a%x)-aug3D
		orgy=minval(a%y)-aug3D
		orgz=minval(a%z)-aug3D
	else if (igridsel==4) then
		write(*,*) "Input the number of grid points in X,Y,Z directions, e.g. 139,59,80"
		write(*,"(a)") " or input grid spacing (Bohr) in X,Y,Z directions, e.g. 0.05,0.08,0.08  (if only input one value, it will be applied to all directions)"
		read(*,"(a)") c80tmp
		if (index(c80tmp,'.')/=0) then
			if (index(c80tmp,',')/=0) then
				read(c80tmp,*) dx,dy,dz
			else
				read(c80tmp,*) tmp
				dx=tmp
				dy=tmp
				dz=tmp
			end if
			nx=molxlen/dx+1
			ny=molylen/dy+1
			nz=molzlen/dz+1
		else
			read(c80tmp,*) nx,ny,nz
			dx=molxlen/(nx-1)
			dy=molylen/(ny-1)
			dz=molzlen/(nz-1)
		end if
		orgx=minval(a%x)-aug3D
		orgy=minval(a%y)-aug3D
		orgz=minval(a%z)-aug3D
	else if (igridsel==5) then
		write(*,*) "Input X,Y,Z coordinate of original point (Bohr), e.g. 0.1,4,-1"
		read(*,*) orgx,orgy,orgz
		write(*,*) "Input grid spacings in X,Y,Z directions (Bohr), e.g. 0.1,0.1,0.15"
		read(*,*) dx,dy,dz
		write(*,*) "Input the number of points in X,Y,Z directions, e.g. 139,59,80"
		read(*,*) nx,ny,nz
	else if (igridsel==6.or.igridsel==7) then
		if (igridsel==6) then
			write(*,*) "Input X,Y,Z coordinate of center (Angstrom)"
			read(*,*) cenx,ceny,cenz
			cenx=cenx/b2a
			ceny=ceny/b2a
			cenz=cenz/b2a
		else if (igridsel==7) then
			write(*,*) "Input index of the two atoms, e.g. 2,5"
			write(*,*) "If the two indices are identical, box center will be placed at the nucleus"
			read(*,*) indatm1,indatm2
			cenx=(a(indatm1)%x+a(indatm2)%x)/2.0D0
			ceny=(a(indatm1)%y+a(indatm2)%y)/2.0D0
			cenz=(a(indatm1)%z+a(indatm2)%z)/2.0D0
		end if
		write(*,*) "Input the number of points in X,Y,Z directions, e.g. 40,40,25"
		read(*,*) nx,ny,nz
		write(*,*) "Input the extended distance in X,Y,Z directions (Bohr), e.g. 4.0,4.0,6.5"
		read(*,*) aug3Dx,aug3Dy,aug3Dz
		orgx=cenx-aug3Dx
		orgy=ceny-aug3Dy
		orgz=cenz-aug3Dz
		dx=aug3Dx*2D0/(nx-1)
		dy=aug3Dy*2D0/(ny-1)
		dz=aug3Dz*2D0/(nz-1)
	else if (igridsel==8) then
		write(*,*) "Input path of a cube file, e.g. C:\wake_up_girls.cub"
		do while(.true.)
			read(*,"(a)") cubefilename
			inquire(file=cubefilename,exist=filealive)
			if (filealive) then
				open(10,file=cubefilename,status="old")
				read(10,*)
				read(10,*)
				read(10,*) nouse,orgx,orgy,orgz
				read(10,*) nx,gridv1
				read(10,*) ny,gridv2
				read(10,*) nz,gridv3
				close(10)
                dx=gridv1(1);dy=gridv2(2);dz=gridv3(3)
				exit
			else
				write(*,*) "Error: File cannot be found, input again"
			end if
		end do
	else if (igridsel==9) then
		call setgrid_for_PBC
	else if (igridsel==10) then
		call setboxGUI
	else if (igridsel==11) then
        write(*,*) "Input index of the atoms to define a fragment, e.g. 2,3,7-10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nselatm,selatm)
        write(*,*) "Input extension distance around the fragment in Bohr, e.g. 7.5"
        write(*,*) "To input in Angstrom, add ""A"" suffix, e.g. 3.8 A"
        read(*,"(a)") c80tmp
        read(c80tmp,*) extdist
        if (index(c80tmp,'A')/=0) extdist=extdist/b2a
        write(*,*) "Input grid spacing in Bohr, e.g. 0.2"
        read(*,*) dx
        dy=dx;dz=dx
        orgx=minval(a(selatm(1:nselatm))%x)-extdist
        orgy=minval(a(selatm(1:nselatm))%y)-extdist
        orgz=minval(a(selatm(1:nselatm))%z)-extdist
        endx=maxval(a(selatm(1:nselatm))%x)+extdist
        endy=maxval(a(selatm(1:nselatm))%y)+extdist
        endz=maxval(a(selatm(1:nselatm))%z)+extdist
        nx=nint((endx-orgx)/dx)
        ny=nint((endy-orgy)/dy)
        nz=nint((endz-orgz)/dz)
	end if
    
    if (igridsel/=9) then
        gridv1=0;gridv1(1)=dx
        gridv2=0;gridv2(2)=dy
        gridv3=0;gridv3(3)=dz
    else
        dx=gridv1(1)
        dy=gridv2(2)
        dz=gridv3(3)
    end if
    call getgridend
	write(*,"(' Coordinate of origin in X,Y,Z is   ',3f12.6,' Bohr')") orgx,orgy,orgz
	write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6,' Bohr')") endx,endy,endz
    if (igridsel/=9) then
	    write(*,"(' Grid spacing in X,Y,Z is',3f12.6,' Bohr')") dx,dy,dz
	    write(*,"(' Number of points in X,Y,Z is',3i5,'   Total:',i12)") nx,ny,nz,nx*ny*nz
    else
	    write(*,"(' Grid vector 1 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv1,dsqrt(sum(gridv1**2))
	    write(*,"(' Grid vector 2 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv2,dsqrt(sum(gridv2**2))
	    write(*,"(' Grid vector 3 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv3,dsqrt(sum(gridv3**2))
	    write(*,"(' Number of points in three directions is',3i5,'  Total:',i12)") nx,ny,nz,nx*ny*nz
    end if
end if
end subroutine


!!------- A subroutine directly define grid for PBC case, embedded by subroutine setgrid, setgridfixspc, etc.
!Note that after employing the following rule, then in e.g. X direction,
!the first grid not only differs from last grid by translation vector, but also by a grid spacing. See Section 3.6 for illustration
subroutine setgrid_for_PBC
use defvar
implicit real*8 (a-h,o-z)
character c80tmp*80
write(*,*) "Now input X,Y,Z of origin in Bohr, e.g. 0.2,0,-5.5"
write(*,*) "You can also input in Angstrom by adding ""A"" suffix, e.g. 0.2,0,-5.5 A"
write(*,*) "If press ENTER button directly, (0,0,0) will be used"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    orgx=0;orgy=0;orgz=0
else
    read(c80tmp,*) orgx,orgy,orgz
    if (index(c80tmp,'A')/=0) then
        orgx=orgx/b2a
        orgy=orgy/b2a
        orgz=orgz/b2a
    end if
end if
write(*,*) "Now input lengths of three dimensions of the box in Bohr, e.g. 8.7,9.1,6.55"
write(*,*) "You can also input in Angstrom by adding ""A"" suffix, e.g. 8.7,9.1,6.55 A"
write(*,"(a)") " If length of a dimension is set to be 0, then box length of that dimension will be equal to cell length"
write(*,*) "Pressing ENTER button directly corresponds to inputting 0,0,0"
read(*,"(a)") c80tmp
v1len=dsqrt(sum(cellv1**2))
v2len=dsqrt(sum(cellv2**2))
v3len=dsqrt(sum(cellv3**2))
if (c80tmp/=" ") then
    read(c80tmp,*) atmp,btmp,ctmp
    if (index(c80tmp,'A')/=0) then
        atmp=atmp/b2a
        btmp=btmp/b2a
        ctmp=ctmp/b2a
    end if
    if (atmp/=0) v1len=atmp
    if (btmp/=0) v2len=btmp
    if (ctmp/=0) v3len=ctmp
end if
write(*,*) "Now input grid spacing in Bohr, e.g. 0.25"
write(*,*) "Hint about grid quality: 0.2 is fine, 0.3 is coarse, 0.4 is very poor"
write(*,"(a)") " Note: The grid spacing will be automatically slightly altered so that number of grids in each direction is integer"
read(*,*) grdspc
nx=nint(v1len/grdspc)
grdspcv1=v1len/nx
gridv1(:)=cellv1(:)/dsqrt(sum(cellv1**2)) * grdspcv1
ny=nint(v2len/grdspc)
grdspcv2=v2len/ny
gridv2(:)=cellv2(:)/dsqrt(sum(cellv2**2)) * grdspcv2
nz=nint(v3len/grdspc)
grdspcv3=v3len/nz
gridv3(:)=cellv3(:)/dsqrt(sum(cellv3**2)) * grdspcv3
end subroutine



!!---- Set up grid setting with fixed grid spacing
!Similar to setgridforbasin, but not so specific and complicated, thus may be used for other subroutines
subroutine setgridfixspc
use defvar
use GUI
implicit real*8 (a-h,o-z)
real*8 molxlen,molylen,molzlen
real*8 :: spclowqual=0.2D0,spcmedqual=0.1D0,spchighqual=0.06D0,spclunaqual=0.04D0
character c80tmp*80,cubefilename*200
do while(.true.)
	orgx=minval(a%x)-aug3D
	orgy=minval(a%y)-aug3D
	orgz=minval(a%z)-aug3D
	endx=maxval(a%x)+aug3D
	endy=maxval(a%y)+aug3D
	endz=maxval(a%z)+aug3D
	molxlen=endx-orgx
	molylen=endy-orgy
	molzlen=endz-orgz
	ntotlow=(nint(molxlen/spclowqual)+1)*(nint(molylen/spclowqual)+1)*(nint(molzlen/spclowqual)+1)
	ntotmed=(nint(molxlen/spcmedqual)+1)*(nint(molylen/spcmedqual)+1)*(nint(molzlen/spcmedqual)+1)
	ntothigh=(nint(molxlen/spchighqual)+1)*(nint(molylen/spchighqual)+1)*(nint(molzlen/spchighqual)+1)
	ntotluna=(nint(molxlen/spclunaqual)+1)*(nint(molylen/spclunaqual)+1)*(nint(molzlen/spclunaqual)+1)
	
	write(*,*) "Please select a method for setting up grid"
	write(*,"(a,f10.5,a)") " -10 Set grid extension distance for mode 1~6, current:",aug3D," Bohr"
	write(*,"(a,f4.2,a,i14)") " 1 Low quality grid, spacing=",spclowqual," Bohr, number of grids:    ",ntotlow
	write(*,"(a,f4.2,a,i14)") " 2 Medium quality grid, spacing=",spcmedqual," Bohr, number of grids: ",ntotmed
	write(*,"(a,f4.2,a,i14)") " 3 High quality grid, spacing=",spchighqual," Bohr, number of grids:   ",ntothigh
	write(*,"(a,f4.2,a,i14)") " 4 Lunatic quality grid, spacing=",spclunaqual," Bohr, number of grids:",ntotluna
	write(*,*) "5 Only input grid spacing, automatically set other parameters"
	write(*,*) "6 Only input the number of points in X,Y,Z, automatically set other parameters"
	write(*,*) "7 Input original point, grid spacings, and number of points"
	write(*,*) "8 Set center position, grid spacing and box length"
	write(*,*) "9 Use grid setting of another cube file"
	write(*,*) "10 Set box of grid data visually using a GUI window"
    if (ifPBC>0) write(*,"(a)") " 11 Use translation vectors of current cell, manually specify origin, box lengths and grid spacing"
	read(*,*) igridsel
	if (igridsel/=-10) then
		exit
	else
		write(*,*) "Input extension distance (Bohr), e.g. 6.5"
		read(*,*) aug3D
	end if
end do

!Note: orgx,orgy,orgz,endx,endy,endz as well as molx/y/zlen for igridsel==1~6 have already been set above
if (igridsel==1.or.igridsel==2.or.igridsel==3.or.igridsel==4.or.igridsel==5) then
	if (igridsel==1) dx=spclowqual
	if (igridsel==2) dx=spcmedqual
	if (igridsel==3) dx=spchighqual
	if (igridsel==4) dx=spclunaqual
	if (igridsel==5) then
		write(*,*) "Input the grid spacing (Bohr), e.g. 0.08"
		read(*,*) dx
	end if
	dy=dx
	dz=dx
	nx=nint(molxlen/dx)+1
	ny=nint(molylen/dy)+1
	nz=nint(molzlen/dz)+1
else if (igridsel==6) then
	write(*,*) "Input the number of grid points in X,Y,Z direction, e.g. 139,59,80"
	read(*,*) nx,ny,nz
	dx=molxlen/(nx-1)
	dy=molylen/(ny-1)
	dz=molzlen/(nz-1)
else if (igridsel==7) then
	write(*,*) "Input X,Y,Z coordinate of original point (Bohr), e.g. 0.1,4,-1"
	read(*,*) orgx,orgy,orgz
	write(*,*) "Input grid spacings in X,Y,Z directions (Bohr), e.g. 0.1,0.1,0.15"
	read(*,*) dx,dy,dz
	write(*,*) "Input the number of points in X,Y,Z directions, e.g. 139,59,80"
	read(*,*) nx,ny,nz
else if (igridsel==8) then
	write(*,*) "Input X,Y,Z coordinate of box center (in Angstrom)"
	write(*,*) "or input such as a8 to take the coordinate of atom 8 as box center"
	write(*,*) "or input such as a3,a7 to take the midpoint of atom 3 and atom 7 as box center"
	read(*,"(a)") c80tmp
	if (c80tmp(1:1)=='a') then
		do ich=1,len_trim(c80tmp)
			if (c80tmp(ich:ich)==',') exit
		end do
		if (ich==len_trim(c80tmp)+1) then
			read(c80tmp(2:),*) itmp
			cenx=a(itmp)%x
			ceny=a(itmp)%y
			cenz=a(itmp)%z
		else
			read(c80tmp(2:ich-1),*) itmp
			read(c80tmp(ich+2:),*) jtmp			
			cenx=(a(itmp)%x+a(jtmp)%x)/2D0
			ceny=(a(itmp)%y+a(jtmp)%y)/2D0
			cenz=(a(itmp)%z+a(jtmp)%z)/2D0
		end if
	else
		read(c80tmp,*) cenx,ceny,cenz
		cenx=cenx/b2a
		ceny=ceny/b2a
		cenz=cenz/b2a
	end if
	write(*,*) "Input the grid spacing (Bohr), e.g. 0.08"
	read(*,*) dx
	dy=dx
	dz=dx
	write(*,*) "Input the box lengths in X,Y,Z direction (Bohr), e.g. 8.0,8.0,13.5"
	read(*,*) molxlen,molylen,molzlen
	orgx=cenx-molxlen/2D0
	orgy=ceny-molylen/2D0
	orgz=cenz-molzlen/2D0
	nx=nint(molxlen/dx)+1
	ny=nint(molylen/dy)+1
	nz=nint(molzlen/dz)+1
else if (igridsel==9) then
	write(*,*) "Input path of a cube file, e.g. C:\oppai.cub"
	do while(.true.)
		read(*,"(a)") cubefilename
		inquire(file=cubefilename,exist=alive)
		if (alive) then
			open(10,file=cubefilename,status="old")
			read(10,*)
			read(10,*)
			read(10,*) nouse,orgx,orgy,orgz
			read(10,*) nx,dx
			read(10,*) ny,rnouse,dy
			read(10,*) nz,rnouse,rnouse,dz
			close(10)
			exit
		else
			write(*,*) "Error: File cannot be found, input again"
		end if
	end do
else if (igridsel==11) then
	call setgrid_for_PBC
end if

if (igridsel==10) call setboxGUI

if (igridsel/=11) then
    gridv1=0;gridv1(1)=dx
    gridv2=0;gridv2(2)=dy
    gridv3=0;gridv3(3)=dz
else
    dx=gridv1(1)
    dy=gridv2(2)
    dz=gridv3(3)
end if
call getgridend !Generate endx,endy,endz
write(*,"(' Coordinate of origin in X,Y,Z is   ',3f12.6,' Bohr')") orgx,orgy,orgz
write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6,' Bohr')") endx,endy,endz
if (igridsel/=11) then
	write(*,"(' Grid spacing in X,Y,Z is',3f12.6,' Bohr')") dx,dy,dz
	write(*,"(' Number of points in X,Y,Z is',3i5,'   Total:',i12)") nx,ny,nz,nx*ny*nz
else
	write(*,"(' Grid vector 1 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv1,dsqrt(sum(gridv1**2))
	write(*,"(' Grid vector 2 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv2,dsqrt(sum(gridv2**2))
	write(*,"(' Grid vector 3 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv3,dsqrt(sum(gridv3**2))
	write(*,"(' Number of points in three directions is',3i5,'  Total:',i12)") nx,ny,nz,nx*ny*nz
end if
end subroutine



!!----- Return differential volume according to grid translation vector
subroutine calc_dvol(dvol)
use defvar
use util
real*8 dvol,mat(3,3)
mat(:,1)=gridv1(:)
mat(:,2)=gridv2(:)
mat(:,3)=gridv3(:)
dvol=abs(detmat(mat))
end subroutine



!!----- Get XYZ coordinate of grid based on grid indiex
subroutine getgridxyz(i,j,k,tmpx,tmpy,tmpz)
use defvar
integer i,j,k
real*8 tmpx,tmpy,tmpz
tmpx=orgx+gridv1(1)*(i-1)+gridv2(1)*(j-1)+gridv3(1)*(k-1)
tmpy=orgy+gridv1(2)*(i-1)+gridv2(2)*(j-1)+gridv3(2)*(k-1)
tmpz=orgz+gridv1(3)*(i-1)+gridv2(3)*(j-1)+gridv3(3)*(k-1)
end subroutine



!!----- Get X coordinate of grid based on grid index
real*8 function getgridx(i,j,k)
use defvar
integer i,j,k
getgridx=orgx+gridv1(1)*(i-1)+gridv2(1)*(j-1)+gridv3(1)*(k-1)
end function
!!----- Get Y coordinate of grid based on grid index
real*8 function getgridy(i,j,k)
use defvar
integer i,j,k
getgridy=orgy+gridv1(2)*(i-1)+gridv2(2)*(j-1)+gridv3(2)*(k-1)
end function
!!----- Get Z coordinate of grid based on grid index
real*8 function getgridz(i,j,k)
use defvar
integer i,j,k
getgridz=orgz+gridv1(3)*(i-1)+gridv2(3)*(j-1)+gridv3(3)*(k-1)
end function



!!----- Get x/y/z of ending point of grid data. endx,endy,endz are global variables
subroutine getgridend
use defvar
endx=orgx+gridv1(1)*(nx-1)+gridv2(1)*(ny-1)+gridv3(1)*(nz-1)
endy=orgy+gridv1(2)*(nx-1)+gridv2(2)*(ny-1)+gridv3(2)*(nz-1)
endz=orgz+gridv1(3)*(nx-1)+gridv2(3)*(ny-1)+gridv3(3)*(nz-1)
end subroutine




!----- Return coordinate of specific vertex of the grid data
!idx is the index of following map
!    5-------8
!   /|      /|
!  6-+-----7 |
!  | |     | |
!  | 1-----+-4
!  |/      |/
!  2-------3
!
!   Z
!   |
!   0---Y    
!  / 
! X
subroutine gridvertex(idx,vertx,verty,vertz)
use defvar
integer idx
real*8 vertx,verty,vertz
if (idx==1) then
    icell=0;jcell=0;kcell=0
else if (idx==2) then
    icell=1;jcell=0;kcell=0
else if (idx==3) then
    icell=1;jcell=1;kcell=0
else if (idx==4) then
    icell=0;jcell=1;kcell=0
else if (idx==5) then
    icell=0;jcell=0;kcell=1
else if (idx==6) then
    icell=1;jcell=0;kcell=1
else if (idx==7) then
    icell=1;jcell=1;kcell=1
else if (idx==8) then
    icell=0;jcell=1;kcell=1
end if
vertx=orgx+icell*(nx-1)*gridv1(1)+jcell*(ny-1)*gridv2(1)+kcell*(nz-1)*gridv3(1)
verty=orgy+icell*(nx-1)*gridv1(2)+jcell*(ny-1)*gridv2(2)+kcell*(nz-1)*gridv3(2)
vertz=orgz+icell*(nx-1)*gridv1(3)+jcell*(ny-1)*gridv2(3)+kcell*(nz-1)*gridv3(3)
end subroutine

!!----- Wrapper of gridvertex, input two vertex indices return two sets of coordinates
subroutine gridvertex2(idx,jdx,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
integer idx,jdx
real*8 vert1x,vert1y,vert1z,vert2x,vert2y,vert2z
call gridvertex(idx,vert1x,vert1y,vert1z)
call gridvertex(jdx,vert2x,vert2y,vert2z)
end subroutine


!!------- Return max x,y,z of all grid vertices
subroutine gridmaxxyz(xmax,ymax,zmax)
use defvar
implicit real*8 (a-h,o-z)
real*8 xmax,ymax,zmax
xmax=-999
ymax=-999
zmax=-999
do idx=1,8
    call gridvertex(idx,vertx,verty,vertz)
    if (vertx>xmax) xmax=vertx
    if (verty>ymax) ymax=verty
    if (vertz>zmax) zmax=vertz
end do
end subroutine

!!------- Return mix x,y,z of all grid vertices
subroutine gridminxyz(xmin,ymin,zmin)
use defvar
implicit real*8 (a-h,o-z)
real*8 xmin,ymin,zmin
xmin=999
ymin=999
zmin=999
do idx=1,8
    call gridvertex(idx,vertx,verty,vertz)
    if (vertx<xmin) xmin=vertx
    if (verty<ymin) ymin=verty
    if (vertz<zmin) zmin=vertz
end do
end subroutine



!!!------ A concise routine specifically for filling up electron density to "rhocub" array
subroutine saverhocub
use defvar
use function
implicit real*8(a-h,o-z)
if (allocated(rhocub)) then
    if (size(rhocub,1)==nx.and.size(rhocub,2)==ny.and.size(rhocub,3)==nz) return !Do need to calculate again
else
    allocate(rhocub(nx,ny,nz))
end if
write(*,*) "Calculating grid data of electron density..."
ifinish=0
!$OMP PARALLEL DO SHARED(rhocub,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,tmprho) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			rhocub(i,j,k)=fdens(tmpx,tmpy,tmpz)
		end do
	end do
	ifinish=ifinish+1
    call showprog(ifinish,nz)
end do
!$OMP END PARALLEL DO
end subroutine


!!----------- If three vectors of the grid data are parallel to X/Y/Z axis, respectively. Return 1 means yes, 0 means no
integer function ifgridortho
use defvar
if (abs(gridv1(2))<1E-10.and.abs(gridv1(3))<1E-10.and.abs(gridv2(1))<1E-10.and.abs(gridv2(3))<1E-10.and.abs(gridv3(1))<1E-10.and.abs(gridv3(2))<1E-10) then
    ifgridortho=1
else
    ifgridortho=0
end if
end function


!!----------- Show grid data information
!itype=1: Grid information
!itype=2: Statistical information
subroutine showgridinfo(itype)
use defvar
implicit real*8 (a-h,o-z)
integer itype
type(content) maxv,minv
if (itype==1) then
	write(*,*)
	write(*,*) "Grid information:"
	write(*,"(' Translation vector:        X           Y           Z     (Bohr)')")
	write(*,"(a20,3f12.6)") "Vector 1: ",gridv1
	write(*,"(a20,3f12.6)") "Vector 2: ",gridv2
	write(*,"(a20,3f12.6)") "Vector 3: ",gridv3
	call calc_dvol(dvol)
    write(*,"(a,f12.6,' Bohr^3')") " Volume of each grid:",dvol
	write(*,"(' The range of x is from ',f12.6,' to ',f12.6,' Bohr,' i5,' points')") ,orgx,endx,nx
	write(*,"(' The range of y is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgy,endy,ny
	write(*,"(' The range of z is from ',f12.6,' to ',f12.6,' Bohr,',i5,' points')") ,orgz,endz,nz
	write(*,"(' Total number of grid points:',i12)") nx*ny*nz
	write(*,"(' This grid data will take up at least',i6,' MB memory')") nx*ny*nz*8/1024/1024
else if (itype==2) then
    maxv%value=cubmat(1,1,1)
    maxv%x=orgx
    maxv%y=orgy
    maxv%z=orgz
    minv%value=cubmat(1,1,1)
    minv%x=orgx
    minv%y=orgy
    minv%z=orgz
    sumuppos=0D0
    sumupneg=0D0
    do k=1,nz
	    do j=1,ny
		    do i=1,nx
			    if (cubmat(i,j,k)>0) sumuppos=sumuppos+cubmat(i,j,k)
			    if (cubmat(i,j,k)<0) sumupneg=sumupneg+cubmat(i,j,k)
			    if (cubmat(i,j,k)>maxv%value) then
				    maxv%value=cubmat(i,j,k)
                    call getgridxyz(i,j,k,maxv%x,maxv%y,maxv%z)
			    end if
			    if (cubmat(i,j,k)<minv%value) then
				    minv%value=cubmat(i,j,k)
                    call getgridxyz(i,j,k,minv%x,minv%y,minv%z)
			    end if
		    end do
	    end do
    end do

	call calc_dvol(dvol)
	write(*,*)
	write(*,*) "Statistical information:"
	write(*,"(' Global minimum value:',E16.8,' at',3f10.5,' Bohr')") minv%value,minv%x,minv%y,minv%z
	write(*,"(' Global maximum value:',E16.8,' at',3f10.5,' Bohr')") maxv%value,maxv%x,maxv%y,maxv%z
	write(*,"(' Differential element:',f15.10,' Bohr^3')") dvol
	write(*,"(' Summing up positive value in grid file:  ',f30.10)") sumuppos
	write(*,"(' After multiplied by differential element:',f30.10)") sumuppos*dvol
	write(*,"(' Summing up negative value in grid file:  ',f30.10)") sumupneg
	write(*,"(' After multiplied by differential element:',f30.10)") sumupneg*dvol
	write(*,"(' Summing up all value in grid file:       ',f30.10)") sumuppos+sumupneg
	write(*,"(' After multiplied by differential element:',f30.10)") (sumuppos+sumupneg)*dvol

    !if (sum(gridv1*gridv2*gridv3)/=0) then
    !    write(*,*)
    !    write(*,*) "   Warning! Warning! Warning! Warning! Warning! Warning! Warning! Warning!"
    !    write(*,"(a)") " The grid is not rectangle, in this case, many outputs (including the ones shown above) and &
    !    analysis results are completely meaningless. &
    !    Only the grid data calculation function in main function 13 must work normally"
    !    write(*,*) "Press ENTER button to continue"
    !    read(*,*)
    !end if
end if
end subroutine