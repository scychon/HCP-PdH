program anal_hmap

        implicit none
        integer :: idx,i,j,k
        integer :: ix,iy
        integer :: dnmax = 200
        real*8  :: drsq, drmaxsq, twosig
        real*8 :: rScale = 0.3479
        integer :: ngridvec = 230, nslice = 31
        integer,dimension(230,230,31) :: num_h_grid ! number of H atoms on this grid point
        real*8, dimension(230,230,31) :: rhoh_sum   ! output htrajectory

        write(*,*) 'start reading hmap.dat'
        open(unit=21,file='hmap.dat', status = 'old')
        do i=1,ngridvec
          do j=1,ngridvec
            read(21,*) ix,iy,num_h_grid(i,j,:)
          enddo
          read(21,*)
        enddo
        close(21)
        write(*,*) 'finished to read hmap.dat'
 
        drmaxsq = dnmax**2 +.1
        twosig = 2.*10./rScale
        write(*,*) 'maxval ', maxval(num_h_grid)
        write(*,*) 'maxval ', minval(num_h_grid)
        open(unit=21,file='error.log')

        rhoh_sum = 0.
        !$OMP PARALLEL &
        !$OMP   PRIVATE (idx,ix,iy,i,j,k,drsq)
        !$OMP DO
        do idx=1,ngridvec*ngridvec
            ix=(idx-1)/ngridvec +1
            iy=idx-(ix-1)*ngridvec
            do i=max(1,ix-dnmax),min(ngridvec,ix+dnmax)
              do j=max(1,iy-dnmax),min(ngridvec,iy+dnmax)
                do k=1,nslice
                  if (num_h_grid(i,j,k) > 0) then
!                    write(*,*) ix,iy,i,j,'test write'
                    drsq = (i-ix)**2+(j-iy)**2
                    if (drsq < drmaxsq) then
!                        write(*,*) i,j,k,ix,iy,'test write'
                        if ((ix<1) .or. (iy<1)) then
                        write(21,*) i,j,k,ix,iy,num_h_grid(i,j,k),'test write'
                        endif
                        rhoh_sum(ix,iy,k) = rhoh_sum(ix,iy,k)+exp(-drsq/(twosig))*num_h_grid(i,j,k)
                    endif
                  endif
                enddo
              enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        open(unit=22,file='rhohsum.dat')
        do i=1,ngridvec
          do j=1,ngridvec
            write(22,*) i,j,rhoh_sum(i,j,:)
          enddo
          write(22,*) ''
        enddo
        close(22)

end program anal_hmap
