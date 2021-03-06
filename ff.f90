subroutine fast(psi,n,del,psi_omega)

implicit real*8 (a-z)
integer*8 :: n,k,l
complex*16, dimension (n) :: psi(n)
complex*16 :: i,y
complex*16, dimension (-n/2:n/2) ::  psi_omega
external fft

call fft(psi,n,1)

do k = n/2+1,n
   f = dble(k-1-n)/dble(n)/del
   psi_omega(k-1-n)=psi(k)
enddo

do k = 1,n/2+1
   f = dble(k-1)/dble(n)/del
   psi_omega(k-1)=psi(k)
enddo

return
end
!
! Subroutine that handles the factor of n that appears in the FFT
! routine of
! Numerical Recipes (see below).



subroutine fft(data,n,sign)
implicit real*8 (a-z)

integer*8 :: i,n,sign
complex*16, dimension (n) :: data
external fourl

call four1(data,n,sign)
fact = 1.0d0
if (sign.eq.-1) fact = 1.0d0/n

do i = 1,n
   data(i) = data(i)*fact
enddo

return
end


subroutine four1(data,nn,isign)
implicit real*8 (a-z)

integer*8 :: isign, nn
real*8, dimension (2*nn) ::  data
integer*8 :: i, istep, j, m, mmax, n
real*8 :: tempi, tempr, theta, wi, wpi, wpr, wr, wtemp

n = 2*nn
j = 1
do i = 1,n,2

   if(j .gt. i)then
      tempr = data(j)
      tempi = data(j+1)
      data(j) = data(i)
      data(j+1) = data(i+1)
      data(i) = tempr
      data(i+1) = tempi
   endif

   m = n/2
1          if((m .ge. 2).and.(j .gt. m))then
      j = j - m
      m = m/2
      goto 1
   endif

   j = j + m
enddo

mmax = 2
2       if(n .gt. mmax)then
   istep = 2*mmax
   theta = 6.28318530717959d0/(isign*mmax)
   wpr = -2.d0*dsin(0.5d0*theta)**2
   wpi = dsin(theta)
   wr = 1.0d0
   wi = 0.0d0

   do m = 1,mmax,2
       do i = m,n,istep
         j = i + mmax
         tempr = wr*data(j) - wi*data(j+1)
         tempi = wr*data(j+1) + wi*data(j)
         data(j) = data(i) - tempr
         data(j+1) = data(i+1) - tempi
         data(i) = data(i) + tempr
         data(i+1) = data(i+1) + tempi
      enddo
      wtemp = wr
      wr = wr*wpr - wi*wpi + wr
      wi = wi*wpr + wtemp*wpi + wi
   enddo
   mmax = istep
   goto 2
endif
return
end
