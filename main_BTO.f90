program natddft
      use convert
    implicit real*8 (a-z)
    include 'header_BTO.f90'
   !  include 'declare_varaible.f90'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    call cpu_time(start_time)
    ii=dcmplx(0.0d0,1.0d0)
    pi=4.0d0*datan(1.0d0)


!CCCCCCCCCCCCCCCCCCCCCCCCC calculate P^lm CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!        gap=2.6d0
       gap=1.49d0

!        open(1,File='./DFT/kpoints.dat')
!        do k=1,n_k
!           read(1,*) k_m(1,k), k_m(2,k), k_m(3,k)
!        end do
!        close(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(1,File='data_BTO/energy_all')
    do k=1,n_k_full
       do i=1,n_bands
       do j=1,n_spin

       read(1,*) energy_full(k,j,i)
  !   write(*,*) i
   end do
        end do
     end do
     close(1)

!     write(*,*) energy_full(1,1,21)-energy_full(1,1,20)
!     stop
!       do i=1,n_k_full
!          do i1=1,n_bands
!          energy_full(i,2,i1)=energy_full(i,1,i1)
!          write(*,*) i,i1,energy_full(i,2,i1)
!          end do
!       end do



     delta_k_x=5.537065703617170E-002-0.0001d0
     delta_k_y=5.537065703617170E-002-0.0001d0
!      delta_k_y=0.0661d0-0.0001d0
     delta_k_z=5.258371988240430E-002-0.0001d0


     open(154,File='numbers')
     do i=1,n_k_full
        i1=k_m_full(1,i)/delta_k_x
        i2=k_m_full(2,i)/delta_k_y
        i3=k_m_full(3,i)/delta_k_z

        do j=1,n_bands
            do m=1,n_spin
               energy_full_3D(i1,i2,i3,m,j)=energy_full(i,m,j)
               write(154,*) i1,i2,i3,energy_full_3D(i1,i2,i3,m,j)
            end do
         end do
     end do
     close(154)


!       do i1=-n_k_1_full_max,n_k_1_full_max
!          do i2=-n_k_2_full_max,n_k_2_full_max
!             do i3=-n_k_3_full_max,n_k_3_full_max
!
!             i1p=i1
!             if(i1.eq.n_k_1_full_max) i1p=-n_k_1_full_max
!
!             i2p=i2
!             if(i2.eq.n_k_2_full_max) i2p=-n_k_2_full_max
!
!             i3p=i3
!             if(i3.eq.n_k_3_full_max) i3p=-n_k_3_full_max
!
!
!          do j=1,n_bands
!              do m=1,n_spin
!                 energy_full_3D(i1,i2,i3,m,j)=energy_full_3D(i1p,i2p,i3p,1,j)
!                 write(*,*)  energy_full_3D(i1,i2,i3,m,j)
!              end do
!           end do
!
!             end do
!          end do
!       end do

     do j=n_c,n_c_max
         do j_spin=1,n_spin
            do i1=-n_k_1_full_max,n_k_1_full_max-1
               do i2=-n_k_2_full_max,n_k_2_full_max-1
                  do i3=-n_k_3_full_max,n_k_3_full_max-1

            grad_vel_c_full(i1,i2,i3,j_spin,j,1)=(a_lat/(2.0d0*pi))*(energy_full_3D(i1+1,i2,i3,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_x
            grad_vel_c_full(i1,i2,i3,j_spin,j,2)=(a_lat/(2.0d0*pi))*(energy_full_3D(i1,i2+1,i3,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_y
            grad_vel_c_full(i1,i2,i3,j_spin,j,3)=(a_lat/(2.0d0*pi))*(energy_full_3D(i1,i2,i3+1,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_z

                   end do
                end do
             end do
          end do
      end do


     open(1,File='data_BTO/velocity_v')
     do j=n_v_min,n_v
         do j_spin=1,n_spin
            do i1=-n_k_1_full_max,n_k_1_full_max-1
               do i2=-n_k_2_full_max,n_k_2_full_max-1
                  do i3=-n_k_3_full_max,n_k_3_full_max-1

            grad_vel_v_full(i1,i2,i3,j_spin,j,1)=-(a_lat/(2.0d0*pi))*(energy_full_3D(i1+1,i2,i3,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_x
            grad_vel_v_full(i1,i2,i3,j_spin,j,2)=-(a_lat/(2.0d0*pi))*(energy_full_3D(i1,i2+1,i3,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_y
            grad_vel_v_full(i1,i2,i3,j_spin,j,3)=-(a_lat/(2.0d0*pi))*(energy_full_3D(i1,i2,i3+1,j_spin,j)-energy_full_3D(i1,i2,i3,j_spin,j))/delta_k_z

            write(1,*) i1,i2,i3, grad_vel_v_full(i1,i2,i3,j_spin,j,1)

                   end do
                end do
             end do
          end do
      end do
      close(1)

      open(1,File='data_BTO/weight.dat')
      do k=1,n_k
         read(1,*) weight(k)
      end do
      close(1)

      open(1,File='data_BTO/kpoints.dat')
      do k=1,n_k
         read(1,*) k_m(1,k), k_m(2,k), k_m(3,k)
      end do
      close(1)



     do j=n_c,n_c_max
         do j_spin=1,n_spin
            do i=1,n_k

            do i1=-n_k_1_full_max,n_k_1_full_max-1
               do i2=-n_k_2_full_max,n_k_2_full_max-1
                  do i3=-n_k_3_full_max,n_k_3_full_max-1

                  k_x_full=i1*delta_k_x
                  k_y_full=i2*delta_k_y
                  k_z_full=i3*delta_k_z


                  delta_k_map=dsqrt((k_x_full-k_m(1,i))**2.0d0+(k_y_full-k_m(2,i))**2.0d0+(k_z_full-k_m(3,i))**2.0d0)
!                   write(*,*) delta_k_map
                  if(delta_k_map.lt.0.01d0) then

                  grad_vel_c(i,j_spin,j,1)=grad_vel_c_full(i1,i2,i3,j_spin,j,1)
                  grad_vel_c(i,j_spin,j,2)=grad_vel_c_full(i1,i2,i3,j_spin,j,2)
                  grad_vel_c(i,j_spin,j,3)=grad_vel_c_full(i1,i2,i3,j_spin,j,3)

                  end if

                     end do
                  end do
               end do
              write(*,*) grad_vel_c(i,j_spin,j,1)


            end do
         end do
      end do


     do j=n_v_min,n_v
         do j_spin=1,n_spin
         do i=1,n_k


            do i1=-n_k_1_full_max,n_k_1_full_max-1
               do i2=-n_k_2_full_max,n_k_2_full_max-1
                  do i3=-n_k_3_full_max,n_k_3_full_max-1

                  k_x_full=i1*delta_k_x
                  k_y_full=i2*delta_k_y
                  k_z_full=i3*delta_k_z

                  delta_k_map=dsqrt((k_x_full-k_m(1,i))**2.0d0+(k_y_full-k_m(2,i))**2.0d0+(k_z_full-k_m(3,i))**2.0d0)
!                    write(*,*) delta_k_map

                  if(delta_k_map.lt.0.01d0) then

                  grad_vel_v(i,j_spin,j,1)=grad_vel_v_full(i1,i2,i3,j_spin,j,1)
                  grad_vel_v(i,j_spin,j,2)=grad_vel_v_full(i1,i2,i3,j_spin,j,2)
                  grad_vel_v(i,j_spin,j,3)=grad_vel_v_full(i1,i2,i3,j_spin,j,3)

                  end if

                     end do
                  end do
               end do

            write(*,*) grad_vel_v(i,j_spin,j,1)

            end do
         end do
      end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do k=1,n_k-1
           delta_kx(k)=k_m(1,k+1)-k_m(1,k)
           delta_ky(k)=k_m(2,k+1)-k_m(2,k)
           delta_kz(k)=k_m(3,k+1)-k_m(3,k)
           delta_k(k)=dsqrt(delta_kx(k)**2.0d0+delta_ky(k)**2.0d0+delta_kz(k)**2.0d0)
           alpha_x(k)=delta_kx(k)/delta_k(k)
           alpha_y(k)=delta_ky(k)/delta_k(k)
           alpha_z(k)=delta_kz(k)/delta_k(k)
         end do


           delta_kx(n_k)=0.0d0
           delta_ky(n_k)=0.0d0
           delta_kz(n_k)=0.0d0
           delta_k(n_k)=1000.0d0
           alpha_x(n_k)=0.0d0
           alpha_z(n_k)=0.0d0


!CCCCCCCCCCCCC read fxc CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      call execute_command_line('mkdir rho')
      call execute_command_line('mkdir J')


      open(1,File='J/E')
      do i=-n_time_min,n_time_max
         time=delta_t*i

!           Efield(i)=Efield_max*dexp(-(time**2.0d0)/(tau_pulse**2.0d0))*dcos(omega0*time)


!             Efield(i)=Efield_max*(1.0d0/(dcosh(time/tau_pulse)**2.0d0))*dcos(omega0*time)

             Efield(i)=0.0d0
             if((time-tau_pulse_c)*(time-tau_pulse_c).le.tau_pulse*tau_pulse)                                                   &
             Efield(i)=Efield_max*(dcos(0.5d0*pi*(time-tau_pulse_c)/tau_pulse)**2.0d0)*dcos(omega0*(time-tau_pulse_c))

!              Efield(i)=Efield_max*dexp(-((time-tau_pulse_c)**2.0d0)/(tau_pulse**2.0d0))*dcos(omega0*(time-tau_pulse_c))

         write(1,*) time, Efield(i)
      end do
      close(1)


     open(1,File='J/A')
     do i=-n_time_min,n_time_max
         time=delta_t*i
        Afield(i)=0.0d0
         do j=-n_time_min,i
           Afield(i)=Afield(i)-Efield(j)*delta_t
         end do
         write(1,*) time, Afield(i)
      end do
      close(1)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       open(1,File='data/fxc_omega_U4_oneloop')
       do k=1,2*n_w_max+1

          read(1,*) omega, x,y
          fxc_0(k-n_w_max-1)=dcmplx(x,y)
       end do
       close(1)


       open(1,File='fxc_time')
       do i=0,2*n_time_max
          time=delta_t*i
          fxc_0_time(i)=0.0d0
          do k=-n_w_max/10,n_w_max/10-1
          omega=delta_omega*k
           fxc_0_time(i)=fxc_0_time(i)+cdexp(-ii*omega*time)*fxc_0(k)*delta_omega


          end do
          fxc_0_time(i)=0.01d0*(1.0d0/(2.0d0*pi))*fxc_0_time(i)
          write(1,*) 0.658d0*time, dreal(fxc_0_time(i)), dimag(fxc_0_time(i))
       end do
       close(1)
!         fxc_0_time(0)=fxc_0_time(0)-1.0d0*(U_eff/delta_t)
!          fxc_0_time(0)=dreal(fxc_0(0))*(1.0d0-0.025*ii)-U_eff
        fxc_0_time(0)=dreal(fxc_0(0))-U_eff
       write(*,*) fxc_0_time(0)
!         stop
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




     open(1,File='data_BTO/energy_all')
     do i=1,n_k
        do j=1,n_bands
           do j_spin=1,n_spin
              read(1,*) energy(i,j_spin,j)
           end do
        end do
     end do
     close(1)


      do j=n_v_min,n_v
         do j_spin=1,n_spin
         do i=1,n_k
            energy_v(i,j_spin,j)=energy(i,j_spin,j)
         end do
         end do
      end do


      do j=n_c,n_c_max
         do j_spin=1,n_spin
         do i=1,n_k
            energy_c(i,j_spin,j)=energy(i,j_spin,j)
         end do
         end do
      end do



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!       do j=n_c,n_c_max
!           do j_spin=1,n_spin
!           do i=1,n_k-1
!              grad_vel_c(i,j_spin,j,1)=alpha_grad*((energy_c(i+1,j_spin,j)-energy_c(i,j_spin,j))/delta_k(i))*alpha_x(i)
!              grad_vel_c(i,j_spin,j,2)=alpha_grad*((energy_c(i+1,j_spin,j)-energy_c(i,j_spin,j))/delta_k(i))*alpha_y(i)
!              grad_vel_c(i,j_spin,j,3)=alpha_grad*((energy_c(i+1,j_spin,j)-energy_c(i,j_spin,j))/delta_k(i))*alpha_z(i)
!           end do
!           end do
!        end do

!       do j=n_c,n_c_max
!           do j_spin=1,n_spin
!              do i1=1,3
!              grad_vel_c(n_k,j_spin,j,i1)=0.0d0
!              end do
!           end do
!        end do




!       do j=n_v_min,n_v
!           do j_spin=1,n_spin
!           do i=1,n_k-1
!              grad_vel_v(i,j_spin,j,1)=-alpha_grad*((energy_v(i+1,j_spin,j)-energy_v(i,j_spin,j))/delta_k(i))*alpha_x(i)
!              grad_vel_v(i,j_spin,j,2)=-alpha_grad*((energy_v(i+1,j_spin,j)-energy_v(i,j_spin,j))/delta_k(i))*alpha_y(i)
!              grad_vel_v(i,j_spin,j,3)=-alpha_grad*((energy_v(i+1,j_spin,j)-energy_v(i,j_spin,j))/delta_k(i))*alpha_z(i)
!           end do
!           end do
!        end do

!       do j=n_v_min,n_v
!           do j_spin=1,n_spin
!              do i1=1,3
!              grad_vel_v(n_k,j_spin,j,i1)=0.0d0
!              end do
!           end do
!        end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!      open(1,File='./G/G_kq.dat')
!
!      do i=1,n_k
!         do m=1,n_c_orb
!            do m1=1,n_spin
!
!               do j=1,n_k
!                  do l=1,n_v_orb
!                     do l1=1,n_spin
!
!               read(1,*) G(m,l,m1,l1,i,j)
!               G(m,l,m1,l1,i,j)=700.0d0*G(m,l,m1,l1,i,j)
!     multiply by volume factor in G (delta function))
!                    G(m,l,m1,l1,i,j)=13820.974d0*G(m,l,m1,l1,i,j)
!                     end do
!                  end do
!               end do
!            end do
!         end do
!      end do
!      close(1)


    do i=1,n_k
       do j=1,n_k

        do m1=n_v_min,n_c_max
           do m_spin_1=1,n_spin

              do m2=n_v_min,n_c_max
                 do m_spin_2=1,n_spin

                    do n1=n_v_min,n_c_max
                       do n_spin_1=1,n_spin

                          do n2=n_v_min,n_c_max
                             do n_spin_2=1,n_spin



!                G(m1,m2,m_spin_1,m_spin_2,n1,n2,n_spin_1,n_spin_2,i,j)=1000000.0d0*(1.0d0+abs(i-j)/n_k)
               G(m1,m2,m_spin_1,m_spin_2,n1,n2,n_spin_1,n_spin_2,i,j)= (-1)*num3 ! -140*0.086


                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do






       do i=1,n_k
          do m=1,n_bands
             do m1=1,n_spin

                do l=1,n_bands
                   do l1=1,n_spin

    dipole_0(m,l,m1,l1,i,1)=0.0d0
    dipole_0(m,l,m1,l1,i,2)=0.0d0
    dipole_0(m,l,m1,l1,i,3)=0.0d0

                   end do
                end do

          end do
       end do
    end do




    open(1,File='data_BTO/dipole_BTO.dat')

       do i=1,n_k
          do m=n_v,n_v
             do m1=1,n_spin

                do l=n_c,n_c
                   do l1=1,n_spin

    read(1,*) dipole_0(m,l,m1,l1,i,1), dipole_0(m,l,m1,l1,i,2), dipole_0(m,l,m1,l1,i,3)
!    read(1,*)
!    read(1,*) 


    dipole_0(m,l,m1,l1,i,1)=0.529d0*ii*dipole_0(m,l,m1,l1,i,1)
    dipole_0(m,l,m1,l1,i,2)=0.529d0*ii*dipole_0(m,l,m1,l1,i,2)
    dipole_0(m,l,m1,l1,i,3)=0.529d0*ii*dipole_0(m,l,m1,l1,i,3)


!         dipole_0(m,l,m1,l1,i,1)=dsqrt( dreal(dipole_0(m,l,m1,l1,i,1))**2.0d0              &
!                                      +dimag(dipole_0(m,l,m1,l1,i,1))**2.0d0)
!         dipole_0(m,l,m1,l1,i,2)=dsqrt( dreal(dipole_0(m,l,m1,l1,i,2))**2.0d0              &
!                                       +dimag(dipole_0(m,l,m1,l1,i,2))**2.0d0)
!         dipole_0(m,l,m1,l1,i,3)=dsqrt( dreal(dipole_0(m,l,m1,l1,i,3))**2.0d0              &
!                                       +dimag(dipole_0(m,l,m1,l1,i,3))**2.0d0)


!        dipole_0(m,l,m1,l1,i,1)=ii*1.0d0
!        dipole_0(m,l,m1,l1,i,2)=ii*1.0d0
!        dipole_0(m,l,m1,l1,i,3)=ii*1.0d0



                   end do
                end do

          end do
       end do
    end do
    close(1)


       do i=1,n_k
             do m1=1,n_spin
                   do l1=1,n_spin



         dipole_0(n_c,n_v,m1,l1,i,1)=dreal(dipole_0(n_v,n_c,m1,l1,i,1))              &
                                      -ii*dimag(dipole_0(n_v,n_c,m1,l1,i,1))

        dipole_0(n_c,n_v,m1,l1,i,2)=dreal(dipole_0(n_v,n_c,m1,l1,i,2))              &
                                      -ii*dimag(dipole_0(n_v,n_c,m1,l1,i,2))

        dipole_0(n_c,n_v,m1,l1,i,3)=dreal(dipole_0(n_v,n_c,m1,l1,i,3))              &
                                      -ii*dimag(dipole_0(n_v,n_c,m1,l1,i,3))

                   end do

          end do
    end do



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      open(12,File='J/I_omega')
      open(14,File='J/HH_11')
      open(201,File='J/j_x')
      open(202,File='J/j_y')
      open(203,File='J/j_z')
      open(71,File='J/density_v')
      open(72,File='J/density_c')
      open(1411,File='J/j_v_x')
      open(1412,File='J/j_v_y')
      open(1413,File='J/j_v_z')
      open(1421,File='J/j_c_x')
      open(1422,File='J/j_c_y')
      open(1423,File='J/j_c_z')
      open(1431,File='J/j_vv_x')
      open(1432,File='J/j_vv_y')
      open(1433,File='J/j_vv_z')
      open(1441,File='J/j_vc_x')
      open(1442,File='J/j_vc_y')
      open(1443,File='J/j_vc_z')
      open(1451,File='J/j_cv_x')
      open(1452,File='J/j_cv_y')
      open(1453,File='J/j_cv_z')
      open(1461,File='J/j_cc_x')
      open(1462,File='J/j_cc_y')
      open(1463,File='J/j_cc_z')
      do iiii=0, n_phi
!       do iiii=n_phi,n_phi
!
      phi=(pi/2.0d0/n_phi)*iiii



      do i=-n_time_min,n_time_max
          E_vec(1,i)=Efield(i)*dsin(phi)
          E_vec(2,i)=Efield(i)*0.0d0
          E_vec(3,i)=Efield(i)*dcos(phi)
      end do



      do i=-n_time_min,n_time_max
          A_vec(1,i)=Afield(i)*dsin(phi)
          A_vec(2,i)=Afield(i)*0.0d0
          A_vec(3,i)=Afield(i)*dcos(phi)
      end do


       do i=1,n_k
          do m=1,n_bands
             do m1=1,n_spin

                do l=1,n_bands
                   do l1=1,n_spin


    dipole(m,l,m1,l1,i,1)=dipole_0(m,l,m1,l1,i,1)*dsin(phi)
    dipole(m,l,m1,l1,i,2)=dipole_0(m,l,m1,l1,i,2)*0.0d0
    dipole(m,l,m1,l1,i,3)=dipole_0(m,l,m1,l1,i,3)*dcos(phi)


                   end do
                end do

          end do
       end do
    end do


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=1,n_k
          do m=n_v_min,n_v
             do m1=1,n_spin

                do l=n_c,n_c_max
                   do l1=1,n_spin

    dipole_vc_total(m,l,m1,l1,i)= dipole(m,l,m1,l1,i,1)                               &
                                 +dipole(m,l,m1,l1,i,2)                               &
                                 +dipole(m,l,m1,l1,i,3)


                   end do
                end do

          end do
       end do
    end do




       do i=1,n_k
          do m=n_c,n_c_max
             do m1=1,n_spin

                do l=n_v_min,n_v
                   do l1=1,n_spin

    dipole_cv_total(m,l,m1,l1,i)= dipole(m,l,m1,l1,i,1)                                &
                                 +dipole(m,l,m1,l1,i,2)                                &
                                 +dipole(m,l,m1,l1,i,3)


                   end do
                end do

          end do
       end do
    end do


       do i=1,n_k
          do m=n_v_min,n_v
             do m1=1,n_spin

                do l=n_v_min,n_v
                   do l1=1,n_spin

    dipole_vv_total(m,l,m1,l1,i)= dipole(m,l,m1,l1,i,1)                                &
                                 +dipole(m,l,m1,l1,i,2)                                &
                                 +dipole(m,l,m1,l1,i,3)


                   end do
                end do

          end do
       end do
    end do


       do i=1,n_k
          do m=n_c,n_c_max
             do m1=1,n_spin

                do l=n_c,n_c_max
                   do l1=1,n_spin

    dipole_cc_total(m,l,m1,l1,i)= dipole(m,l,m1,l1,i,1)                               &
                                 +dipole(m,l,m1,l1,i,2)                               &
                                 +dipole(m,l,m1,l1,i,3)

                   end do
                end do

          end do
       end do
    end do

!CCCCCCCCCCCCCCCCC define initial occupations CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


             do i=n_c,n_c_max
                do m=1,n_spin
                   do j=n_v_min,n_v
                      do l=1,n_spin
                         do k=1,n_k

             rho_cv_1_m1_1(i,j,m,l,k)=0.0d0
             rho_cv_m1_1(i,j,m,l,k)=0.0d0
             B_cv(i,j,m,l,k)=0.0d0

                         end do
                      end do
                   end do
                end do
             end do


             do i=n_v_min,n_v
                do m=1,n_spin
                   do j=n_c,n_c_max
                      do l=1,n_spin
                         do k=1,n_k

             rho_vc_1_m1_1(i,j,m,l,k)=0.0d0
             rho_vc_m1_1(i,j,m,l,k)=0.0d0
             B_vc(i,j,m,l,k)=0.0d0

                         end do
                      end do
                   end do
                end do
             end do




             do i=n_c,n_c_max
                do m=1,n_spin
                   do j=n_c,n_c_max
                      do l=1,n_spin
                         do k=1,n_k

             rho_cc_1_m1_1(i,j,m,l,k)=0.0d0
             rho_cc_m1_1(i,j,m,l,k)=0.0d0
             B_cc(i,j,m,l,k)=0.0d0

                         end do
                      end do
                   end do
                end do
             end do



             do i=n_v_min,n_v
                do m=1,n_spin
                   do j=n_v_min,n_v
                      do l=1,n_spin
                         do k=1,n_k

             rho_vv_1_m1_1(i,j,m,l,k)=0.0d0
             rho_vv_m1_1(i,j,m,l,k)=0.0d0
             B_vv(i,j,m,l,k)=0.0d0

                         end do
                      end do
                   end do
                end do
             end do


             do i=n_c,n_c_max
                do m=1,n_spin
                   do k=1,n_k

             rho_c_1_m1_1(i,m,k)=0.0d0
             rho_c_m1_1(i,m,k)=0.0d0

                   end do
                end do
             end do


             do i=n_v_min,n_v
                do m=1,n_spin
                   do k=1,n_k

                      rho_vv_1_m1_1(i,i,m,m,k)=1.0d0
                      rho_vv_m1_1(i,i,m,m,k)=1.0d0


                      rho_v_1_m1_1(i,m,k)=1.0d0
                      rho_v_m1_1(i,m,k)=1.0d0

                   end do
                end do
             end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    open(1,File='rho/rho_vc')
    open(2,File='rho/rho_cv')
    open(3,File='rho/rho_cc')
    open(4,File='rho/rho_vv')
    open(5,File='rho/rho_c')
    open(6,File='rho/rho_v')
    do m1=-n_time_min+1,n_time_max

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   rho_vc, rho_cv
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do j=n_c,n_c_max
                do j_spin=1,n_spin
                   do k=1,n_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        sum_1_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_1=sum_1_1+a1*ii*Afield(m1-1)*( dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,j,l_spin,j_spin,k)           &
                                                         -dipole_cc_total(l,j,l_spin,j_spin,k)*rho_vc_m1_1(i,l,i_spin,l_spin,k))                &
                                                    +( B_vc(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,j,l_spin,j_spin,k)                               &
                                                      -B_cc(l,j,l_spin,j_spin,k)*rho_vc_m1_1(i,l,i_spin,l_spin,k))

                           end do
                        end do


!                        write(*,*) sum_1_1


                        sum_2_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_1=sum_2_1+a1*ii*Afield(m1-1)*( dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,j,l_spin,j_spin,k)-          &
                                                         -dipole_vc_total(l,j,l_spin,j_spin,k)*rho_vv_m1_1(i,l,i_spin,l_spin,k))                &
                                                        +( B_vv(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,j,l_spin,j_spin,k)                           &
                                                          -B_vc(l,j,l_spin,j_spin,k)*rho_vv_m1_1(i,l,i_spin,l_spin,k))

                           end do
                        end do





            F_vc_1=(energy_v(k,i_spin,i)-energy_c(k,j_spin,j)+gap+Ueff-ii*(1.0d0/T_2))*rho_vc_m1_1(i,j,i_spin,j_spin,k)        &
                    +sum_1_1+sum_2_1                                                                                           &
                    +(rho_c_m1_1(j,j_spin,k)-rho_v_m1_1(i,i_spin,k))*ii*Afield(m1-1)*dipole_vc_total(i,j,i_spin,j_spin,k)      &
                    +(rho_c_m1_1(j,j_spin,k)-rho_v_m1_1(i,i_spin,k))*B_vc(i,j,i_spin,j_spin,k)                                 &
                    +(B_vv(i,i,i_spin,i_spin,k)-B_cc(j,j,j_spin,j_spin,k))*rho_vc_m1_1(i,j,i_spin,j_spin,k)


            rho_vc_1_m1(i,j,i_spin,j_spin,k)=rho_vc_m1_1(i,j,i_spin,j_spin,k)-ii*delta_t*F_vc_1
            rho_cv_1_m1(j,i,j_spin,i_spin,k)=dcmplx(dreal(rho_vc_1_m1(i,j,i_spin,j_spin,k)),-dimag(rho_vc_1_m1(i,j,i_spin,j_spin,k)))

!             write(*,*) rho_vc_1(i,j,i_spin,j_spin,k,m1), rho_cv_1(j,i,j_spin,i_spin,k,m1)

                    end do
                end do
             end do
          end do
       end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC rho_cc_1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=n_c,n_c_max
          do i_spin=1,n_spin
             do j=n_c,n_c_max
                do j_spin=1,n_spin
                   do k=1,n_k

                        sum_1_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_1=sum_1_1+a1*ii*Afield(m1-1)*dipole_cc_total(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,j,l_spin,j_spin,k) &
                                   +B_cc(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,j,l_spin,j_spin,k)




                           end do
                        end do



                        sum_2_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_2_1=sum_2_1-a1*ii*Afield(m1-1)*dipole_cc_total(l,j,l_spin,j_spin,k)*rho_cc_m1_1(i,l,i_spin,l_spin,k) &
                                   -B_cc(l,j,l_spin,j_spin,k)*rho_cc_m1_1(i,l,i_spin,l_spin,k)




                           end do
                        end do


                        sum_3_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_3_1=sum_3_1+a1*ii*Afield(m1-1)*( dipole_cv_total(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,j,l_spin,j_spin,k)-          &
                                                         -dipole_vc_total(l,j,l_spin,j_spin,k)*rho_cv_m1_1(i,l,i_spin,l_spin,k))                &
                                          +( B_cv(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,j,l_spin,j_spin,k)-                                        &
                                            -B_vc(l,j,l_spin,j_spin,k)*rho_cv_m1_1(i,l,i_spin,l_spin,k))

                           end do
                        end do


                        F_cc_1= (energy_c(k,i_spin,i)-energy_c(k,j_spin,j)-ii*(1.0d0/T_2))*rho_cc_m1_1(i,j,i_spin,j_spin,k)             &
                               +sum_1_1+sum_2_1+sum_3_1                                                                                 &
                               +(rho_c_m1_1(j,j_spin,k)-rho_c_m1_1(i,i_spin,k))*ii*Afield(m1-1)*dipole_cc_total(i,j,i_spin,j_spin,k)    &
                               +(rho_c_m1_1(j,j_spin,k)-rho_c_m1_1(i,i_spin,k))*B_cc(i,j,i_spin,j_spin,k)


                        rho_cc_1_m1(i,j,i_spin,j_spin,k)=rho_cc_m1_1(i,j,i_spin,j_spin,k)-ii*delta_t*F_cc_1

                        if(dreal(rho_cc_1_m1(i,j,i_spin,j_spin,k))**2.0d0+dimag(rho_cc_1_m1(i,j,i_spin,j_spin,k))**2.0d0.gt.1.0d0) &
                        rho_cc_1_m1(i,j,i_spin,j_spin,k)=1.0d0
!              write(*,*) i,j,i_spin,j_spin,k,m1,rho_cc(i,j,i_spin,j_spin,k,m1)

                   end do
                end do
             end do
          end do
       end do


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC rho_vv_1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do j=n_v_min,n_v
                do j_spin=1,n_spin
                   do k=1,n_k

                        sum_1_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_1=sum_1_1+a1*ii*Afield(m1-1)*dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vv_m1_1(l,j,l_spin,j_spin,k) &
                                   +B_vv(i,l,i_spin,l_spin,k)*rho_vv_m1_1(l,j,l_spin,j_spin,k)
                           end do
                        end do



                        sum_2_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                         a1=1.0d0

                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0

                           sum_2_1=sum_2_1-a1*ii*Afield(m1-1)*dipole_vv_total(l,j,l_spin,j_spin,k)*rho_vv_m1_1(i,l,i_spin,l_spin,k)   &
                                   -B_vv(l,j,l_spin,j_spin,k)*rho_vv_m1_1(i,l,i_spin,l_spin,k)



                           end do
                        end do


                        sum_3_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_3_1=sum_3_1+a1*ii*Afield(m1-1)*( dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cv_m1_1(l,j,l_spin,j_spin,k)   &
                                                            -dipole_cv_total(l,j,l_spin,j_spin,k)*rho_vc_m1_1(i,l,i_spin,l_spin,k))     &
                                          +( B_vc(i,l,i_spin,l_spin,k)*rho_cv_m1_1(l,j,l_spin,j_spin,k)                                 &
                                            -B_cv(l,j,l_spin,j_spin,k)*rho_vc_m1_1(i,l,i_spin,l_spin,k))

                           end do
                        end do


                       F_vv_1= (energy_v(k,i_spin,i)-energy_v(k,j_spin,j)-ii*(1.0d0/T_2))*rho_vv_m1_1(i,j,i_spin,j_spin,k)              &
                              +sum_1_1+sum_2_1+sum_3_1                                                                                  &
                              +(rho_v_m1_1(j,j_spin,k)-rho_v_m1_1(i,i_spin,k))*ii*Afield(m1-1)*dipole_vv_total(i,j,i_spin,j_spin,k)     &
                              +(rho_v_m1_1(j,j_spin,k)-rho_v_m1_1(i,i_spin,k))*B_vv(i,j,i_spin,j_spin,k)


            rho_vv_1_m1(i,j,i_spin,j_spin,k)= rho_vv_m1_1(i,j,i_spin,j_spin,k)-ii*delta_t*F_vv_1
            if(dreal(rho_vv_1_m1(i,j,i_spin,j_spin,k))**2.0d0+dimag(rho_vv_1_m1(i,j,i_spin,j_spin,k))**2.0d0.gt.1.0d0)    &
                        rho_vv_1_m1(i,j,i_spin,j_spin,k)=1.0d0

!              write(*,*) i,j,i_spin,j_spin,k,m1,rho_vv(i,j,i_spin,j_spin,k,m1)

                  end do
                end do
             end do
          end do
       end do
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCC rho_c_1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=n_c,n_c_max
          do i_spin=1,n_spin
             do k=1,n_k

                        sum_1_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_1=sum_1_1+a1*ii*Afield(m1-1)*dipole_cc_total(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,i,l_spin,i_spin,k) &
                                   +B_cc(i,l,i_spin,l_spin,k)*rho_cc_m1_1(l,i,l_spin,i_spin,k)

                           end do
                        end do

                        sum_2_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_1=sum_2_1+a1*ii*Afield(m1-1)*dipole_cv_total(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,i,l_spin,i_spin,k) &
                                   +B_cv(i,l,i_spin,l_spin,k)*rho_vc_m1_1(l,i,l_spin,i_spin,k)

                           end do
                        end do




                       F_c_1=2.0d0*dimag(sum_1_1)+2.0d0*dimag(sum_2_1)+2.0d0*dimag(B_cc(i,i,i_spin,i_spin,k)*rho_c_m1_1(i,i_spin,k))


                       rho_c_1_m1(i,i_spin,k)= rho_c_m1_1(i,i_spin,k)+delta_t*F_c_1
                       if(dreal(rho_c_1_m1(i,i_spin,k))**2.0d0+dimag(rho_c_1_m1(i,i_spin,k))**2.0d0.gt.1.0d0)            &
                       rho_c_1_m1(i,i_spin,k)=1.0d0


!              write(*,*) i,i_spin,k,m1,rho_c(i,i_spin,k,m1)
             end do
          end do
       end do


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC rho_v_1

       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do k=1,n_k

                        sum_1_1=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0

                           sum_1_1=sum_1_1+a1*ii*Afield(m1-1)*dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vv_m1_1(l,i,l_spin,i_spin,k)  &
                                 +B_vv(i,l,i_spin,l_spin,k)*rho_vv_m1_1(l,i,l_spin,i_spin,k)

                           end do
                        end do

                        sum_2_1=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_1=sum_2_1+a1*ii*Afield(m1-1)*dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cv_m1_1(l,i,l_spin,i_spin,k)    &
                                   +B_vc(i,l,i_spin,l_spin,k)*rho_cv_m1_1(l,i,l_spin,i_spin,k)


                           end do
                        end do


                       F_v_1=2.0d0*dimag(sum_1_1)+2.0d0*dimag(sum_2_1)+2.0d0*dimag(B_vv(i,i,i_spin,i_spin,k)*rho_v_m1_1(i,i_spin,k))


                       rho_v_1_m1(i,i_spin,k)= rho_v_m1_1(i,i_spin,k)+delta_t*F_v_1
                       if(dreal(rho_v_1_m1(i,i_spin,k))**2.0d0+dimag(rho_v_1_m1(i,i_spin,k))**2.0d0.gt.1.0d0)         &
                        rho_v_1_m1(i,i_spin,k)=1.0d0

!              write(*,*) i,i_spin,k,m1,rho_v(i,i_spin,k,m1)
             end do
          end do
       end do

       call sub_B(n_v_min,n_v,n_c,n_c_max,n_spin,n_k,fxc_0_time,G,weight,rho_vc_1_m1,rho_cv_1_m1,rho_cc_1_m1,rho_vv_1_m1,rho_c_1_m1,rho_v_1_m1,n_time_min,n_time_max,n_memory,m1,delta_t,delta_t_step,B_vc_1,B_cv_1,B_cc_1,B_vv_1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! rho_vc_2, rho_cv_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do j=n_c,n_c_max
                do j_spin=1,n_spin
                   do k=1,n_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        sum_1_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_2=sum_1_2+a1*ii*Afield(m1)*( dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,j,l_spin,j_spin,k)                          &
                                                             -dipole_cc_total(l,j,l_spin,j_spin,k)*rho_vc_1_m1(i,l,i_spin,l_spin,k))                         &
                                             +( B_vc_1(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,j,l_spin,j_spin,k)-                                                &
                                               -B_cc_1(l,j,l_spin,j_spin,k)*rho_vc_1_m1(i,l,i_spin,l_spin,k))




                           end do
                        end do

                        sum_2_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_2=sum_2_2+a1*ii*Afield(m1)*(dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,j,l_spin,j_spin,k)                             &
                                                         -dipole_vc_total(l,j,l_spin,j_spin,k)*rho_vv_1_m1(i,l,i_spin,l_spin,k))                               &
                                                +( B_vv_1(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,j,l_spin,j_spin,k)-                                               &
                                                  -B_vc_1(l,j,l_spin,j_spin,k)*rho_vv_1_m1(i,l,i_spin,l_spin,k))

                           end do
                        end do


                       F_vc_2=(energy_v(k,i_spin,i)-energy_c(k,j_spin,j)+gap+Ueff-ii*(1.0d0/T_2))*rho_vc_1_m1(i,j,i_spin,j_spin,k)                       &
                    +sum_1_2+sum_2_2                                                                                                                     &
                    +(rho_c_1_m1(j,j_spin,k)-rho_v_1_m1(i,i_spin,k))*ii*Afield(m1)*dipole_vc_total(i,j,i_spin,j_spin,k)                                  &
                    +(rho_c_1_m1(j,j_spin,k)-rho_v_1_m1(i,i_spin,k))*B_vc_1(i,j,i_spin,j_spin,k)                                                         &
                    +(B_vv_1(i,i,i_spin,i_spin,k)-B_cc_1(j,j,j_spin,j_spin,k))*rho_vc_1_m1(i,j,i_spin,j_spin,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                 rho_vc_m1(i,j,i_spin,j_spin,k)=rho_vc_m1_1(i,j,i_spin,j_spin,k)-0.5d0*ii*delta_t*(F_vc_1+F_vc_2)

            rho_cv_m1(j,i,j_spin,i_spin,k)=dcmplx(dreal(rho_vc_m1(i,j,i_spin,j_spin,k)),-dimag(rho_vc_m1(i,j,i_spin,j_spin,k)))

             write(1,*) rho_vc_m1(i,j,i_spin,j_spin,k)
             write(2,*) rho_cv_m1(j,i,j_spin,i_spin,k)

                   end do
                end do
             end do
          end do
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!! rho_cc_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=n_c,n_c_max
          do i_spin=1,n_spin
             do j=n_c,n_c_max
                do j_spin=1,n_spin
                   do k=1,n_k

                        sum_1_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_2=sum_1_2+a1*ii*Afield(m1)*dipole_cc_total(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,j,l_spin,j_spin,k)     &
                                   +B_cc_1(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,j,l_spin,j_spin,k)

                           end do
                        end do



                        sum_2_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_2_2=sum_2_2-a1*ii*Afield(m1)*dipole_cc_total(l,j,l_spin,j_spin,k)*rho_cc_1_m1(i,l,i_spin,l_spin,k)    &
                                   -B_cc_1(l,j,l_spin,j_spin,k)*rho_cc_1_m1(i,l,i_spin,l_spin,k)

                           end do
                        end do


                        sum_3_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_3_2=sum_3_2+a1*ii*Afield(m1)*(dipole_cv_total(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,j,l_spin,j_spin,k)   &
                                                         -dipole_vc_total(l,j,l_spin,j_spin,k)*rho_cv_1_m1(i,l,i_spin,l_spin,k))     &
                                   +  ( B_cv_1(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,j,l_spin,j_spin,k)                                 &
                                       -B_vc_1(l,j,l_spin,j_spin,k)*rho_cv_1_m1(i,l,i_spin,l_spin,k))
                           end do
                        end do


                        F_cc_2=(energy_c(k,i_spin,i)-energy_c(k,j_spin,j)-ii*(1.0d0/T_2))*rho_cc_1_m1(i,j,i_spin,j_spin,k)           &
                               +sum_1_2+sum_2_2+sum_3_2                                                                              &
                               +(rho_c_1_m1(j,j_spin,k)-rho_c_1_m1(i,i_spin,k))*ii*Afield(m1)*dipole_cc_total(i,j,i_spin,j_spin,k)   &
                               +(rho_c_1_m1(j,j_spin,k)-rho_c_1_m1(i,i_spin,k))*B_cc_1(i,j,i_spin,j_spin,k)

                        rho_cc_m1(i,j,i_spin,j_spin,k)=rho_cc_m1_1(i,j,i_spin,j_spin,k)-0.5d0*ii*delta_t*(F_cc_1+F_cc_2)
                        if(dreal(rho_cc_m1(i,j,i_spin,j_spin,k))**2.0d0+dimag(rho_cc_m1(i,j,i_spin,j_spin,k))**2.0d0.gt.1.0d0)      &
                        rho_cc_m1(i,j,i_spin,j_spin,k)=1.0d0


             write(3,*) rho_cc_m1(i,j,i_spin,j_spin,k)
!              write(*,*)
!              i,j,i_spin,j_spin,k,m1,rho_cc(i,j,i_spin,j_spin,k,m1)

                   end do
                end do
             end do
          end do
       end do


!!!!!!!!!!!!!!!!!!!! rho_vv_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do j=n_v_min,n_v
                do j_spin=1,n_spin
                   do k=1,n_k

                        sum_1_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.j) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_2=sum_1_2+a1*ii*Afield(m1)*dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vv_1_m1(l,j,l_spin,j_spin,k)   &
                                   +B_vv_1(i,l,i_spin,l_spin,k)*rho_vv_1_m1(l,j,l_spin,j_spin,k)

                           end do
                        end do



                        sum_2_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                         a1=1.0d0

                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0

                           sum_2_2=sum_2_2-a1*ii*Afield(m1)*dipole_vv_total(l,j,l_spin,j_spin,k)*rho_vv_1_m1(i,l,i_spin,l_spin,k)   &
                                   -B_vv_1(l,j,l_spin,j_spin,k)*rho_vv_1_m1(i,l,i_spin,l_spin,k)

                           end do
                        end do


                        sum_3_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_3_2=sum_3_2+a1*ii*Afield(m1)*( dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cv_1_m1(l,j,l_spin,j_spin,k) &
                                                          -dipole_cv_total(l,j,l_spin,j_spin,k)*rho_vc_1_m1(i,l,i_spin,l_spin,k))   &
                                   +( B_vc_1(i,l,i_spin,l_spin,k)*rho_cv_1_m1(l,j,l_spin,j_spin,k)                                  &
                                     -B_cv_1(l,j,l_spin,j_spin,k)*rho_vc_1_m1(i,l,i_spin,l_spin,k))
                           end do
                        end do


                       F_vv_2=(energy_v(k,i_spin,i)-energy_v(k,j_spin,j)-ii*(1.0d0/T_2))*rho_vv_1_m1(i,j,i_spin,j_spin,k)           &
                              +sum_1_2+sum_2_2+sum_3_2                                                                              &
                              +(rho_v_1_m1(j,j_spin,k)-rho_v_1_m1(i,i_spin,k))*ii*Afield(m1)*dipole_vv_total(i,j,i_spin,j_spin,k)   &
                              +(rho_v_1_m1(j,j_spin,k)-rho_v_1_m1(i,i_spin,k))*B_vv_1(i,j,i_spin,j_spin,k)

            rho_vv_m1(i,j,i_spin,j_spin,k)=rho_vv_m1_1(i,j,i_spin,j_spin,k)-0.5d0*ii*delta_t*(F_vv_1+F_vv_2)
            if(dreal(rho_vv_m1(i,j,i_spin,j_spin,k))**2.0d0+dimag(rho_vv_m1(i,j,i_spin,j_spin,k))**2.0d0.gt.1.0d0)                  &
                        rho_vv_m1(i,j,i_spin,j_spin,k)=1.0d0

             write(4,*) rho_vv_m1(i,j,i_spin,j_spin,k)
!              write(*,*)
!              i,j,i_spin,j_spin,k,m1,rho_vv(i,j,i_spin,j_spin,k,m1)

                  end do
                end do
             end do
          end do
       end do


!!!!!!!!!!!!!!! rho_c_2
!!!!!!!!!!!!!!!!!!!!!!!
       do i=n_c,n_c_max
          do i_spin=1,n_spin
             do k=1,n_k

                        sum_1_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.j.and.l_spin.eq.j_spin) a1=0.0d0

                           sum_1_2=sum_1_2+a1*ii*Afield(m1)*dipole_cc_total(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,i,l_spin,i_spin,k)  &
                                   +B_cc_1(i,l,i_spin,l_spin,k)*rho_cc_1_m1(l,i,l_spin,i_spin,k)

                           end do
                        end do

                        sum_2_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_2=sum_2_2+a1*ii*Afield(m1)*dipole_cv_total(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,i,l_spin,i_spin,k) &
                                   +B_cv_1(i,l,i_spin,l_spin,k)*rho_vc_1_m1(l,i,l_spin,i_spin,k)

                           end do
                        end do


                       F_c_2=2.0d0*dimag(sum_1_2)+2.0d0*dimag(sum_2_2)+2.0d0*dimag(B_cc_1(i,i,i_spin,i_spin,k)*rho_c_1_m1(i,i_spin,k))

                       rho_c_m1(i,i_spin,k)=rho_c_m1_1(i,i_spin,k)+0.5d0*delta_t*(F_c_1+F_c_2)
                       if(dreal(rho_c_m1(i,i_spin,k))**2.0d0+dimag(rho_c_m1(i,i_spin,k))**2.0d0.gt.1.0d0)            &
                       rho_c_m1(i,i_spin,k)=1.0d0

             write(5,*) rho_c_m1(i,i_spin,k)
!              write(*,*) i,i_spin,k,m1,rho_c(i,i_spin,k,m1)
             end do
          end do
       end do

!!!!!!!!!!!!!!!!!!! rho_v_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=n_v_min,n_v
          do i_spin=1,n_spin
             do k=1,n_k

                        sum_1_2=0.0d0
                        do l=n_v_min,n_v
                           do l_spin=1,n_spin

                           a1=1.0d0
                           if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0

                           sum_1_2=sum_1_2+a1*ii*Afield(m1)*dipole_vv_total(i,l,i_spin,l_spin,k)*rho_vv_1_m1(l,i,l_spin,i_spin,k)   &
                                   +a1*B_vv_1(i,l,i_spin,l_spin,k)*rho_vv_1_m1(l,i,l_spin,i_spin,k)

                           end do
                        end do

                        sum_2_2=0.0d0
                        do l=n_c,n_c_max
                           do l_spin=1,n_spin

                           a1=1.0d0
!                            if(l.eq.i) a1=0.0d0
!                            if(l.eq.i.and.l_spin.eq.i_spin) a1=0.0d0
                           sum_2_2=sum_2_2+a1*ii*Afield(m1)*dipole_vc_total(i,l,i_spin,l_spin,k)*rho_cv_1_m1(l,i,l_spin,i_spin,k)  &
                                   +a1*B_vc_1(i,l,i_spin,l_spin,k)*rho_cv_1_m1(l,i,l_spin,i_spin,k)
                           end do
                        end do


                       F_v_2=2.0d0*dimag(sum_1_2)+2.0d0*dimag(sum_2_2)


                       rho_v_m1(i,i_spin,k)= rho_v_m1_1(i,i_spin,k)+0.5d0*delta_t*(F_v_1+F_v_2)+2.0d0*dimag(B_vv_1(i,i,i_spin,i_spin,k)*rho_v_1_m1(i,i_spin,k))
                       if(dreal(rho_v_m1(i,i_spin,k))**2.0d0+dimag(rho_v_m1(i,i_spin,k))**2.0d0.gt.1.0d0)                                      &
                       rho_v_m1(i,i_spin,k)=1.0d0


             write(6,*) rho_v_m1(i,i_spin,k)
!              write(*,*) i,i_spin,k,m1,rho_v(i,i_spin,k,m1)
             end do
          end do
       end do
       call sub_B(n_v_min,n_v,n_c,n_c_max,n_spin,n_k,fxc_0_time,G,weight,rho_vc_m1,rho_cv_m1,rho_cc_m1,rho_vv_m1,rho_c_m1,rho_v_m1,n_time_min,n_time_max,n_memory,m1,delta_t,delta_t_step,B_vc,B_cv,B_cc,B_vv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            calculate densities, intra currents and polarizabilities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           density_v(m1)=0.0d0
           do i=n_v_min,n_v
              do i_spin=1,n_spin
                   do k=1,n_k
                    density_v(m1)=density_v(m1)+rho_v_m1(i,i_spin,k)*weight(k)
                   end do
              end do
           end do


           density_c(m1)=0.0d0
           do i=n_c,n_c_max
              do i_spin=1,n_spin
                   do k=1,n_k
                    density_c(m1)=density_c(m1)+rho_c_m1(i,i_spin,k)*weight(k)
                   end do
              end do
           end do





           do i1=1,3
           do j=n_v_min,n_v
              do j_spin=1,n_spin
                    j_intra_v(j,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     j_intra_v(j,j_spin,i1,m1)=j_intra_v(j,j_spin,i1,m1)+grad_vel_v(k,j_spin,j,i1)*(1.0d0-rho_v_m1(j,j_spin,k))*weight(k)

                   end do
              end do
           end do
           end do




           do i1=1,3
           do j=n_c,n_c_max
              do j_spin=1,n_spin
                    j_intra_c(j,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     j_intra_c(j,j_spin,i1,m1)=j_intra_c(j,j_spin,i1,m1)+grad_vel_c(k,j_spin,j,i1)*rho_c_m1(j,j_spin,k)*weight(k)

                   end do
              end do
           end do
           end do




           do i=n_v_min,n_v
              do i_spin=1,n_spin

              do j=n_v_min,n_v
                 do j_spin=1,n_spin
                    do i1=1,3
                    P_inter_vv(i,j,i_spin,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     P_inter_vv(i,j,i_spin,j_spin,i1,m1)=P_inter_vv(i,j,i_spin,j_spin,i1,m1)                               &
                                               +rho_vv_m1(i,j,i_spin,j_spin,k)*dipole(i,j,i_spin,j_spin,k,i1)*weight(k)

                   end do

                     P_inter_vv(i,j,i_spin,j_spin,i1,m1)=2.0d0*P_inter_vv(i,j,i_spin,j_spin,i1,m1)


                       end do
                    end do
                 end do
              end do
           end do



        do i=n_v_min,n_v
           do i_spin=1,n_spin

              do j=n_c,n_c_max
                 do j_spin=1,n_spin
                    do i1=1,3
                    P_inter_vc(i,j,i_spin,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     P_inter_vc(i,j,i_spin,j_spin,i1,m1)=P_inter_vc(i,j,i_spin,j_spin,i1,m1)                             &
                                               +rho_vc_m1(i,j,i_spin,j_spin,k)*dipole(i,j,i_spin,j_spin,k,i1)*weight(k)
                   end do

                   P_inter_vc(i,j,i_spin,j_spin,i1,m1)=2.0d0*P_inter_vc(i,j,i_spin,j_spin,i1,m1)


                    end do
                 end do
              end do
!                write(*,*) 'j',i1,m1,j_curr(i1,m1)
              end do
           end do



        do i=n_c,n_c_max
           do i_spin=1,n_spin

              do j=n_v_min,n_v
                 do j_spin=1,n_spin
                    do i1=1,3
                    P_inter_cv(i,j,i_spin,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     P_inter_cv(i,j,i_spin,j_spin,i1,m1)=P_inter_cv(i,j,i_spin,j_spin,i1,m1)                                   &
                                               +rho_cv_m1(i,j,i_spin,j_spin,k)*dipole(i,j,i_spin,j_spin,k,i1)*weight(k)

                   end do

                   P_inter_cv(i,j,i_spin,j_spin,i1,m1)=2.0d0*P_inter_cv(i,j,i_spin,j_spin,i1,m1)

                   end do
              end do
              end do
!                write(*,*) 'j',i1,m1,j_curr(i1,m1)
              end do
           end do



        do i=n_c,n_c_max
           do i_spin=1,n_spin

              do j=n_c,n_c_max
                 do j_spin=1,n_spin
                    do i1=1,3
                    P_inter_cc(i,j,i_spin,j_spin,i1,m1)=0.0d0

                   do k=1,n_k

                     P_inter_cc(i,j,i_spin,j_spin,i1,m1)=P_inter_cc(i,j,i_spin,j_spin,i1,m1)                                &
                                               +rho_cc_m1(i,j,i_spin,j_spin,k)*dipole(i,j,i_spin,j_spin,k,i1)*weight(k)

                   end do

                   P_inter_cc(i,j,i_spin,j_spin,i1,m1)=2.0d0*P_inter_cc(i,j,i_spin,j_spin,i1,m1)

                   end do
              end do
              end do
!                write(*,*) 'j',i1,m1,j_curr(i1,m1)
              end do
           end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                             redefine rho before moving back to the next time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             do i=n_c,n_c_max
                do m=1,n_spin
                   do j=n_v_min,n_v
                      do l=1,n_spin
                         do k=1,n_k

             rho_cv_1_m1_1(i,j,m,l,k)=rho_cv_1_m1(i,j,m,l,k)
             rho_cv_m1_1(i,j,m,l,k)=rho_cv_m1(i,j,m,l,k)


                         end do
                      end do
                   end do
                end do
             end do


             do i=n_v_min,n_v
                do m=1,n_spin
                   do j=n_c,n_c_max
                      do l=1,n_spin
                         do k=1,n_k

             rho_vc_1_m1_1(i,j,m,l,k)=rho_vc_1_m1(i,j,m,l,k)
             rho_vc_m1_1(i,j,m,l,k)=rho_vc_m1(i,j,m,l,k)

                         end do
                      end do
                   end do
                end do
             end do



             do i=n_c,n_c_max
                do m=1,n_spin
                   do j=n_c,n_c_max
                      do l=1,n_spin
                         do k=1,n_k

             rho_cc_1_m1_1(i,j,m,l,k)=rho_cc_1_m1(i,j,m,l,k)
             rho_cc_m1_1(i,j,m,l,k)=rho_cc_m1(i,j,m,l,k)

                         end do
                      end do
                   end do
                end do
             end do



             do i=n_v_min,n_v
                do m=1,n_spin
                   do j=n_v_min,n_v
                      do l=1,n_spin
                         do k=1,n_k

             rho_vv_1_m1_1(i,j,m,l,k)= rho_vv_1_m1_1(i,j,m,l,k)
             rho_vv_m1_1(i,j,m,l,k)=rho_vv_m1(i,j,m,l,k)

                         end do
                      end do
                   end do
                end do
             end do



             do i=n_c,n_c_max
                do m=1,n_spin
                   do k=1,n_k

             rho_c_1_m1_1(i,m,k)=rho_c_1_m1(i,m,k)
             rho_c_m1_1(i,m,k)=rho_c_m1(i,m,k)

                   end do
                end do
             end do


             do i=n_v_min,n_v
                do m=1,n_spin
                   do k=1,n_k

                      rho_v_1_m1_1(i,m,k)=rho_v_1_m1(i,m,k)
                      rho_v_m1_1(i,m,k)=rho_v_m1(i,m,k)

                   end do
                end do
             end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     write(*,*) m1
    end do


    close(6)
    close(5)
    close(4)
    close(3)
    close(2)
    close(1)


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


        do m1=-n_time_min+1,n_time_max
           write(71,*) m1*delta_t,density_v(m1)
         end do


        do m1=-n_time_min+1,n_time_max
           write(72,*) m1*delta_t,density_c(m1)
         end do



        do m1=-n_time_min+1,n_time_max
           do i1=1,3
              j_tot(i1,m1)=0.0d0
           end do
        end do



        do m1=-n_time_min+1,n_time_max
           do i1=1,3
              Ax=0.0d0
           do j=n_v_min,n_v
              do j_spin=1,n_spin
                   Ax=Ax+j_intra_v(j,j_spin,i1,m1)

              end do
           end do
!                write(*,*) 'j',i1,m1,j_curr(i1,m1)
               j_tot(i1,m1)=j_tot(i1,m1)+0.0d0*Ax
               j_tot_intra_v(i1,m1)=Ax
!                write(*,*) m1, 'j_intra_v'
           end do
           write(1411,*) m1*delta_t, j_tot_intra_v(1,m1)
           write(1412,*) m1*delta_t, j_tot_intra_v(2,m1)
           write(1413,*) m1*delta_t, j_tot_intra_v(3,m1)
         end do



        do m1=-n_time_min+1,n_time_max
           do i1=1,3
              Ax=0.0d0
           do j=n_c,n_c_max
              do j_spin=1,n_spin
                   Ax=Ax+j_intra_c(j,j_spin,i1,m1)

              end do
           end do
              j_tot(i1,m1)=j_tot(i1,m1)+1.0d0*Ax


               j_tot_intra_c(i1,m1)=Ax
!                write(*,*) m1, 'j_intra_v'
           end do
           write(1421,*) m1*delta_t, j_tot_intra_c(1,m1)
           write(1422,*) m1*delta_t, j_tot_intra_c(2,m1)
           write(1423,*) m1*delta_t, j_tot_intra_c(3,m1)

        end do


!         write(*,*) 'j_intra ok'



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do m1=-n_time_min+1,n_time_max-1
           do i1=1,3
              Ax=0.0d0
              do i=n_v_min,n_v
                 do i_spin=1,n_spin

                    do j=n_v_min,n_v
                       do j_spin=1,n_spin

                     j_inter_vv(i,j,i_spin,j_spin,i1,m1)=(P_inter_vv(i,j,i_spin,j_spin,i1,m1+1)-P_inter_vv(i,j,i_spin,j_spin,i1,m1))/delta_t

                     Ax=Ax+j_inter_vv(i,j,i_spin,j_spin,i1,m1)

                       end do
                    end do

                 end do
              end do
              j_tot(i1,m1)=j_tot(i1,m1)+0.0d0*Ax
               j_tot_inter_vv(i1,m1)=Ax
!                write(*,*) m1, 'j_intra_v'
           end do
           write(1431,*) m1*delta_t, j_tot_inter_vv(1,m1)
           write(1432,*) m1*delta_t, j_tot_inter_vv(2,m1)
           write(1433,*) m1*delta_t, j_tot_inter_vv(3,m1)


!            write(*,*) m1, 'j_inter_vv'
        end do


        do m1=-n_time_min+1,n_time_max-1
           do i1=1,3
           Ax=0.0d0


        do i=n_v_min,n_v
           do i_spin=1,n_spin

              do j=n_c,n_c_max
                 do j_spin=1,n_spin


                     j_inter_vc(i,j,i_spin,j_spin,i1,m1)=(P_inter_vc(i,j,i_spin,j_spin,i1,m1+1)-P_inter_vc(i,j,i_spin,j_spin,i1,m1))/delta_t

                     Ax=Ax+j_inter_vc(i,j,i_spin,j_spin,i1,m1)
                    end do
                 end do
              end do
              end do
              j_tot(i1,m1)=j_tot(i1,m1)+0.0d0*Ax
!                write(*,*) 'j',i1,m1,j_curr(i1,m1)

!                write(*,*) m1, 'j_intra_v'
               j_tot_inter_vc(i1,m1)=Ax

           end do
           if(iiii.eq.n_phi) then
           write(1441,*) m1*delta_t, j_tot_inter_vc(1,m1)
           end if
           write(1442,*) m1*delta_t, j_tot_inter_vc(2,m1)
           write(1443,*) m1*delta_t, j_tot_inter_vc(3,m1)


!            write(*,*) m1, 'j_inter_vc'
        end do




        do m1=-n_time_min+1,n_time_max-1
           do i1=1,3
           Ax=0.0d0


        do i=n_c,n_c_max
           do i_spin=1,n_spin

              do j=n_v_min,n_v
                 do j_spin=1,n_spin


                     j_inter_cv(i,j,i_spin,j_spin,i1,m1)=(P_inter_cv(i,j,i_spin,j_spin,i1,m1+1)-P_inter_cv(i,j,i_spin,j_spin,i1,m1))/delta_t
                     Ax=Ax+j_inter_cv(i,j,i_spin,j_spin,i1,m1)

                    end do
                 end do
              end do
              end do
              j_tot(i1,m1)=j_tot(i1,m1)+0.0d0*Ax

               j_tot_inter_cv(i1,m1)=Ax

           end do
           write(1451,*) m1*delta_t, j_tot_inter_cv(1,m1)
           write(1452,*) m1*delta_t, j_tot_inter_cv(2,m1)
           write(1453,*) m1*delta_t, j_tot_inter_cv(3,m1)



!            write(*,*) m1, 'j_inter_cv'
        end do


        do m1=-n_time_min+1,n_time_max-1
           do i1=1,3
           Ax=0.0d0

        do i=n_c,n_c_max
           do i_spin=1,n_spin

              do j=n_c,n_c_max
                 do j_spin=1,n_spin


                     j_inter_cc(i,j,i_spin,j_spin,i1,m1)=(P_inter_cc(i,j,i_spin,j_spin,i1,m1+1)-P_inter_cc(i,j,i_spin,j_spin,i1,m1))/delta_t
                     Ax=Ax+j_inter_cc(i,j,i_spin,j_spin,i1,m1)

                    end do
                 end do
              end do
              end do
              j_tot(i1,m1)=j_tot(i1,m1)+0.0d0*Ax
!               write(*,*) m1,i1,j_tot(i1,m1)

               j_tot_inter_cc(i1,m1)=Ax

           end do
           write(1461,*) m1*delta_t, j_tot_inter_cc(1,m1)
           write(1462,*) m1*delta_t, j_tot_inter_cc(2,m1)
           write(1463,*) m1*delta_t, j_tot_inter_cc(3,m1)

        end do







        write(20,*) 'phi', (phi/pi)*180.0d0
        do m1=-n_time_min,n_time_max-1

           write(201,*) m1*delta_t, j_tot(1,m1)
           write(202,*) m1*delta_t, j_tot(2,m1)
           write(203,*) m1*delta_t, j_tot(3,m1)

        end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          do i1=1,3
              do m1=-n_time_min,n_time_max
                 current(m1)=j_tot(i1,m1)
              end do
                 call fast(current,n_time_max,delta_t,j_omega)

              do i=-n_time_max/2,n_time_max/2
                 j_tot_omega(i1,i)=j_omega(i)
              end do
          end do

!go to 100
          open(123321,File='J/j_tot_omega')
          do i=-n_time_max/2,n_time_max/2
                 write(123321,*) 4.0d0*pi*delta_omega*i/omega0,dreal(j_tot_omega(1,i))**2.0d0+dimag(j_tot_omega(1,i))**2.0d0, dreal(j_tot_omega(3,i))**2.0d0+dimag(j_tot_omega(3,i))**2.0d0
          end do
          close(123321)

!100 continue

!           do i=1,n_omega_max
          do i=2,n_omega_max/100
             omega=delta_omega*i

!               do i1=1,3
!
!               j_tot_omega(i1,i)=0.0d0
!
!                do m1=-n_time_min+1,n_time_max-1
!                   time=delta_t*m1
!
!                j_tot_omega(i1,i)=j_tot_omega(i1,i)+cdexp(ii*omega*time)*j_tot(i1,m1)*delta_t
!
!                end do
!             end do

!             I_int(i)= dreal(j_tot_omega(1,i)*dsin(phi)+j_tot_omega(2,i)*0.0d0+j_tot_omega(3,i)*dcos(phi))**2.0d0             &
!                      +dimag(j_tot_omega(1,i)*dsin(phi)+j_tot_omega(2,i)*0.0d0+j_tot_omega(3,i)*dcos(phi))**2.0d0

            I_int(i)= dreal(j_tot_omega(1,i))**2.0d0+dimag(j_tot_omega(1,i))**2.0d0               &
                     +dreal(j_tot_omega(2,i))**2.0d0+dimag(j_tot_omega(2,i))**2.0d0               &
                     +dreal(j_tot_omega(3,i))**2.0d0+dimag(j_tot_omega(3,i))**2.0d0




         write(12,*) 4.0d0*pi*omega/omega0,(phi/pi)*180.0d0,I_int(i)

!             write(*,*) 'w',omega,'phi',phi,I_int(i)
        end do


        i=4*omega0/(4.0d0*pi*delta_omega)

        i4=5*omega0/(4.0d0*pi*delta_omega)

        write(14,*) (phi/pi)*180.0d0, I_int(i), I_int(i4), 6*omega0/(4.0d0*pi*delta_omega)
        write(*,*) (phi/pi)*180.0d0
     end do


     close(1463)
     close(1462)
     close(1461)


     close(1453)
     close(1452)
     close(1451)

     close(1443)
     close(1442)
     close(1441)

     close(1433)
     close(1432)
     close(1431)

     close(1423)
     close(1422)
     close(1421)

     close(1413)
     close(1412)
     close(1411)
     close(72)
     close(71)
     close(203)
     close(202)
     close(201)
     close(14)
     close(12)

     include 'saving.f90'

    print *, Efield_max, n_time_max, T_2, num3
    print *, char(7)  !will make a sound after the completion of running
    call cpu_time(stop_time)
    print *, "Original loop time:", stop_time - start_time, "seconds", (stop_time - start_time)/60, "minutes"
    stop
    end
