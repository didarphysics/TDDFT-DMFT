subroutine sub_B(n_v_min,n_v,n_c,n_c_max,n_spin,n_k,fxc_0_time,G,weight,rho_vc,rho_cv,rho_cc,rho_vv,rho_c,rho_v,n_time_min,n_time_max,n_memory,m1,delta_t,delta_t_step,B_vc,B_cv,B_cc,B_vv)

implicit real*8 (a-z)

integer*8 :: n_v_min,n_v,n_c,n_c_max,n_spin,n_k,n_time_min,n_time_max,n_memory,delta_t_step,m1

integer*8 :: i, i_1, i_spin, i_spin_1, j, j_1, j_spin, j_spin_1,k, p

real*8, dimension (n_k) :: weight

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_vc
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_cv
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_vv
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_cc
complex*16, dimension (n_v_min:n_v,n_spin,n_k) :: rho_v
complex*16, dimension (n_c:n_c_max,n_spin,n_k) :: rho_c

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: U_vc
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: U_cv
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: U_vv
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: U_cc

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: B_vc
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: B_cv
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: B_vv
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: B_cc

complex*16, dimension (0:2*n_time_max) :: fxc_0_time
complex*16, dimension(n_v_min:n_c_max,n_v_min:n_c_max,n_spin,n_spin,n_v_min:n_c_max,n_v_min:n_c_max,n_spin,n_spin,n_k,n_k) :: G


complex*16 :: ii, Ax, Bx


external sub_U
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

call sub_U(n_v_min,n_v,n_c,n_c_max,n_spin,n_k,fxc_0_time,rho_vc,rho_cv,rho_cc,rho_vv,rho_c,rho_v,n_time_min,n_time_max,n_memory,m1,delta_t,delta_t_step,U_vc,U_cv,U_cc,U_vv)



    do i=n_v_min,n_v
       do i_spin=1,n_spin
          do j=n_c,n_c_max
             do j_spin=1,n_spin
                do k=1,n_k

         B_vc(i,j,i_spin,j_spin,k)=0.0d0

         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vc(i,j,i_spin,j_spin,k)=B_vc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do

         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vc(i,j,i_spin,j_spin,k)=B_vc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do

         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vc(i,j,i_spin,j_spin,k)=B_vc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do

         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vc(i,j,i_spin,j_spin,k)=B_vc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    do i=n_c,n_c_max
       do i_spin=1,n_spin
          do j=n_v_min,n_v
             do j_spin=1,n_spin
                do k=1,n_k

         B_cv(i,j,i_spin,j_spin,k)=0.0d0

         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cv(i,j,i_spin,j_spin,k)=B_cv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do




         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cv(i,j,i_spin,j_spin,k)=B_cv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do



         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cv(i,j,i_spin,j_spin,k)=B_cv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do



         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cv(i,j,i_spin,j_spin,k)=B_cv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    do i=n_c,n_c_max
       do i_spin=1,n_spin
          do j=n_c,n_c_max
             do j_spin=1,n_spin
                do k=1,n_k

         B_cc(i,j,i_spin,j_spin,k)=0.0d0

         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cc(i,j,i_spin,j_spin,k)=B_cc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do




         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cc(i,j,i_spin,j_spin,k)=B_cc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do



         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cc(i,j,i_spin,j_spin,k)=B_cc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do



         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_cc(i,j,i_spin,j_spin,k)=B_cc(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    do i=n_v_min,n_v
       do i_spin=1,n_spin
          do j=n_v_min,n_v
             do j_spin=1,n_spin
                do k=1,n_k

         B_vv(i,j,i_spin,j_spin,k)=0.0d0

         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vv(i,j,i_spin,j_spin,k)=B_vv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do




         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vv(i,j,i_spin,j_spin,k)=B_vv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do


         do i_1=n_c,n_c_max
            do i_spin_1=1,n_spin
               do j_1=n_c,n_c_max
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vv(i,j,i_spin,j_spin,k)=B_vv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_cc(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

                      end do
                   end do
                end do
             end do
          end do



         do i_1=n_v_min,n_v
            do i_spin_1=1,n_spin
               do j_1=n_v_min,n_v
                   do j_spin_1=1,n_spin
                      do p=1,n_k

            B_vv(i,j,i_spin,j_spin,k)=B_vv(i,j,i_spin,j_spin,k)+G(i,j,i_spin,j_spin,i_1,j_1,i_spin_1,j_spin_1,k,p)*U_vv(i_1,j_1,i_spin_1,j_spin_1,p)*weight(p)

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


return
end


subroutine sub_U(n_v_min,n_v,n_c,n_c_max,n_spin,n_k,fxc_0_time,rho_vc,rho_cv,rho_cc,rho_vv,rho_c,rho_v,n_time_min,n_time_max,n_memory,m1,delta_t,delta_t_step,U_vc,U_cv,U_cc,U_vv)

implicit real*8 (a-z)

integer*8 :: n_v_min,n_v,n_c,n_c_max,n_spin,n_k,n_time_min,n_time_max,n_memory,delta_t_step,m1

integer*8 :: i,i1,j,l,m,k,i_time


complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_vc
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_cv
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_vv
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_cc
complex*16, dimension (n_v_min:n_v,n_spin,n_k) :: rho_v
complex*16, dimension (n_c:n_c_max,n_spin,n_k) :: rho_c

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) ::  U_vc
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) ::  U_cv
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) ::  U_vv
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) ::  U_cc

complex*16, dimension (0:2*n_time_max) :: fxc_0_time

complex*16 :: ii, Ax, Bx



        do i=n_c,n_c_max
           do m=1,n_spin
              do k=1,n_k

                    rho_cc(i,i,m,m,k)=rho_c(i,m,k)

              end do
           end do
        end do




        do i=n_v_min,n_v
           do m=1,n_spin
              do k=1,n_k

                    rho_vv(i,i,m,m,k)=rho_v(i,m,k)

              end do
           end do
        end do







         do i=n_v_min,n_v
            do m=1,n_spin
               do j=n_c,n_c_max
                  do l=1,n_spin
                     do k=1,n_k


               U_vc(i,j,m,l,k)=fxc_0_time(0)*rho_vc(i,j,m,l,k)


                     end do
                  end do
               end do
            end do
         end do



         do i=n_c,n_c_max
            do m=1,n_spin
               do j=n_v_min,n_v
                  do l=1,n_spin
                     do k=1,n_k


               U_cv(i,j,m,l,k)=fxc_0_time(0)*rho_cv(i,j,m,l,k)


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


               U_cc(i,j,m,l,k)=fxc_0_time(0)*rho_cc(i,j,m,l,k)


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

               occup_in=0.0d0
               if(i.eq.j.and.m.eq.l) occup_in=1.0d0


               U_vv(i,j,m,l,k)=fxc_0_time(0)*(rho_vv(i,j,m,l,k)-occup_in)


                     end do
                  end do
               end do
            end do
         end do


return
end
