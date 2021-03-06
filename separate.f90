     program separate

       implicit none
        integer :: i

        real*8, dimension (4,20):: k_m

        open(1,File='d_J_E_0.5/HH_20_17_t_100_p_6_w_1.0_T2_10.0_G_0.0')
        do i=1,20
           read(1,*) k_m(1,i), k_m(2,i), k_m(3,i), k_m(4,i)
        end do
        close(1)
        open(2,file='HH_6')

        write(2,'(10(2I4))')  k_m
       close(2)
 

        end
