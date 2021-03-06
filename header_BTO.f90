!varaible declaration

integer*8 :: i,i1,i1p,i2,i2p,i3,i3p,i4,iiii, j,j1,j2,j3,k,k1,k2,l,l1,m,m1,m2, p, q, q_bar, a, b, i_c,i_v, i_c1, i_v1,                        &
             i_orb,i_spin,j_spin, l_spin

integer*8, parameter :: n_spin=1, n_time_max=2**17,  n_omega_max=n_time_max/2, n_time_min=-1, n_w_max=1000,      &
                        n_k=288, n_k_full=288, n_v_orb=1, n_c_orb=1,                                              &
                        delta_t_step=1,                                                                          &
                        n_bands=24, n_v=20, n_c=21, n_phi=9,                                                     &
                        n_v_min=n_v-n_v_orb+1, n_c_max=n_c+n_c_orb-1,                                            &
                        n_k_1_full_max=5,n_k_2_full_max=8, n_k_3_full_max=4


real*8, parameter :: time_max=100.0d0, delta_t=time_max/n_time_max,                                &
                     omega_max=0.5d0*1.0d0/delta_t,                                               &
                     delta_omega=omega_max/n_time_max,                                            &
                      U_eff=0.0d0,                                                &
                     tau_pulse=time_max/6.0d0, tau_pulse_c=time_max/2, omega0=1.0d0, alpha_grad=1.0d0, a_lat=7.565d0


real*8, dimension (-n_time_min:n_time_max) ::  Efield, Afield
real*8, dimension (3,-n_time_min:n_time_max) ::  E_vec, A_vec

real*8, dimension (3,-n_time_min:n_time_max) ::  j_tot
real*8, dimension (3,-n_time_min:n_time_max) ::  j_tot_intra_v, j_tot_intra_c,                                &
                                                 j_tot_inter_vv, j_tot_inter_vc,j_tot_inter_cv, j_tot_inter_cc
complex*16, dimension (-n_time_min:n_time_max) ::  current

real*8, dimension (n_v_min:n_v,n_spin,3,-n_time_max:n_time_max) ::  j_intra_v
real*8, dimension (n_c:n_c_max,n_spin,3,-n_time_max:n_time_max) ::  j_intra_c

real*8, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,3,-n_time_min:n_time_max) :: P_inter_vc
real*8, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,3,-n_time_min:n_time_max) :: P_inter_cv
real*8, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,3,-n_time_min:n_time_max) :: P_inter_vv
real*8, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,3,-n_time_min:n_time_max) :: P_inter_cc

real*8, dimension  (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,3,-n_time_min:n_time_max) :: j_inter_vc
real*8, dimension  (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,3,-n_time_min:n_time_max) :: j_inter_cv
real*8, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,3,-n_time_min:n_time_max) :: j_inter_vv
real*8, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,3,-n_time_min:n_time_max) :: j_inter_cc

complex*16, dimension (3,-n_omega_max:n_omega_max) ::  j_tot_omega
complex*16, dimension (-n_time_max/2:n_time_max/2) :: j_omega
real*8, dimension (-100*n_omega_max:100*n_omega_max) :: I_int

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_vc_m1,   rho_vc_m1_1, rho_vc_1_m1, rho_vc_1_m1_1
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_cv_m1, rho_cv_m1_1, rho_cv_1_m1, rho_cv_1_m1_1
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: rho_vv_m1, rho_vv_m1_1, rho_vv_1_m1, rho_vv_1_m1_1
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: rho_cc_m1, rho_cc_m1_1, rho_cc_1_m1, rho_cc_1_m1_1
complex*16, dimension (n_v_min:n_v,n_spin,n_k) :: rho_v_m1, rho_v_m1_1, rho_v_1_m1, rho_v_1_m1_1
complex*16, dimension (n_c:n_c_max,n_spin,n_k) :: rho_c_m1, rho_c_m1_1, rho_c_1_m1, rho_c_1_m1_1


real*8, dimension (-n_time_min:n_time_max) :: density_v, density_c


real*8, dimension (n_k,n_spin,n_bands):: energy
real*8, dimension (n_k_full,n_spin,n_bands):: energy_full
real*8, dimension (-n_k_1_full_max:n_k_1_full_max,-n_k_2_full_max:n_k_2_full_max,-n_k_3_full_max:n_k_3_full_max,n_spin,n_bands) :: energy_full_3D

real*8, dimension (n_k,n_spin,n_v_min:n_v):: energy_v
real*8, dimension (n_k,n_spin,n_c:n_c_max):: energy_c

real*8, dimension (3,n_k):: k_m
real*8, dimension (3,n_k_full):: k_m_full
real*8, dimension (n_k):: weight

real*8, dimension (n_k) :: delta_kx, delta_ky, delta_kz,                     &
                           delta_k, alpha_x, alpha_y, alpha_z

real*8, dimension (n_k,n_spin,n_c:n_c_max,3) :: grad_vel_c
real*8, dimension (n_k,n_spin,n_v_min:n_v,3) :: grad_vel_v
real*8, dimension (-n_k_1_full_max:n_k_1_full_max-1,-n_k_2_full_max:n_k_2_full_max-1,-n_k_3_full_max:n_k_3_full_max-1,n_spin,n_c:n_c_max,3) :: grad_vel_c_full
real*8, dimension (-n_k_1_full_max:n_k_1_full_max-1,-n_k_2_full_max:n_k_2_full_max-1,-n_k_3_full_max:n_k_3_full_max-1,n_spin,n_v_min:n_v,3) :: grad_vel_v_full



complex*16, dimension (n_bands,n_bands,n_spin,n_spin,n_k,3) :: dipole_0,dipole
complex*16, dimension (n_bands,n_bands,n_spin,n_spin,n_k) :: dipole_vc_total, dipole_cv_total,   &
                                                             dipole_vv_total, dipole_cc_total

complex*16, dimension (n_v_min:n_v,n_c:n_c_max,n_spin,n_spin,n_k) :: B_vc, B_vc_1
complex*16, dimension (n_c:n_c_max,n_v_min:n_v,n_spin,n_spin,n_k) :: B_cv, B_cv_1
complex*16, dimension (n_v_min:n_v,n_v_min:n_v,n_spin,n_spin,n_k) :: B_vv, B_vv_1
complex*16, dimension (n_c:n_c_max,n_c:n_c_max,n_spin,n_spin,n_k) :: B_cc, B_cc_1

complex*16, dimension (-n_w_max:n_w_max) :: fxc_0

complex*16, dimension (0:2*n_time_max) :: fxc_0_time
complex*16, dimension(n_v_min:n_c_max,n_v_min:n_c_max,n_spin,n_spin,n_v_min:n_c_max,n_v_min:n_c_max,n_spin,n_spin,n_k,n_k) :: G


!     complex*16, dimension (n_bands,n_bands,n_spin,n_spin,n_k,n_k) :: G


complex*16 :: ii, Ax, Bx, dx,dy,dz, sum_1_1, sum_1_2, sum_2_1, sum_2_2, sum_3_1, sum_3_2, sum_4_1, sum_4_2, grad, grad_1, grad_2,                   &
              F_vc_1, F_vc_2, F_cv_1, F_cv_2, F_cc_1, F_cc_2, F_vv_1, F_vv_2,  F_c_1, F_c_2,  F_v_1, F_v_2

character*100 :: text1,text2
real*8             :: Efield_max, T_2
CHARACTER(100) :: num1char
CHARACTER(100) :: num2char
CHARACTER(100) :: num3char

REAL*8 :: num1
REAL*8 :: num2
REAL*8 :: num3

external sub_B, fast

IF(COMMAND_ARGUMENT_COUNT().NE.3)THEN
  WRITE(*,*)'ERROR, TWO COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  STOP
ENDIF

CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
CALL GET_COMMAND_ARGUMENT(2,num2char)
CALL GET_COMMAND_ARGUMENT(3,num3char)

READ(num1char,*)num1                    !then, convert them to REALs
READ(num2char,*)num2
READ(num3char,*)num3

Efield_max=num1
T_2 = num2
