!this section saves the files

call execute_command_line('cp J/I_omega  J/I_c_omega_10_p_6_w_1.0_t_100_T2_'// adjustl(trim( dirname_real(T_2) ) //'_G_'// adjustl(trim( dirname_real(num3) ))))

!call execute_command_line('cp J/HH_11  J/HH_20_17_t_100_p_6_w_1.0_T2_'// adjustl(trim( dirname_real(T_2) ) //'_G_'// adjustl(trim( dirname_real(num3) ))))

!call execute_command_line('cp J/j_tot_omega  J/j_tot_omega_0_p_6_w_1.0_t_100_T2_'// adjustl(trim( dirname_real(T_2) ) //'_G_'// adjustl(trim( dirname_real(num3) ))))

!call execute_command_line('cp J/j_z J/j_tot_z_0_p_6_w_1.0_t_100_T2_'// adjustl(trim( dirname_real(T_2) )//'_G_'// adjustl(trim( dirname_real(num3) ))))


call execute_command_line('mkdir -p d_J_E_'// adjustl(trim( dirname_real(Efield_max) )))
!call execute_command_line('mkdir -p d_J_E_'//adjustl(trim(dirname_real(Efield_max) )//'_G_'// adjustl(trim(dirname_real(num3) ))))
call execute_command_line('cp -r J/* d_J_E_'// adjustl(trim( dirname_real(Efield_max) )))
!call execute_command_line('cp -r J/* d_J_E_'// adjustl(trim(dirname_real(Efield_max) )//'_G_'// adjustl(trim( dirname_real(num3) ))))
call  execute_command_line('rm -r J')
call  execute_command_line('rm -r rho')
