module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   real(dp), parameter :: MB_K1 = 1.30d0, MB_K2 = 0.0506d0, MB_M = 0.2177d0

   real(dp) :: mb_B_surf = 0d0
   real(dp) :: mb_lever2 = 1d0
   real(dp) :: mb_torque = 0d0
   real(dp) :: mb_J_removed = 0d0
   real(dp) :: mb_mdot_rot = 0d0
   real(dp) :: mb_mdot_target = 0d0
   real(dp) :: mdot_rot_cgs = 0d0

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      s% other_torque => fossil_braking
      s% other_adjust_mdot => smooth_critical_mdot

      s% extras_startup => extras_startup
      s% extras_start_step => extras_start_step
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns
      s% how_many_extra_history_header_items => how_many_extra_history_header_items
      s% data_for_extra_history_header_items => data_for_extra_history_header_items
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
   end subroutine extras_controls

   subroutine fossil_braking(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: B_ms, R_ms, f_dip, t_B, frac
      real(dp) :: R, M, mdot, omega_s, omega_c, B_dip, v_esc, v_mod
      real(dp) :: lever2, torque, J_region, m_out
      integer :: k, k_in, nz

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      s% extra_jdot(:) = 0d0
      s% extra_omegadot(:) = 0d0
      mb_torque = 0d0
      mb_lever2 = 1d0
      mb_B_surf = 0d0

      if (.not. s% rotation_flag) return

      B_ms = s% x_ctrl(1)
      R_ms = s% x_ctrl(2)
      f_dip = s% x_ctrl(3)
      t_B = s% x_ctrl(4)
      frac = s% x_ctrl(5)
      if (B_ms <= 0d0 .or. R_ms <= 0d0 .or. f_dip <= 0d0) return

      omega_c = star_surface_omega_crit(id, ierr)
      if (ierr /= 0) return

      nz = s% nz
      R = s% r(1)
      M = s% mstar
      mdot = abs(s% mstar_dot)
      omega_s = s% omega_avg_surf

      if (s% x_ctrl(10) >= 0.5d0) mdot = max(mdot - mdot_rot_cgs, 0d0)

      mb_B_surf = B_ms*pow2(R_ms*Rsun/R)
      if (t_B > 0d0) mb_B_surf = mb_B_surf*exp(-s% star_age/t_B)
      B_dip = f_dip*mb_B_surf

      if (mdot <= 0d0 .or. omega_s <= 0d0) return

      v_esc = sqrt(2d0*standard_cgrav*M/R)
      v_mod = sqrt(pow2(MB_K2*v_esc) + pow2(omega_s*R))
      lever2 = MB_K1*MB_K1*pow(pow2(B_dip*R)/(mdot*v_mod), 2d0*MB_M)
      lever2 = max(lever2, 1d0)
      mb_lever2 = lever2

      torque = mdot*omega_s*R*R*(lever2 - 1d0)
      if (torque <= 0d0) return

      k_in = nz
      if (frac > 0d0 .and. frac < 1d0) then
         m_out = 0d0
         do k = 1, nz
            m_out = m_out + s% dm_bar(k)
            if (m_out >= frac*M) then
               k_in = k
               exit
            end if
         end do
      end if

      J_region = dot_product(s% j_rot(1:k_in), s% dm_bar(1:k_in))
      if (J_region <= 0d0) return

      if (s% dt > 0d0) torque = min(torque, 0.1d0*J_region/s% dt)

      do k = 1, k_in
         s% extra_jdot(k) = -torque*s% j_rot(k)/J_region
      end do
      mb_torque = torque

   end subroutine fossil_braking

   subroutine smooth_critical_mdot(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: w, w0, p, fac, f_lag, R, M, L, t_kh, x, mdot_tgt, prev, tau, mdot_rot
      real(dp) :: mdot_floor

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      mb_mdot_rot = 0d0
      mb_mdot_target = 0d0
      mdot_rot_cgs = 0d0
      if (.not. s% rotation_flag) return

      w0 = s% x_ctrl(6)
      p = s% x_ctrl(7)
      fac = s% x_ctrl(8)
      f_lag = s% x_ctrl(9)
      if (w0 <= 0d0 .or. w0 >= 1d0 .or. p <= 0d0 .or. fac <= 0d0) return

      R = s% r(1)
      M = s% mstar
      L = s% L(1)
      if (L <= 0d0 .or. R <= 0d0) return

      t_kh = standard_cgrav*M*M/(R*L)
      mdot_floor = 1d-16*Msun/secyer

      w = s% w_div_w_crit_avg_surf
      mdot_tgt = 0d0
      if (w > w0) then
         x = min((w - w0)/(1d0 - w0), 1d0)
         mdot_tgt = fac*(M/t_kh)*pow(x, p)
      end if

      mdot_rot = mdot_tgt
      if (f_lag > 0d0 .and. s% dt > 0d0) then
         prev = max(s% xtra(1), 0d0)
         tau = f_lag*min(t_kh, max(s% star_age*secyer, s% dt))
         mdot_rot = prev + (mdot_tgt - prev)*(1d0 - exp(-s% dt/tau))
      end if
      if (mdot_rot < mdot_floor) mdot_rot = 0d0

      s% mstar_dot = s% mstar_dot - mdot_rot
      mdot_rot_cgs = mdot_rot
      mb_mdot_rot = mdot_rot/(Msun/secyer)
      mb_mdot_target = mdot_tgt/(Msun/secyer)

   end subroutine smooth_critical_mdot

   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (restart) then
         mb_J_removed = s% xtra(2)
      else
         mb_J_removed = 0d0
         s% xtra(1) = 0d0
         s% xtra(2) = 0d0
      end if
   end subroutine extras_startup

   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_start_step = 0
   end function extras_start_step

   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
   end function extras_check_model

   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 6
   end function how_many_extra_history_columns

   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      names(1) = 'mb_B_surf'
      vals(1) = mb_B_surf
      names(2) = 'mb_lever2'
      vals(2) = mb_lever2
      names(3) = 'mb_torque'
      vals(3) = mb_torque
      names(4) = 'mb_J_removed'
      vals(4) = mb_J_removed
      names(5) = 'mb_mdot_rot'
      vals(5) = mb_mdot_rot
      names(6) = 'mb_mdot_target'
      vals(6) = mb_mdot_target
   end subroutine data_for_extra_history_columns

   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
   end function how_many_extra_profile_columns

   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine data_for_extra_profile_columns

   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 0
   end function how_many_extra_history_header_items

   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine data_for_extra_history_header_items

   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 0
   end function how_many_extra_profile_header_items

   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine data_for_extra_profile_header_items

   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
      mb_J_removed = mb_J_removed + mb_torque*s% dt
      s% xtra(1) = mdot_rot_cgs
      s% xtra(2) = mb_J_removed
   end function extras_finish_step

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine extras_after_evolve

end module run_star_extras
