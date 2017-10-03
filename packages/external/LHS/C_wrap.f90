!   _______________________________________________________________________
!
!   LHS (Latin Hypercube Sampling) wrappers for C clients.
!   Copyright (c) 2006, Sandia National Laboratories.
!   This software is distributed under the GNU Lesser General Public License.
!   For more information, see the README file in the LHS directory.
!
!   NOTE: this "C wrapper layer" is NOT a part of the original LHS source
!   code.  It was added by the DAKOTA team to allow C clients to easily
!   link with the LHS f90 routines without having to assume the burden
!   of managing the "mixed-language string translations" themselves.
!
!   BMA (20160315): Changed to use Fortran 2003 ISO C bindings.  Could
!   go further to process the null-terminated C string here in the
!   Fortran and allow variable-length input.
!   _______________________________________________________________________
!
C These Fortran wrappers circumvent problems with implicit string sizes
C in f90.

C ----- 
C Convert a fixed-length C string to a fixed-length Fortran string
C TODO: 
C  * Could instead iterate until C_NULL_CHAR and allow variable length
C  * Could determine the length of the string from len(fort_str)
C -----
      subroutine lhs_cstr_to_fortran( c_str, fort_str, num_char )

      use iso_c_binding, only: C_CHAR
      integer, intent(in) :: num_char
      character(C_CHAR), intent(in) :: c_str(num_char)
      character fort_str(num_char)

      loop_str: do i=1, num_char
        fort_str(i:i) = c_str(i)
      end do loop_str

      end

C ---------------------------
C Wrapper for LHS lhs_options
C ---------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_options2
      subroutine lhs_options2( lhsreps, lhspval, lhsopts_in, ierror )
     1           bind(C) 

C Fix the string size and always call lhs_options2 from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension (32) :: lhsopts_in

      character*32 lhsopts
      integer      lhsreps, lhspval, ierror

      call lhs_cstr_to_fortran(lhsopts_in, lhsopts, 32)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_options( lhsreps, lhspval, lhsopts, ierror )

      end

C ------------------------
C Wrapper for LHS lhs_dist
C ------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_dist2
      subroutine lhs_dist2( namvar_in, iptflag, ptval, distype_in, 
     1                      aprams, numprms, ierror, idistno, ipvno )
     2           bind(C) 

C Fix the string size and always call lhs_dist2 from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(16) :: namvar_in
      character (kind=C_CHAR, len=1), dimension(32) :: distype_in

      character*16     namvar
      character*32     distype
      integer          iptflag, numprms, ierror, idistno, ipvno
      double precision ptval, aprams(numprms)

      call lhs_cstr_to_fortran(namvar_in, namvar, 16)
      call lhs_cstr_to_fortran(distype_in, distype, 32)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_dist( namvar, iptflag, ptval, distype, aprams,
     1               numprms, ierror, idistno, ipvno )

      end

C -------------------------
C Wrapper for LHS lhs_udist
C -------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_udist2
      subroutine lhs_udist2( namvar_in, iptflag, ptval, distype_in, 
     1                       numpts, xval, yval, ierror, idistno, 
     2                       ipvno )
     3           bind(C) 

C Fix the string size and always call lhs_udist2 from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(16) :: namvar_in
      character (kind=C_CHAR, len=1), dimension(32) :: distype_in

      character*16     namvar
      character*32     distype
      integer          iptflag, numpts, ierror, idistno, ipvno
      double precision ptval, xval(1), yval(1)

      call lhs_cstr_to_fortran(namvar_in, namvar, 16)
      call lhs_cstr_to_fortran(distype_in, distype, 32)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_udist( namvar, iptflag, ptval, distype, numpts,
     1                xval, yval, ierror, idistno, ipvno )

      end

C -------------------------
C Wrapper for LHS lhs_const
C -------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_const2
      subroutine lhs_const2( namvar_in, ptval, ierror, ipvno )
     1           bind(C) 

C Fix the string size and always call lhs_const2 from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(16) :: namvar_in

      character*16     namvar
      integer          ierror, ipvno
      double precision ptval

      call lhs_cstr_to_fortran(namvar_in, namvar, 16)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_const( namvar, ptval, ierror, ipvno )

      end

C ------------------------
C Wrapper for LHS lhs_corr
C ------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_corr2
      subroutine lhs_corr2( nam1_in, nam2_in, corrval, ierror )
     1           bind(C)

C Fix the string size and always call lhs_corr2 from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(16) :: nam1_in, nam2_in

      character*16     nam1, nam2
      integer          ierror
      double precision corrval

      call lhs_cstr_to_fortran(nam1_in, nam1, 16)
      call lhs_cstr_to_fortran(nam2_in, nam2, 16)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_corr( nam1, nam2, corrval, ierror )

      end

C -----------------------
C Wrapper for LHS lhs_run
C -----------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_run2
      subroutine lhs_run2( max_var, max_obs, max_names, ierror, 
     1                     dist_names_in, name_order, pt_vals,
     2                     num_names, sample_matrix, num_vars, 
     3                     rank_matrix, rflag )
     4           bind(C) 

      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(16) :: dist_names_in

      integer          max_var, max_obs, max_names, num_names, num_vars
      integer          rflag, ierror, name_order(1) 
      character*16     dist_names(1) 
      double precision pt_vals(1), sample_matrix(1),rank_matrix(1)

! BMA: there may be an incorrect declaration of this above, or maybe
! the following will fail, but this call isn't currently used
      call lhs_cstr_to_fortran(dist_names_in, dist_names, 16)

      call lhs_run( max_var, max_obs, max_names, ierror, 
     1              dist_names, name_order, pt_vals, num_names,
     2              sample_matrix, num_vars, rank_matrix, rflag )

      end

C -------------------------
C Wrapper for LHS lhs_files
C -------------------------
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::lhs_files2
      subroutine lhs_files2( lhsout_in, lhsmsg_in, lhstitl_in, 
     1                       lhsopts_in, ierror )
     2           bind(C) 

C Fix the string size and always call lhs_files from C++ with strings of
C length 32
      use iso_c_binding, only: C_CHAR
      character (kind=C_CHAR, len=1), dimension(32) :: 
     1           lhsout_in, lhsmsg_in, lhstitl_in, lhsopts_in

      character*32 lhsout, lhsmsg, lhstitl, lhsopts
      integer      ierror

      call lhs_cstr_to_fortran(lhsout_in, lhsout, 32)
      call lhs_cstr_to_fortran(lhsmsg_in, lhsmsg, 32)
      call lhs_cstr_to_fortran(lhstitl_in, lhstitl, 32)
      call lhs_cstr_to_fortran(lhsopts_in, lhsopts, 32)

C Since calling from F90 now, the implicit string size passing should work
      call lhs_files( lhsout, lhsmsg, lhstitl, lhsopts, ierror )
      
      end
