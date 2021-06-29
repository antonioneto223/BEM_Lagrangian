! COPYRIGHT (c) 2007 Council for the Central Laboratory
!               of the Research Councils
! Original date 21 August 2007. Version 1.0.0.
! Purpose: Given a dense unsymmetric n by n matrix, hsl_ma74
! performs partial factorizations and solutions of corresponding sets of
! equations. It is suitable for use in a frontal or multifrontal solver.
! Eliminations are limited to the leading p .le. n columns.
! Can also (partially) solve transpose systems.

! Version 1.5.0
! For version history see ChangeLog

module hsl_ma74_double

  implicit none
  EXTERNAL dger, dgemm, dtrsm, dtrsv, dscal, dswap, dgemv

  integer, parameter, private :: wp = kind(1.0d0) ! Precision parameter.
  real(wp), parameter, private :: one = 1.0_wp
  real(wp), parameter, private :: zero = 0.0_wp

!****************************************************************
! Derived type to control action
      type ma74_control
      integer :: pivoting = 1 !  controls pivoting
!         Possible values are:
!         1  : threshold partial pivoting
!         2  : threshold diagonal pivoting
!         3  : threshold rook pivoting
      real(wp) :: small = 1e-20_wp !
!         If, during the factorization, all the entries in a
!         column of the reduced matrix are of modulus less than
!         or equal to small, all the entries in the column
!         are replaced by zero. Every pivot must also be
!         of absolute value greater than abs(small).
      real(wp) :: static = zero !
!         If static > 0 and fewer than p pivots
!         have been found, choose pivots that do not satisfy threshold
!         test (may replace best pivot candidate by static).
!         Restriction: static = 0.0 or static >= small
      real(wp) :: u = 0.01 ! u is the pivoting threshold
!          parameter. An entry of the leading pxp submatrix is only considered
!          suitable for use as a pivot if it is of absolute value at least as
!          large as u times the entry of largest absolute value in its
!          column (and also its row for rook pivoting).
      end type ma74_control
!****************************************************************
! Derived type for information
      type ma74_info
      integer :: flag = 0 ! Error flag
      integer :: detsign = 1 ! Used to hold the product of the
!        sign of the determinant of  P  and the sign of the
!        determinant of D  and the sign of the
!        determinant of Q, or zero if the  determinant is zero.
      real(wp) :: detlog = zero ! Used to hold the natural logarithm of
!          the absolute value of the determinant of D or zero if the
!          determinant is zero.
      integer :: num_diag = 0 ! number of pivots selected from the diagonal
      integer :: num_nothresh = 0 ! holds number of pivots that do
!          not satisfy threshold condition for the user-supplied u
      integer :: num_perturbed = 0 ! holds number of pivots that were
!          replaced by control%static
       integer :: num_zero = 0 ! holds number zeros on the diagonal of D.
      real(wp) :: spivot = huge(zero) ! on exit holds smallest non-zero pivot
      real(wp) :: usmall = zero ! if num_perturbed = 0, on exit
!          usmall holds threshold parameter that was used and
!          if pivoting = 1 or 2 and static = 0.0 and q < p holds value of u
!          that would have allowed another pivot to be chosen.
      integer :: col_search = 0 ! holds number of cols searched
      integer :: row_search = 0 ! holds number of rows searched
      end type ma74_info
!****************************************************************
contains
      subroutine MA74_factor(n,p,nb,ap,ldap,q,pperm,qperm, &
                 control,info,s)

! Partially factorize a square unsymmetric matrix A with pivoting. The
! factorization has the form
!                 P A Q = (L   ) (D   ) (U  V )
!                         (M  I) (   S) (   I )
! where P and Q are permutation matrices that alters only the first p rows and
! columns, L is unit lower triangular of order q,
! U is unit upper triangular of order q,
! D is diagonal of order q.

      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: p ! number of fully summed rows/cols
      integer, intent(in)  :: nb ! block size used during factorization
      integer, intent(in)  :: ldap ! leading dim. of ap
      real(wp), intent(inout) :: ap(ldap,*)
!        Holds matrix to be factored.
!        On exit, the row/cols that have been factored are overwritten
!        by the factors. On diag., hold the inverse pivot
!        (or zero if pivot is zero). The remainder of the matrix is updated
!        (Schur complement update).
      integer, intent(out)  :: q ! number of eliminations performed
      integer, intent(out)  :: pperm(p) ! pperm(i) holds index of the row of
!        A that is permuted to row i, i = 1,..., p.
      integer, intent(out)  :: qperm(p) ! qperm(j) holds index of the col. of
!        A that is permuted to col j, j = 1,..., p.
      type(ma74_control), intent(in) :: control
      type(ma74_info), intent(out) :: info
      integer, optional, intent(in) :: s ! columns 1:s not suitable for pivoting
!         on until some eliminations performed.

      integer :: idamax
!
      integer :: i
      integer :: j1 ! points to end of previous pivot block
      integer :: jswap
      integer :: kbest ! row index of best pivot candidate
      integer :: kmax ! row index of largest entry in col. m1
      integer :: kpivro ! row index of pivot
      integer :: l
      integer :: lbest ! col. index of best pivot candidate
      integer :: lpivco ! col. index of pivot
      integer :: lcol
      integer :: m
      integer :: m1 ! index of col. currently being tested for pivot
      integer :: m2 ! index of rightmost col. that has been tested since
!                     last major update of the matrix
      integer :: m3 !  no. of rows updated since major update
      integer :: nup ! number of updates applied to all cols (1:p)
!                     (in rook pivoting, the updates have also been applied
!                      to rows (1:p))
      integer :: nzero_col ! no. of cols with all entries .le. control%small
      integer :: p1 ! if p1 < p, cols p1+1:p are zero cols
      integer :: pivoting ! indicates the pivoting being used.
      logical :: chosen ! set to true when pivot found
      logical :: search_row

      real(wp) :: au
      real(wp) :: cmax ! Largest entry in pivot col.
      real(wp) :: pivot
      real(wp) :: ratio ! ratio of largest entry in fully summed part
!                 of col. m1 to largest entry in column m1.
      real(wp) :: ratio_max ! maximum ratio when searching for pivot
      real(wp) :: ratio_c
      real(wp) :: ratio_r
      real(wp) :: rmax ! Largest entry in pivot row.
      real(wp) :: small ! set to abs(control%small).
      real(wp) :: static ! set to control%static.
      real(wp) :: u ! threshold for pivot testing
      real(wp) :: uinitial ! original threshold for pivot testing (u may be
!                 altered if static pivoting in use)

!      real(wp) :: ddot

! Initialise info
      info%flag = 0
      info%num_nothresh = 0
      info%num_perturbed = 0
      info%num_zero = 0
      info%num_diag = 0
      info%col_search = 0
      info%row_search = 0
      info%detsign = 1
      info%detlog = zero
      info%usmall = zero
      info%spivot = huge(zero)

! Ensure q is set
      q = 0
! Check for errors
      if (n < 0) then
        info%flag = -1
      else if (p < 0) then
        info%flag = -2
      else if (p > n) then
        info%flag = -3
      else if (nb < 1) then
        info%flag = -4
      else if (control%static < abs(control%small) .and. &
               control%static /= zero) then
        info%flag = -10
      else if (control%pivoting < 1 .or. control%pivoting > 3) then
        info%flag = -11
      else if (ldap < n) then
        info%flag = -12
      end if
      if (info%flag < 0) return

! nothing we can do if p = 0 or n = 0
      if (p == 0 .or. n == 0) return

! Initialise
      u = min(control%u,one)
      u = max(u,zero)
      uinitial = u
      info%usmall = u
      pivoting = control%pivoting
! ensure small is not equal to zero (otherwise, if u = 0.0
! we could choose a zero pivot, without having checked row
! and col. identically equal to zero)
      small = max(tiny(one),abs(control%small))
      static = control%static
      do l = 1,p
        pperm(l) = l
        qperm(l) = l
      end do

      if (n == 1) then
        pivot = ap(1,1)
        if (abs(pivot) <= small) then
          ap(1,1) = zero
          q = 1
          info%num_zero = 1
        else
          info%detlog = log(abs(pivot))
          ap(1,1) = one/pivot
          q = 1
        end if
        return
      end if

      if (present(s)) then
        if (s > 0 .and. s < p) then
! swap cols 1:min(s,p-s) to end of fully summed block
          do l = 1,min(s,p-s)
! swap col. l with col. p-l+1
            call dswap(n,ap(1,l),1,ap(1,p-l+1),1)
            info%detsign = -info%detsign
            jswap = qperm(l)
            qperm(l) = qperm(p-l+1)
            qperm(p-l+1) = jswap
          end do
        end if
      end if
      m1 = 0
      m2 = 0
      j1 = 0
      nup = 0
      p1 = p
      nzero_col = 0
      q = 0

      kbest = 0; lbest = 0

      if (pivoting == 3) go to 100

! Main loop for threshold partial and diagonal pivoting
  main: do
! Search for pivot among fully summed cols.
          chosen = .false.
          ratio_max = zero
! This loop looks for one pivot
          do m = q+1,p

! m1 controls cyclic search of the candidate columns
! (to avoid repeatedly searching cols that have just been rejected).
            m1 = m1 + 1
            if (m1 > p) then
              m1 = q + 1
              nup = q ! cols 1:p have been updated with q pivots chosen so far
            else if (nup < q) then
! update column m1 before searching it for a pivot candidate.
! First update rows nup+1:q of col m1 using strsv with unit diagonal
              call dtrsv('l','n','u',q-nup,ap(nup+1,nup+1),ldap, &
                   ap(nup+1,m1),1)
! Then update rows q+1:n of col m1 using sgemv
              call dgemv('n',n-q,q-nup,-one,ap(q+1,nup+1),ldap, &
                   ap(nup+1,m1),1,one,ap(q+1,m1),1)
            end if
            m2 = max(m1,m2)  ! m2 does not exceed p

!          write (6,*) 'm1,m2,q,j1,p1,p,qperm',m1,m2,q,j1,p1,p,qperm(m1)
! cycle if col. m1 has been found to be a column of zeros
            if (qperm(m1) < 0) cycle

! In case of diagonal pivoting and u = 0.0, we want to do no
! searching and so we accept diagonal as pivot if it is at least small
! and does not exceed huge.
            if (u == zero .and. pivoting == 2) then
               kpivro = m1
               pivot = ap(kpivro,kpivro)
               if (abs(pivot) <= small) then
                 if (static == zero) then
                   info%flag = -13
                   return
                 end if
! Note: if static > 0 then we will carry on searching
               else if (abs(pivot) > huge(zero)) then
                   info%flag = -14
                   return
               else
! diagonal entry > small so use as a pivot (no comparisons as u = 0.0)
                 chosen = .true.
                 info%num_diag = info%num_diag + 1
                 exit
               end if
            end if

! Find row index of largest entry in col m1
            info%col_search = info%col_search + 1
            kmax = q + idamax(n-q,ap(q+1,m1),1)
            cmax = abs(ap(kmax,m1))
!       write (6,'(a,2es12.4,6i3)') 'cmax,small,q,m1,m2',cmax,small,q,m1,m2
            if (cmax > huge(zero)) then
              info%flag = -14
              return
            end if
! Is every entry in col. m1 too small?
            if (cmax <= small) then
! Set all entries in the col. to zero and swap to end of fully summed block
              ap(q+1:n,m1) = zero
!   write (6,'(a,2es12.4,6i3)') 'cmax,small,q,m1,m2,p1',cmax,small,q,m1,m2,p1
              if (m1 /= p1) then
                call dswap(n,ap(1,m1),1,ap(1,p1),1)
                jswap = qperm(m1)
                qperm(m1) = qperm(p1)
                qperm(p1) = jswap
              end if
! flag col. p1 as being all zeros (so we don't search it again)
              qperm(p1) = -qperm(p1)
              nzero_col = nzero_col + 1
              p1 = p1 - 1
! we want to search the new col. m1 next. To ensure this, decrement m1.
              if (m < p) m1 = m1 - 1
!             m1 = m1 - 1
!             m1 = max(m1,q+1)
              cycle
            end if
            if (pivoting == 1) then
! Jump if pivot found
              if (kmax <= p) then
                kpivro = kmax
                chosen = .true.
                exit
              end if
! Largest entry was not fully-summed.
! Find the largest entry in the fully-summed part of column m1.
              kpivro = q + idamax(p-q,ap(q+1,m1),1)
            else
! diagonal pivoting. Jump if pivot found
              kpivro = m1
              if (kmax == kpivro) then
                chosen = .true.
                exit
              end if
            end if
            pivot = ap(kpivro,m1)
            au = max(cmax*u,small)
!       write (6,'(a,2es12.4,6i3)') 'candidate,cmax',pivot,cmax
            if (abs(pivot) >= au) then
              chosen = .true.
! Check whether pivot satisfies original test
              if (abs(pivot) < max(cmax*uinitial,small)) & 
                  info%num_nothresh = info%num_nothresh + 1
              exit
            else
! Store the row/col indices of the best pivot candidate found during search
              ratio = abs(pivot)/cmax
              if (ratio > ratio_max) then
                ratio_max = ratio
                kbest = kpivro
                lbest = m1
              end if
            end if

          end do

          if (chosen) then
             lpivco = m1
          else if (q == j1) then

! All updates have been applied to pivots chosen so far.
! no further pivots can be chosen that satisfy threshold test
            if (p1 < p) then
! p - p1 cols of all zeros were found. Check to see if
! any of the final p - p1 rows are all zero. If so, a zero
! pivot may be chosen
              do m = p1+1,p
! Search row m for its largest entry.
                info%col_search = info%col_search + 1
                lcol = q + idamax(n-q,ap(m,q+1),ldap)
                rmax = abs(ap(m,lcol))
! Is every entry in row m too small?
                if (rmax <= small) then
! Set all entries in the row and pivot to zero
                  kpivro = m
                  lpivco = p1 + 1
                  ap(kpivro,q+1:n) = zero
                  info%num_zero = info%num_zero + 1
                  info%detlog = zero
                  info%detsign = 0
                  chosen = .true.
                  exit
                end if
              end do
              if (chosen) p1 = p1 + 1
            end if
          end if

          if (.not.chosen .and. q == j1) then
            if (static == zero) then
! we are done (not forcing pivots).
              if (p == n .and. pivoting == 2) then
! we want to pick n pivots but cannot choose them from the diagonal
! so we will switch to off-diagonal pivots
                pivoting = 1
                info%flag = 1
                cycle main
              end if
              info%usmall = ratio_max
              exit main
            else if (lbest > 0) then
! Use the best available pivot (perturb to be at least static)
              kpivro = kbest
              lpivco = lbest
              if (abs(ap(kpivro,lpivco)) < static) then
                ap(kpivro,lpivco) = sign(static,ap(kpivro,lpivco))
                info%num_perturbed = info%num_perturbed + 1
              else
                info%usmall = min(info%usmall,ratio_max)
                u = info%usmall
              end if
              info%num_nothresh = info%num_nothresh + 1
              chosen = .true.
            else
! no pivot on hold.
              if (p == n .and. pivoting == 2) then
! we want to pick n pivots but cannot choose them from the diagonal
! so we will switch to off-diagonal pivots
                pivoting = 1
                info%flag = 1
                cycle main
              end if
! All entries in fully summed part of col. are 0
!!!! we may want to consider perturbing pivots.
              exit main ! we are done as no pivots we can force.
            end if
          end if

        if (chosen) then
! Pivot chosen.
           pivot = ap(kpivro,lpivco)
!          write (6,*) 'pivot',kpivro,lpivco, pivot
           if (lpivco == kpivro) info%num_diag = info%num_diag + 1
           if (pivot /= zero) info%spivot = min(info%spivot,abs(pivot))
           q = q + 1
           if (info%num_zero == 0) then
             if (pivot < zero) info%detsign = -info%detsign
             info%detlog = info%detlog + log(abs(pivot))
           end if

! Swap pivotal row/col with row/col q of ap.
           if (q /= kpivro) then
             call dswap(n,ap(q,1),ldap,ap(kpivro,1),ldap)
             info%detsign = -info%detsign
             jswap = pperm(q)
             pperm(q) = pperm(kpivro)
             pperm(kpivro) = jswap
           end if
           if (q /= lpivco) then
             call dswap(n,ap(1,q),1,ap(1,lpivco),1)
             info%detsign = -info%detsign
             jswap = qperm(q)
             qperm(q) = qperm(lpivco)
             qperm(lpivco) = jswap
           end if

           if (pivot /= zero) then
! Scale pivotal column below diagonal.
             if (q < n) then
               call dscal(n-q,one/pivot,ap(q+1,q),1)
! Update the candidate pivot cols q+1:m1 and rows q+1:n
! (this is so that all cols that have been tested have had the same
! number of updates applied)
               call dger(n-q,m1-q,-one,ap(q+1,q),1,ap(q,q+1),ldap, &
                    ap(q+1,q+1),ldap)
             end if

! If either we have chosen nb pivots since last major update of matrix or we
! have been requested number of pivots, ready for major update
             if (q < p .and. (q/nb)*nb /= q) cycle main

           else
! pivot = zero
             if (p1 == p) then
               exit main
             else
! more zero cols to check
              j1 = q; m1 = q; m2 = q; nup = q; kbest = 0; lbest = 0
              cycle main
            end if
          end if

        end if

!        write (6,*) 'j1,nup,q,m1,m2',j1,nup,q,m1,m2

! Perform update of rest of matrix (called a major update)
       if (nup == j1) then
         if (m1 < n) then
! Update rows nup+1:q and cols m1+1:n (dtrsm)
           call dtrsm('l','l','n','u',q-nup,n-m1,one,ap(nup+1,nup+1),ldap,   &
                ap(nup+1,m1+1),ldap)
! and then update rows q+1 :n  and cols m1+1 :n (gemm)
           if (q < n) &
             call dgemm('n','n',n-q,n-m1,q-nup,-one,                       &
                   ap(q+1,nup+1),ldap,ap(nup+1,m1+1),ldap,one,ap(q+1,m1+1),ldap)
         end if
       else

! Note: if nup > j1 then m2 = p and cols m1+1:p have
!       been updated already for nup pivots.
! Must first update rows nup+1:q of cols m1+1:p for q-nup pivots
          if (m1 < p) then
            call dtrsm('l','l','n','u',q-nup,p-m1,one,ap(nup+1,nup+1),ldap,   &
                 ap(nup+1,m1+1),ldap)
! Then update rows q+1:n of cols m1+1:p (dgemm)
            if (q < n) &
              call dgemm('n','n',n-q,p-m1,q-nup,-one,ap(q+1,nup+1),ldap,      &
                   ap(nup+1,m1+1),ldap,one,ap(q+1,m1+1),ldap)
          end if
          if (p < n) then
! Update rows j1+1:q and cols p+1:n (dtrsm)
            call dtrsm('l','l','n','u',q-j1,n-p,one,ap(j1+1,j1+1),ldap,       &
                 ap(j1+1,p+1),ldap)
! and finally update rows q+1:n  and cols p+1:n (dgemm)
            if (q < n) &
              call dgemm('n','n',n-q,n-p,q-j1,-one,ap(q+1,j1+1),ldap,         &
                   ap(j1+1,p+1),ldap,one,ap(q+1,p+1),ldap)
          end if

        end if

! If q = p pivots already chosen then we have finished
        if (q == p) exit main

! reset ready for next block
        do i = max(1,j1),p
          qperm(i) = abs(qperm(i))
        end do
        p1 = p

        j1 = q;  m1 = q;  m2 = q; nup = q
        kbest = 0; lbest = 0

      end do main

! Store the inverse of non-zero pivots and scale rows that have been pivoted on
      do l = 1,q
        pivot = ap(l,l)
        if (pivot /= zero) then
          ap(l,l) = one/pivot
          if (l < n) call dscal(n-l,one/pivot,ap(l,l+1),ldap)
        end if
      end do

      if (nzero_col > 0) then
        do i = 1,p
          qperm(i) = abs(qperm(i))
        end do
      end if

!          do i = 1,n
!             write (6,'(6es11.4)') ap(i,1:n)
!          end do

      if (p == q) then
        if (info%num_nothresh == 0) then
          info%usmall = control%u
        else if (info%num_perturbed > 0) then
          info%usmall = zero
        end if
      end if

      return  !  threshold partial and diagonal pivoting done

!***************
  100 continue

      m3 = 0
! Main loop for rook pivoting (pivoting =3)

  rook: do
! Search for pivot among fully summed cols.
          chosen = .false.
          ratio_max = zero
! This loop is looking for one pivot.
          do m = q+1,p
! m1 controls cyclic search of the candidate columns
! (to avoid repeatedly searching cols that have just been rejected).
            m1 = m1 + 1
            if (m1 > p) then
              m1 = q + 1
! Entries in cols 1:p have been updated with q pivots chosen so far
! but have to ensure that rows m3+1 :p, cols p+1 :n are also updated
              if (nup < q .and. m3 < p .and. p < n) &
               call dgemm('n','n',p-m3,n-p,q-nup,-one,ap(m3+1,nup+1),    &
                    ldap,ap(nup+1,p+1),ldap,one,ap(m3+1,p+1),ldap)

! reset nup and m3 since columns 1:p and rows 1:p have had q updates
              nup = q
              m3 = q
            else if (nup < q) then
! update column m1 before searching it for a pivot candidate.
! Rows up to row m3 already updated (and m3 >= q)
! so must update rows m3+1,..., n of col m1 using dgemv
              call dgemv('n',n-m3,q-nup,-one,ap(m3+1,nup+1),ldap, &
                   ap(nup+1,m1),1,one,ap(m3+1,m1),1)
            end if

            m2 = max(m1,m2)  ! m2 does not exceed p
  !      write (6,*) 'm,m1,m2,q',m,m1,m2,q
! cycle if col. m1 has been found to be a column of zeros
            if (qperm(m1) < 0) cycle

! Find row index of largest entry in col m1
            info%col_search = info%col_search + 1
            kmax = q + idamax(n-q,ap(q+1,m1),1)
            cmax = abs(ap(kmax,m1))
            if (cmax > huge(zero)) then
              info%flag = -14
              return
            end if
! Is every entry in col. m1 too small?
            if (cmax <= small) then
! Set all entries in the col. to zero and swap to end of fully summed block
              ap(q+1:n,m1) = zero
              if (m1 /= p1) then
                call dswap(n,ap(1,m1),1,ap(1,p1),1)
                jswap = qperm(m1)
                qperm(m1) = qperm(p1)
                qperm(p1) = jswap
              end if
! flag col. p1 as being all zeros (so we don't search it again)
              qperm(p1) = -qperm(p1)
              nzero_col = nzero_col + 1
              p1 = p1 - 1
! we want to search the new col. m1 next. To ensure this, decrement m1.
              if (m < p) m1 = m1 - 1
              cycle
            end if
! Rook pivoting involves searching the candidate pivot row
! If we may have to force pivots then we need the best pivot
! candidate to be on hold and so we must search candidate pivot row.
            search_row = .false.
            if (kmax <= p) then
              kpivro = kmax
              search_row = .true.
            else
! Largest entry was not fully-summed.
! Find the largest entry in the fully-summed part of column m1.
              au = max(cmax*u,small)
              kpivro = q + idamax(p-q,ap(q+1,m1),1)
              pivot = ap(kpivro,m1)
              if (abs(pivot) >= au) search_row = .true.
            end if
            if (static > zero) search_row = .true.
            if (search_row) then
  !       write (6,*) 'kpivro,m3,nup,q,j1',kpivro,m3,nup,q,j1
! must now search kpivro for its largest entry.
! If necessary, must first update row kpivro
              if (kpivro > m3) then
! swap row kpivro with row m3+1 and then update
                if (kpivro /= m3+1) then
                  call dswap(n,ap(m3+1,1),ldap,ap(kpivro,1),ldap)
                  info%detsign = -info%detsign
                  jswap = pperm(m3+1)
                  pperm(m3+1) = pperm(kpivro)
                  pperm(kpivro) = jswap
                  kpivro = m3 + 1
                end if
! update row kpivro columns m1+1:n
! (note: found time very similar using either dgemv or ddot in a loop)
                if (nup < q .and. m1 < n) &
                  call dgemv('t',q-nup,n-m1,-one,ap(nup+1,m1+1),ldap, &
                       ap(kpivro,nup+1),ldap,one,ap(kpivro,m1+1),ldap)
                m3 = m3 + 1
              end if
              info%row_search = info%row_search + 1
              lcol = q + idamax(n-q,ap(kpivro,q+1),ldap)
              rmax = abs(ap(kpivro,lcol))
              if (rmax > huge(zero)) then
                info%flag = -14
                return
              end if
! Is every entry in pivot row too small?
              if (rmax <= small) then
! Set all entries in the row to zero
                ap(kpivro,q+1:n) = zero
                cycle
              end if
              pivot = ap(kpivro,m1)
              if (abs(pivot) >= max(cmax,rmax)*u) then
! threshold condition satisfied so we have a pivot chosen
                chosen = .true.
! Check whether pivot satisfies original test
              if (abs(pivot) < max(cmax,rmax)*uinitial) & 
                  info%num_nothresh = info%num_nothresh + 1
                exit
              end if
! swap the column that had the largest entry in kpivro to be the next 
! one to be searched for a pivot
              if (lcol > m1+1 .and. lcol < p1) then
                call dswap(n,ap(1,m1+1),1,ap(1,lcol),1)
                jswap = qperm(m1+1)
                qperm(m1+1) = qperm(lcol)
                qperm(lcol) = jswap
              end if
            end if
   !     write (6,*) 'rejected,search_row,kpivro,m1,j1,nup,q,m1,m2,m3', &
   !      search_row,kpivro,m1,j1,nup,q,m1,m2,m3
! pivot not chosen.
            if (static > zero) then
! Store the col. index of the best pivot candidate found during search
! (possible row swaps later in search mean no use holding row index).
              ratio = (abs(pivot)/cmax)*(abs(pivot)/rmax)
              if (ratio > ratio_max) then
                ratio_max = ratio
                ratio_c = abs(pivot)/cmax;   ratio_r = abs(pivot)/rmax
                lbest = m1
              end if
            end if

          end do

  !   write (6,*) 'chosen,q,j1,kpivro,m1,m2,m3',chosen,q,j1,kpivro,m1,m2,m3
          if (chosen) then
             lpivco = m1
          else if (q == j1) then
! All updates have been applied to pivots chosen so far.
! no further pivots can be chosen that satisfy threshold test
            if (p1 < p) then
! p - p1 cols of all zeros were found. Check to see if
! any of the final p - p1 rows are all zero. If so, a zero
! pivot may be chosen
              do m = p1+1,p
! Search row m for its largest entry.
                info%row_search = info%row_search + 1
                lcol = q + idamax(n-q,ap(m,q+1),ldap)
                rmax = abs(ap(m,lcol))
! Is every entry in row m too small?
                if (rmax <= small) then
! Set all entries in the row and pivot to zero
                  kpivro = m
                  lpivco = p1 + 1
      !        write (6,*) 'zero pivot',m,p1+1,q
                  ap(kpivro,q+1:n) = zero
                  info%num_zero = info%num_zero + 1
                  info%detlog = zero
                  info%detsign = 0
                  chosen = .true.
                  exit
                end if
              end do
              if (chosen) p1 = p1 + 1
            end if
          end if

          if (.not.chosen .and. q == j1) then
            if (static == zero) then
              exit rook ! we are done (not forcing pivots).
            else if (lbest > 0) then
! Use the best available pivot
              lpivco = lbest
! The row index is the row index of the largest entry in the
! fully-summed part of column lpivco.
              info%col_search = info%col_search + 1
              kpivro = q + idamax(p-q,ap(q+1,lpivco),1)
! Ensure pivot is at least static
              pivot = ap(kpivro,lpivco)
!             write (6,*) 'best ', ap(kpivro,lpivco),static
              if (abs(pivot) < static) then
                ap(kpivro,lpivco) = sign(static,pivot)
                info%num_perturbed = info%num_perturbed + 1
              else
                ratio = min(ratio_c,ratio_r)
                info%usmall = min(info%usmall,ratio)
                u = info%usmall
              end if
              info%num_nothresh = info%num_nothresh + 1
              chosen = .true.
            else
              exit rook ! no pivots we can force.
            end if
          end if

         if (chosen) then
! Pivot chosen.
           pivot = ap(kpivro,lpivco)
           if (lpivco == kpivro) info%num_diag = info%num_diag + 1
           if (pivot /= zero) info%spivot = min(info%spivot,abs(pivot))
           q = q + 1
  !            write (6,*) 'q,kpivro,lpivco,pivot',q,kpivro,lpivco,pivot
           if (info%num_zero == 0) then
             if (pivot < zero) info%detsign = -info%detsign
             info%detlog = info%detlog + log(abs(pivot))
           end if

! Swap pivotal row/col with row/col q of ap.
           if (q /= kpivro) then
             call dswap(n,ap(q,1),ldap,ap(kpivro,1),ldap)
             info%detsign = -info%detsign
             jswap = pperm(q)
             pperm(q) = pperm(kpivro)
             pperm(kpivro) = jswap
           end if
           if (q /= lpivco) then
             call dswap(n,ap(1,q),1,ap(1,lpivco),1)
             info%detsign = -info%detsign
             jswap = qperm(q)
             qperm(q) = qperm(lpivco)
             qperm(lpivco) = jswap
           end if
   !           write (6,*) q,kpivro,lpivco,pperm(kpivro),qperm(lpivco),pivot

           if (pivot /= zero) then
             if (q < n) then
! Scale pivotal column below diagonal.
               call dscal(n-q,one/pivot,ap(q+1,q),1)
! Update the candidate pivot cols q+1:m1, rows q+1: n
               call dger(n-q,m1-q,-one,ap(q+1,q),1,ap(q,q+1),ldap, &
                    ap(q+1,q+1),ldap)
             end if
             if (q < m3 .and. m1 < n) &
! Also update rows q+1:m3, cols m1+1:n
             call dger(m3-q,n-m1,-one,ap(q+1,q),1,ap(q,m1+1),ldap,  &
                  ap(q+1,m1+1),ldap)

! If either we have chosen nb pivots since last major update of matrix or we
! have been chosen sufficient pivots, ready for major update
             if (q < p .and. (q/nb)*nb /= q) cycle rook
           else
! pivot = zero
             if (p1 == p) then
               exit rook
             else
! more zero cols to check
              j1 = q; m1 = q; m2 = q; m3 = q; nup = q
              kbest = 0; lbest = 0
              cycle rook
            end if
          end if

        end if

! Update rest of matrix (call this a major update)
!         write (6,*) 'j1,nup,q,m1,m2,m3',j1,nup,q,m1,m2,m3
        if (nup == j1) then
! columns m1+1:m2 have to be updated (they may have been partially updated
! since the last major update but rows m3+1 :n
! of these columns must be updated)
          if (m1 < m2) &
            call dgemm('n','n',n-m3,m2-m1,q-j1,-one,ap(m3+1,j1+1),ldap,    &
                 ap(j1+1,m1+1),ldap,one,ap(m3+1,m1+1),ldap)
! Update rows m3+1 :n  and cols m2+1 :n (gemm)
          if (m2 < n) &
            call dgemm('n','n',n-m3,n-m2,q-j1,-one,ap(m3+1,j1+1),ldap, &
                 ap(j1+1,m2+1),ldap,one,ap(m3+1,m2+1),ldap)
        else

! nup > j1. Must have m2=p at this point.
! nup updates have been applied to rows 1:p and cols 1:p.
! Update rows m3+1: p (=m2) so that q updates have been applied to
! rows 1:p and to cols 1:p
          if (nup < q) then
            if (m1 < m2) &
              call dgemm('n','n',n-m3,m2-m1,q-nup,-one,ap(m3+1,nup+1),ldap,  &
                   ap(nup+1,m1+1),ldap,one,ap(m3+1,m1+1),ldap)
! Then update rows m3+1:m2 and cols m2+1:n with q-nup pivots
            if (m2 < n .and. m3 < m2) &
              call dgemm('n','n',m2-m3,n-m2,q-nup,-one,ap(m3+1,nup+1),ldap,   &
                   ap(nup+1,m2+1),ldap,one,ap(m3+1,m2+1),ldap)
          end if

! Finally, update rows m2+1 :n  and cols m2+1 :n (gemm)
          if (m2 < n) &
            call dgemm('n','n',n-m2,n-m2,q-j1,-one,ap(m2+1,j1+1),ldap,       &
                 ap(j1+1,m2+1),ldap,one,ap(m2+1,m2+1),ldap)

        end if



! If q = p pivots already chosen then we have finished
        if (q == p) exit rook

! reset ready for next block
        do i = max(1,j1),p
          qperm(i) = abs(qperm(i))
        end do
        p1 = p;  j1 = q;  m1 = q;  m2 = q;  m3 = q;  nup = q
        kbest = 0; lbest = 0

!           do i = 1,n
!              write (6,'(6es11.4)') ap(i,1:n)
!           end do

      end do rook

! Store the inverse of non-zero pivots and scale rows that have been pivoted on
      do l = 1,q
        pivot = ap(l,l)
        if (pivot /= zero) then
          ap(l,l) = one/pivot
          if (l < n) call dscal(n-l,one/pivot,ap(l,l+1),ldap)
        end if
      end do

      if (nzero_col > 0) then
        do i = 1,p
          qperm(i) = abs(qperm(i))
        end do
      end if

      if (p == q) then
        if (info%num_nothresh == 0) then
          info%usmall = control%u
        else if (info%num_perturbed > 0) then
          info%usmall = zero
        end if
      end if

      end subroutine ma74_factor

!*************************************************
      subroutine MA74_solveL1(n,q,nb,ap,ldap,b,flag)
! Partial forward substitution with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return


      do j = 1,q,nb
! All blocks except possibly the final one are order nb
        jb = min(q-j+1,nb)
        call dtrsv('l','n','u',jb,ap(j,j),ldap,b(j),1)
        if (n >= j+jb) &
          call dgemv('n',n-jb-j+1,jb,-one,ap(j+jb,j),ldap,b(j),1,one,b(j+jb),1)
      end do

      end subroutine ma74_solveL1
!*************************************************
      subroutine MA74_solveD1(n,q,nb,ap,ldap,b,flag)
! Partial diagonal solve with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = 1,q
        b(j) = b(j)*ap(j,j)
      end do

      end subroutine MA74_solveD1
!*************************************************
      subroutine MA74_solveDU1(n,q,nb,ap,ldap,b,flag)
! Partial diagonal + back substitution with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = 1,q
        b(j) = b(j)*ap(j,j)

      end do

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemv('n',jb,n-j,-one,ap(j1,j+1),ldap,b(j+1),1,one,b(j1),1)
        call dtrsv('u','n','u',jb,ap(j1,j1),ldap,b(j1),1)
      end do

      end subroutine MA74_solveDU1
!*************************************************
      subroutine MA74_solveU1(n,q,nb,ap,ldap,b,flag)
! Partial back substitution with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemv('n',jb,n-j,-one,ap(j1,j+1),ldap,b(j+1),1,one,b(j1),1)
        call dtrsv('u','n','u',jb,ap(j1,j1),ldap,b(j1),1)
      end do

      end subroutine MA74_solveU1
!*************************************************
      subroutine MA74_solveL2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial forward substitution with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveL1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = 1,q,nb
! All blocks except possibly the final one are order nb
        jb = min(q-j+1,nb)
        call dtrsm('l','l','n','u',jb,nrhs,one,ap(j,j),ldap,b(j,1),ldb)
        if (n >= j+jb) &
          call dgemm('n','n',n-jb-j+1,nrhs,jb,-one,ap(j+jb,j),ldap, &
                      b(j,1),ldb,one,b(j+jb,1),ldb)
      end do

      end subroutine ma74_solveL2
!*************************************************
      subroutine MA74_solveD2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial diagonal solve with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      real(wp) :: atemp

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      do j = 1,q
        atemp = ap(j,j)
        b(j,1:nrhs) = b(j,1:nrhs)*atemp
      end do

      end subroutine MA74_solveD2
!*************************************************
      subroutine MA74_solveDU2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial diagonal + back substitution with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1
      real(wp) :: atemp

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveDU1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = 1,q
        atemp = ap(j,j)
        b(j,1:nrhs) = b(j,1:nrhs)*atemp
      end do

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemm('n','n',jb,nrhs,n-j,-one,ap(j1,j+1),ldap, &
                      b(j+1,1),ldb,one,b(j1,1),ldb)
        call dtrsm('l','u','n','u',jb,nrhs,one,ap(j1,j1),ldap,b(j1,1),ldb)
      end do

      end subroutine MA74_solveDU2
!*************************************************
      subroutine MA74_solveU2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial back substitution with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveU1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemm('n','n',jb,nrhs,n-j,-one,ap(j1,j+1),ldap, &
                      b(j+1,1),ldb,one,b(j1,1),ldb)
        call dtrsm('l','u','n','u',jb,nrhs,one,ap(j1,j1),ldap,b(j1,1),ldb)
      end do

      end subroutine MA74_solveU2

!*************************************************
      subroutine MA74_solveUT1(n,q,nb,ap,ldap,b,flag)
! Partial forward substitution for transpose system with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = 1,q,nb
! All blocks except possibly the final one are order nb
        jb = min(q-j+1,nb)
        call dtrsv('u','t','u',jb,ap(j,j),ldap,b(j),1)
        if (n >= j+jb) &
          call dgemv('t',jb,n-jb-j+1,-one,ap(j,j+jb),ldap,b(j),1,one,b(j+jb),1)
      end do

      end subroutine ma74_solveUT1

!*************************************************
      subroutine MA74_solveDLT1(n,q,nb,ap,ldap,b,flag)
! Partial diagonal + back substitution for transpose system
! with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = 1,q
        b(j) = b(j)*ap(j,j)
      end do

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemv('t',n-j,jb,-one,ap(j+1,j1),ldap,b(j+1),1,one,b(j1),1)
        call dtrsv('l','t','u',jb,ap(j1,j1),ldap,b(j1),1)
      end do

      end subroutine MA74_solveDLT1
!*************************************************
      subroutine MA74_solveLT1(n,q,nb,ap,ldap,b,flag)
! Partial back substitution for transpose system with single right-hand side
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      real(wp), intent(inout), dimension(:) :: b(n)
!        on entry, holds right-hand side. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0) return

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemv('t',n-j,jb,-one,ap(j+1,j1),ldap,b(j+1),1,one,b(j1),1)
        call dtrsv('l','t','u',jb,ap(j1,j1),ldap,b(j1),1)
      end do

      end subroutine MA74_solveLT1
!*************************************************
      subroutine MA74_solveUT2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial forward substitution for transpose system
! with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveUT1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = 1,q,nb
! All blocks except possibly the final one are order nb
        jb = min(q-j+1,nb)
        call dtrsm('l','u','t','u',jb,nrhs,one,ap(j,j),ldap,b(j,1),ldb)
        if (n >= j+jb) &
          call dgemm('t','n',n-jb-j+1,nrhs,jb,-one,ap(j,j+jb),ldap, &
                      b(j,1),ldb,one,b(j+jb,1),ldb)
      end do

      end subroutine ma74_solveUT2

!*************************************************
      subroutine MA74_solveDLT2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial diagonal + back substitution for transpose system
! with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1
      real(wp) :: atemp

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveDLT1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = 1,q
        atemp = ap(j,j)
        b(j,1:nrhs) = b(j,1:nrhs)*atemp
      end do

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemm('t','n',jb,nrhs,n-j,-one,ap(j+1,j1),ldap, &
                      b(j+1,1),ldb,one,b(j1,1),ldb)
        call dtrsm('l','l','t','u',jb,nrhs,one,ap(j1,j1),ldap,b(j1,1),ldb)
      end do

      end subroutine MA74_solveDLT2
!*************************************************
      subroutine MA74_solveLT2(n,q,nb,nrhs,ap,ldap,b,ldb,flag)
! Partial back substitution for transpose system
! with multiple right-hand sides
      integer, intent(in)  :: n ! order of matrix
      integer, intent(in)  :: q ! number of pivots
      integer, intent(in)  :: nb ! block size
      integer, intent(in)  :: nrhs ! no. of right-hand sides
      integer, intent(in)  :: ldap ! leading dimension ap
      real(wp), intent(in) :: ap(ldap,*)  ! Must hold factors.
      integer, intent(in)  :: ldb ! leading dimension of b
      real(wp), intent(inout), dimension(:) :: b(ldb,nrhs)
!        on entry, holds right-hand sides. Overwritten by solution.
      integer, intent(out)  :: flag ! error flag

      integer :: j ! Column index
      integer :: jb ! Number of columns in current block column of ap
      integer :: j1

! Check for errors
      flag = 0
      if (n < 0) then
        flag = -1
      else if (nrhs < 0) then
        flag = -5
      else if (ldb < n) then
        flag = -6
      else if (q < 0) then
        flag = -8
      else if (q > n) then
        flag = -9
      else if (nb < 1) then
        flag = -4
      else if (ldap < n) then
        flag = -12
      end if
      if (flag < 0 .or. n == 0 .or. q == 0 .or. nrhs == 0) return

      if (nrhs < 4) then
        do j = 1, nrhs
          call MA74_solveLT1(n,q,nb,ap,ldap,b(1,j),flag)
        end do
        return
      end if

      do j = q,1,-nb
! All blocks except possibly the final one are order nb
        jb = min(j,nb)
        j1 = j - jb + 1
        if (j < n) &
          call dgemm('t','n',jb,nrhs,n-j,-one,ap(j+1,j1),ldap, &
                      b(j+1,1),ldb,one,b(j1,1),ldb)
        call dtrsm('l','l','t','u',jb,nrhs,one,ap(j1,j1),ldap,b(j1,1),ldb)
      end do

      end subroutine MA74_solveLT2

end module hsl_ma74_double
