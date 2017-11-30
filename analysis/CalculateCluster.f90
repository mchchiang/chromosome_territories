      PROGRAM main
      implicit none

      character(40)  filename
      integer N
      integer maxdata
      parameter (maxdata=2097150)
      real RX(maxdata), RY(maxdata), RZ(maxdata)
      integer i,kode
      real rcl, lx,ly,lz
      integer it, nit
      integer nroc, sizeofclusters(maxdata), list(maxdata)

      READ(UNIT=*, FMT=*) filename 
!      write(unit=*,fmt=*)"read file name: ",filename
      read(unit=*, fmt=*)rcl
      read(*,*)lx,ly,lz
      
      !check that the box is cubic
      if (lx .ne. ly .or. lx .ne. lz .or. ly .ne. lz)then
         write(*,*)'only cubic boxes allowed at the moment!'
         stop
      end if
      ! scale the radius to internall units
      rcl = rcl/lx
      
      OPEN(UNIT=15, FILE=filename, STATUS='OLD', IOSTAT=KODE) 
      IF(KODE .NE. 0) THEN 
          WRITE(UNIT=*,FMT=*)filename, ' cannot be opened' 
          stop
      END IF 

      i = 1
      do while (i<= maxdata)
         read(15, *, iostat=kode)RX(i), RY(i), RZ(i)
!         write(unit=*,fmt=*)kode, rx(i), i
         if(kode<0)then
            exit
         else
            i = i + 1
         end if
      end do
      
      N = i-1
!      write(unit=*,fmt=*)N, ' lines read'

      ! place the coordinates in to box of unit length and centered at origin
      do i=1,N
         rx(i) = rx(i) - nint(rx(i)/Lx)*Lx
         rx(i) = rx(i)/Lx

         ry(i) = ry(i) - nint(ry(i)/Ly)*Ly
         ry(i) = ry(i)/Ly

         rz(i) = rz(i) - nint(rz(i)/Lz)*Lz
         rz(i) = rz(i)/Lz
!         write(*,*)rx(i),ry(i),rz(i)
      end do

      it = 1
      call mygang(rcl, it, nit, n, rx, ry, rz,list)
      
     write(*,*)nit,' colloids in cluster containing colloid ',it

     call number_of_clusters(N,list,nroc,sizeofclusters)
     write(*,*)'using cut of radius ',rcl*lx, ' found ',nroc,' clusters'
     do i = 1, nroc
        write(*,*)'cluster ',i,' has ',sizeofclusters(i),' colloids'
     end do
     
    end PROGRAM main

    SUBROUTINE MYGANG ( RCL, IT, NIT, N, RX, RY, RZ, L )

!      COMMON / BLOCK1 / RX, RY, RZ

!    *******************************************************************
!    ** ROUTINE TO IDENTIFY ATOM CLUSTERS IN A CONFIGURATION.         **
!    **                                                               **
!    ** THIS ROUTINE SORTS N ATOMS INTO CLUSTERS DEFINED BY A         **
!    ** CRITICAL CLUSTER RADIUS, AND COUNTS THE NUMBER OF ATOMS IN    **
!    ** THE CLUSTER CONTAINING THE ATOM 'IT'.  THE ATOMS ARE IN A     **
!    ** BOX OF UNIT LENGTH CENTRED AT THE ORIGIN                      **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** STODDARD J COMP PHYS, 27, 291, 1977.                          **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                   NUMBER OF ATOMS                   **
!    ** INTEGER IT                  AN ATOM IN A PARTICULAR CLUSTER   **
!    ** INTEGER L(N)                A LINKED-LIST                     **
!    ** INTEGER NIT                 NUMBER OF ATOMS IN THAT CLUSTER   **
!    ** REAL    RX(N),RY(N),RZ(N)   POSITIONS                         **
!    ** REAL    RCL                 CRITICAL CLUSTER DISTANCE         **
!    *******************************************************************
      implicit none
      INTEGER     N
            
      REAL        RX(N), RY(N), RZ(N)
      REAL        RCL
      INTEGER     IT, NIT
      
      REAL        RCLSQ, RXJK, RYJK, RZJK
      REAL        RJKSQ, RXJ, RYJ, RZJ
      INTEGER     I, J, K, LK, LIT, L(N)
      
!     ****************************************************************

      RCLSQ = RCL * RCL

!     ** SET UP THE SORTING ARRAY **

      DO 10 I = 1, N
         
         L(I) = I
         
 10   CONTINUE

!     ** SORT THE CLUSTERS **
      
      DO 50 I = 1, N - 1
         
         IF ( I .EQ. L(I) ) THEN
            
            J   = I
            RXJ = RX(J)
            RYJ = RY(J)
            RZJ = RZ(J)
            
            DO 20 K = I + 1, N
               
               LK = L(K)
               IF ( LK .EQ. K ) THEN
                  
                  RXJK  = RXJ - RX(K)
                  RYJK  = RYJ - RY(K)
                  RZJK  = RZJ - RZ(K)
                  RXJK  = RXJK - ANINT ( RXJK )
                  RYJK  = RYJK - ANINT ( RYJK )
                  RZJK  = RZJK - ANINT ( RZJK )
                  RJKSQ = RXJK * RXJK + RYJK * RYJK + RZJK * RZJK
                  
                  IF ( RJKSQ .LE. RCLSQ ) THEN
                     
                     L(K) = L(J)
                     L(J) = LK
                     
                  ENDIF

                  
               ENDIF
               
 20         CONTINUE
            
            J   = L(J)
            RXJ = RX(J)
            RYJ = RY(J)
            RZJ = RZ(J)
            
 30         IF ( J .NE. I ) THEN
               
               DO 40 K = I + 1, N
                  
                  LK = L(K)
                  
                  IF ( LK .EQ. K ) THEN
                     
                     RXJK  = RXJ - RX(K)
                     RYJK  = RYJ - RY(K)
                     RZJK  = RZJ - RZ(K)
                     RXJK  = RXJK - ANINT ( RXJK )
                     RYJK  = RYJK - ANINT ( RYJK )
                     RZJK  = RZJK - ANINT ( RZJK )
                     RJKSQ = RXJK * RXJK + RYJK * RYJK + RZJK * RZJK
                     
                     IF ( RJKSQ .LE. RCLSQ ) THEN
                        
                        L(K) = L(J)
                        L(J) = LK
                        
                     ENDIF
                     
                  ENDIF

 40            CONTINUE

               J   = L(J)
               RXJ = RX(J)
               RYJ = RY(J)
               RZJ = RZ(J)
               
               GO TO 30
               
            ENDIF
            
         ENDIF
         
 50   CONTINUE
      
!     **  COUNT THE NUMBER IN A CLUSTER CONTAINING ATOM IT **
      
      NIT = 1
      LIT = L(IT)
 60   IF ( LIT .NE. IT ) THEN
         
         NIT = NIT + 1
         LIT = L(LIT)
         
         GO TO 60
         
      ENDIF
      
      RETURN
      END
      
      subroutine number_of_clusters(N,L,nroc,sizeofclusters)
        implicit none

        integer N
        integer L(N)
        integer nroc

        integer used(N),sizeofclusters(N)
        integer i,lit,nit,it,isused

        used(1:N)=0
        nroc=0
        nit=1
        sizeofclusters(1:N)=0
        do it=1,N
           isused=0
           
           do i = 1,N
              if (it .eq. used(i))then
                 isused = 1
              end if
           end do

           if (isused==0)then
              lit = L(it)
              used(it) = it
              do while(lit .ne. it)
                 used(lit) = lit
                 nit = nit + 1
                 lit = L(lit)
              end do
              nroc = nroc + 1
              sizeofclusters(nroc)=nit
              nit = 1
           end if
        end do

      end subroutine number_of_clusters
