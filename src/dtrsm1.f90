      SUBROUTINE DTRSM4 ( N, M, A, B ,DIAG)
      implicit none
      INTEGER            N, M
      DOUBLE PRECISION   A( N, * ), B( N ), diag( N )
      INTEGER            I, J, K
      DOUBLE PRECISION   TEMP
!
!
!           Form  B := inv( A' )*B.
!
         DO I = 1, M
             TEMP = B( I )
             DO K = 1, I - 1
                TEMP = TEMP - A( K, I )*B( K )
             ENDDO
             TEMP = TEMP * diag( i )
             B( I ) = TEMP
         ENDDO
!
!
!           Form  B := inv( A )*B.
!
         DO K = M, 1, -1
            B( K ) = B( K ) * diag( k )
            DO I = 1, K - 1
               B( I ) = B( I ) - B( K )*A( I, K )
            ENDDO
         ENDDO


      RETURN
!
!     End of DTRSM .
!
      END
