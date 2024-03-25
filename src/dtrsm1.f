      SUBROUTINE DTRSM3 ( M, N, A, B )
      implicit none
      include '../include/sp.inc'
      INTEGER            M, N
      DOUBLE PRECISION   A( nnmod, * ), B( nnmod, * ), diag( nnmod )
      INTEGER            I, INFO, J, K
      DOUBLE PRECISION   TEMP
*
	do i=1,m
         diag(i) = 1.d0/a(i,i)
      enddo
c
c
c$omp parallel share( M, N, A, B, diag),
c$&            local(i,j,k,temp)
c$omp do schedule(runtime)
      DO J = 1, N
*
*           Form  B := inv( A' )*B.
*
         DO I = 1, M
             TEMP = B( I, J )
             DO K = 1, I - 1
                TEMP = TEMP - A( K, I )*B( K, J )
             ENDDO
             TEMP = TEMP * diag( i )
             B( I, J ) = TEMP
         ENDDO
*
*
*           Form  B := inv( A )*B.
*
         DO K = M, 1, -1
            B( K, J ) = B( K, J ) * diag( k )
            DO I = 1, K - 1
               B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
            ENDDO
         ENDDO
c
      ENDDO
*
      RETURN
*
*     End of DTRSM .
*
      END
      SUBROUTINE DTRSM4 ( M, A, B ,DIAG)
      implicit none
      include '../include/sp.inc'
      INTEGER            M
      DOUBLE PRECISION   A( nnmod, * ), B( nnmod ), diag( nnmod )
      INTEGER            I, J, K
      DOUBLE PRECISION   TEMP
*
*
*           Form  B := inv( A' )*B.
*
         DO I = 1, M
             TEMP = B( I )
             DO K = 1, I - 1
                TEMP = TEMP - A( K, I )*B( K )
             ENDDO
             TEMP = TEMP * diag( i )
             B( I ) = TEMP
         ENDDO
*
*
*           Form  B := inv( A )*B.
*
         DO K = M, 1, -1
            B( K ) = B( K ) * diag( k )
            DO I = 1, K - 1
               B( I ) = B( I ) - B( K )*A( I, K )
            ENDDO
         ENDDO
c
*
      RETURN
*
*     End of DTRSM .
*
      END
