C *** This is a stand-alone program that solves a simple boundary value
C *** problem as a demonstration of the generalized Henyey solution code
C *** (consisting of subroutines SOLVE, DERIV, and DIFFEQ).  To address
C *** other problems, it is necessary to specify the appropriate difference
C *** equations and boundary conditions in a particular form (see DIFFEQ)
C *** and to modify array sizes (as noted below).  All derivatives are
C *** calculated numerically.  With relatively little work, these routines
C *** can be modified to address a wide range of boundary value problems.
C *** 
C *** This code solves the pair of first-order non-linear differential
C *** equations: dP/dr = - (G M rho)/r^2,   dM/dr = 4 pi r^2 rho
C *** subject to the boundary conditions: M = 4/3 pi r^3 rho (at the
C *** first grid point) and log P = 12.2 (at the final grid point)
C *** using the generalized Henyey code consisting of SOLVE, DERIV, and
C *** DIFFEQ.  For the purposes of this simple exercise, the density (rho)
C *** is assumed to be a function of pressure only:
C ***               rho = 7.463 x 10^(-12) P^(0.76)
C *** A grid consisting of 21 points in log r together with a trial
C *** solution to the problem is given in the DATA statement.
C ***
C *** To run this demonstration fortran code, it must be compiled,
C *** for example as
C ***      gfortran -o henyey demo_henyey.f
C *** and the run as
C ***      ./henyey
C *** Output is sent to unit 6 (monitor) - convergence details and the
C *** final converged "model".
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R(21),P(21),M(21),ERT(2)
      COMMON /MAJOR/ AS(21,2,2),H(21,3),DH(21,3),DLQ,NMESH,NDE,NBC
      DATA R/9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.,10.1,10.2,
     1   10.3,10.4,10.5,10.6,10.65,10.7,10.75,10.8/,P/17.11,17.105,
     2   17.10,17.095,17.09,17.075,17.06,17.03,17.0,16.92,16.85,
     3   16.65,16.5,16.3,15.9,15.2,14.7,14.1,13.7,12.8,12.2/,M/29.5,
     4   29.8,30.1,30.4,30.7,30.9,31.3,31.5,31.8,32.1,32.4,32.6,
     5   32.8,32.95,33.15,33.22,33.3,33.315,33.33,33.335,33.338/
 1000 FORMAT(' ITERATION:',I3,2X,'dP_max =',F8.5,2X,'dM_max =',
     1   F8.5,2X,'AVG. CORRECTION =',F8.5)
 1001 FORMAT(' CONVERGENCE HAS NOT BEEN OBTAINED')
 1002 FORMAT(//,' CONVERGED VALUES ARE',//,3X,'I',7X,'R',9X,'P',9X,'M')
 1003 FORMAT(I4,3F10.4)
C *** NMESH = number of grid points
      NMESH=21
C *** NDE = number of differential (or difference) equations being solved
      NDE=2
C *** NBC = number of boundary conditions at the first grid point
      NBC=1
C *** The following lines define the arrays containing the trial solution
C *** to the problem [in this case, P(I) = log P values and M(I) = log M
C *** values - the dependent variables], as a function of the independent
C *** variable [in this case, R(I) = log r values].  The corrections to
C *** the dependent variables [the DH(I,J) arrays] are initially set to
C *** zero: at convergence, DH(I,1) and DH(I,2) contain the corrections
C *** to the H(I,1) and H(I,2) arrays so that H(I,1)+DH(I,1) give the
C *** log P values that satisfy the difference equations (etc. for log M).
C *** DH(I,3) contains the stepwidths in the independent variable between
C *** adjacent grid points. 
      DO 10 I=1,NMESH
      H(I,1)=P(I)
      H(I,2)=M(I)
      H(I,3)=R(I)
      DH(I,1)=0.D0
      DH(I,2)=0.D0
      IF(I.EQ.1) GO TO 10
      DH(I,3)=R(I)-R(I-1)
   10 CONTINUE
      DH(1,3)=DH(2,3)
      DO 30 NTER=1,10
C *** The generalized Henyey solution routines are called
      CALL SOLVE
C *** This section calculates the maximum and average corrections to
C *** the dependent variables obtained by SOLVE to monitor how the
C *** convergence is proceeding.
      ERR=0.D0
      DO 20 J=1,NDE
      X=0.D0
      DO 15 K=1,NMESH
      Y=AS(K,J,1)
      IF(DABS(Y).LE.DABS(X)) GO TO 15
      X=Y
   15 ERR=ERR+DABS(Y)
      ERT(J)=X
   20 CONTINUE
      ERR=ERR/DFLOAT(NDE*NMESH)
      FAC=DMIN1(0.9D0,0.1D0/ERR)
      WRITE(6,1000) NTER,ERT(1),ERT(2),ERR
      DO 25 K=1,NDE
      DO 25 J=1,NMESH
   25 DH(J,K)=DH(J,K)+FAC*AS(J,K,1)
      IF(ERR.LT.1.D-4) GO TO 35
   30 CONTINUE
      WRITE(6,1001)
      GO TO 45
   35 WRITE(6,1002)
C *** The converged values of the dependent variables are given as a
C *** function of the independent variable, and the execution stops. 
      DO 40 I=1,NMESH
      P(I)=P(I)+DH(I,1)
      M(I)=M(I)+DH(I,2)
      WRITE(6,1003) I,R(I),P(I),M(I)
   40 CONTINUE
   45 STOP
      END
      SUBROUTINE SOLVE                                                  
C *** This subroutine provides the main component of a generalized
C *** Henyey code for solving systems of first-order differential
C *** equations subject to boundary conditions at the first and last
C *** grid points.  All necessary information for the execution of
C *** this subroutine is contained in common block /MAJOR/: the input
C *** quantities are:
C *** NDE = number of difference equations being solved (= number of
C ***       dependent variables
C *** NBC = number of boundary conditions are the first (I=1) grid point
C *** NMESH = number of grid points at which the independent and
C ***       dependent variables are defined
C *** DLQ = the difference in the independent variable between adjacent
C ***       grid points
C *** The only output quantities are the corrections to the J=1,NDE
C *** dependent variables at the I=1,NMESH grid points: these are stored
C *** in the array AS(I,J,1).  The arrays that are referenced by SOLVE
C *** must have at least the following dimensions:
C ***        AS(NMESH,NDE,N7)     AA(N5,NDE)      Y(N5,NDE)
C ***           Y1(N3,NDE)        Y2(N3,NDE)      Z(N5,NDE)
C ***           Z1(N3,NDE)        Z2(N3,NDE)      BC(N3,NBC)
C ***             A(N3)             B(NDE)         BB(NDE)
C ***             W(NDE)            V(NDE)
C *** (The parameters N3, N5, and N7 are defined in the code.)  A DATA
C *** statement defines the "weights" W(I) that multiply the (I+1)st
C *** terms of the difference equations: for most purposes the best
C *** value to use is 0.5.   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*8 AA(4,2),Y(4,2),Z(4,2),Y1(3,2),Y2(3,2),Z1(3,2),Z2(3,2),
     1   BC(3,1),A(3),B(2),BB(2),W(2),V(2)        
      COMMON /MAJOR/ AS(21,2,2),H(21,3),DH(21,3),DLQ,NMESH,NDE,NBC
      DATA W/0.5D0,0.5D0/
      N3=NDE+1                                                          
      N4=NBC+1                                                          
      N5=NDE+NDE-NBC+1                                                  
      N6=NDE-NBC                                                        
      N7=N6+1                                                           
      DO 10 J=1,NDE                                                     
   10 V(J)=1.D0-W(J)                                                    
      I=1                                                               
      IBC=NBC                                                           
      CALL DERIV(Y,Y1,Z1,I,IBC,N3)                                      
      CALL SIMEQ2(Y,N5,N3,NBC)                                          
      DO 15 M=N4,N3                                                     
      DO 15 K=1,NBC                                                     
      BC(M,K)=-Y(M,K)                                                   
   15 Y(M,K)=BC(M,K)                                                    
      IBC=0                                                             
   20 I=I+1                                                             
      IF(I.NE.NMESH) GO TO 25                                           
      IBC=N6                                                            
   25 CALL DERIV(Z,Y2,Z2,I,IBC,N3)                                      
      DO 50 J=1,NDE                                                     
      WT=W(J)*DLQ                                                       
      VT=V(J)*DLQ                                                       
      DO 30 M=1,NDE                                                     
      A(M)=VT*Y1(M,J)+Z1(M,J)                                           
   30 AA(M+N6,J)=WT*Y2(M,J)-Z2(M,J)                                     
      A(N3)=-Z2(N3,J)+WT*Y2(N3,J)+Z1(N3,J)+VT*Y1(N3,J)                  
      DO 40 M=1,N6                                                      
      K=M+NBC                                                           
      DO 35 L=1,NBC                                                     
   35 A(K)=A(K)+Y(K,L)*A(L)                                             
   40 AA(M,J)=A(K)                                                      
      DO 45 L=1,NBC                                                     
   45 A(N3)=A(N3)+Y(N3,L)*A(L)                                          
   50 AA(N5,J)=A(N3)                                                    
      CALL SIMEQ2(AA,N5,N5,NDE)                                         
      DO 60 L=1,N7                                                      
      K=L+NDE                                                           
      DO 55 J=1,NDE                                                     
   55 AS(I,J,L)=-AA(K,J)                                                
      K=L+NBC                                                           
      DO 60 J=1,NBC                                                     
      M=J+N6                                                            
   60 Y(K,J)=AS(I,M,L)                                                  
      DO 65 J=1,N3                                                      
      DO 65 L=1,NDE                                                     
      Y1(J,L)=Y2(J,L)                                                   
   65 Z1(J,L)=Z2(J,L)                                                   
      IF(I.LT.NMESH) GO TO 20                                           
      DO 70 K=N7,NDE                                                    
      DO 70 M=1,N3                                                      
      L=M+N6                                                            
      Z(M,K)=AA(L,K)                                                    
      IF(M.GT.NBC.OR.K.EQ.L) GO TO 70                                   
      Z(M,K)=0.D0                                                       
   70 CONTINUE                                                          
      CALL SIMEQ2(Z,N5,N3,NDE)                                          
      DO 75 L=1,N6                                                      
      K=L+NBC                                                           
   75 B(L)=-Z(N3,K)                                                     
   80 DO 85 J=1,NDE                                                     
      BB(J)=AS(I,J,N7)                                                  
      DO 85 L=1,N6                                                      
   85 BB(J)=BB(J)+AS(I,J,L)*B(L)                                        
      DO 90 J=1,NBC                                                     
      L=J+N6                                                            
   90 AS(I,J,1)=BB(L)                                                   
      DO 95 J=N4,NDE                                                    
      L=J-NBC                                                           
      AS(I,J,1)=B(L)                                                    
   95 B(L)=BB(L)                                                        
      I=I-1                                                             
      IF(I.GT.1) GO TO 80                                               
      DO 100 J=N4,NDE                                                   
      L=J-NBC                                                           
  100 AS(I,J,1)=B(L)                                                    
      DO 110 J=1,NBC                                                    
      VX=BC(N3,J)                                                       
      DO 105 L=1,N6                                                     
      K=L+NBC                                                           
  105 VX=VX+BC(K,J)*B(L)                                                
  110 AS(I,J,1)=VX                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE SIMEQ2(A,NA,NC,NR)                                     
C *** This subroutine implements the Gaussian elimination method of 
C *** diagonalizing a coefficient matrix.  The input array A(NA,J)
C *** is replaced by the diagonalized array (though it should be noted
C *** that the matrix elements above the main diagonal are not reset
C *** to zero value, even though they have been reduced to zero value).
C *** NA gives the number of columns in the coefficient matrix. 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*8 A(NA,1)                                                    
 1000 FORMAT(/,' STOP29: (SIMEQ2)  SINGULAR MATRIX, BIGA =',1P,D12.5)   
      DO 25 J=1,NR                                                      
      J1=J+1                                                            
      BIGA=0.D0                                                         
      DO 10 I=J,NR                                                      
      IF(DABS(BIGA).GE.DABS(A(J,I))) GO TO 10                           
      IA=I                                                              
      BIGA=A(J,I)                                                       
   10 CONTINUE                                                          
      IF(DABS(BIGA).GT.1.D-66) GO TO 15                                 
      WRITE(6,1000) BIGA                                               
      STOP29                                                            
   15 DO 20 K=J,NC                                                      
      SAVE=A(K,J)                                                       
      A(K,J)=A(K,IA)                                                    
      A(K,IA)=SAVE                                                      
   20 A(K,J)=A(K,J)/BIGA                                                
      IF(J.EQ.NR) GO TO 30                                              
      DO 25 I=J1,NR                                                     
      DO 25 K=J1,NC                                                     
   25 A(K,I)=A(K,I)-A(K,J)*A(J,I)                                       
   30 J1=NR+1                                                           
      DO 35 I=2,NR                                                      
      M=NR+2-I                                                          
      IA=M-1                                                            
      DO 35 J=1,IA                                                      
      DO 35 K=J1,NC                                                     
   35 A(K,J)=A(K,J)-A(K,M)*A(M,J)                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE DERIV(B,R,L,K,KBC,N3)                                  
C *** This subroutine calculates the derivatives that are required by
C *** the Henyey solution procedure.
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*8 B(4,2),R(3,2),L(3,2),B1(1),R1(2),L1(2),B2(1),R2(2),L2(2)
      COMMON /MAJOR/ AS(21,2,2),H(21,3),DH(21,3),DLQ,NMESH,NDE,NBC   
      CALL DIFFEQ(B1,R1,L1,0,K)                                       
      DO 10 J=1,NDE                                                     
      R(N3,J)=R1(J)                                                     
   10 L(N3,J)=L1(J)                                                     
      IF(KBC.EQ.0) GO TO 20                                             
      DO 15 J=1,KBC                                                     
   15 B(N3,J)=B1(J)                                                     
   20 DO 35 M=1,NDE                                                     
      VX=DH(K,M)                                                        
      DX=DSIGN(1.D-5*DMAX1(1.D0,DABS(H(K,M)+VX)),VX)                    
      XD=1.D0/DX                                                        
      DH(K,M)=VX+DX                                                     
      CALL DIFFEQ(B2,R2,L2,M,K)                                       
      DO 25 J=1,NDE                                                     
      R(M,J)=XD*(R2(J)-R1(J))                                           
   25 L(M,J)=XD*(L2(J)-L1(J))                                           
      IF(KBC.EQ.0) GO TO 35                                             
      DO 30 J=1,KBC                                                     
   30 B(M,J)=XD*(B2(J)-B1(J))                                           
   35 DH(K,M)=VX                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE DIFFEQ(BVF,RHS,LHS,L,K)
C *** This subroutine specifies the terms in the difference equations that
C *** are solved by the Henyey technique.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 BVF(1),RHS(2),LHS(2),VAR(2)
      COMMON /MAJOR/ AS(21,2,2),H(21,3),DH(21,3),DLQ,NMESH,NDE,NBC
      DATA CX/2.30258509299D0/,PI/3.1415926535D0/
      IF(L.GE.0) GO TO 10
      DO 9 J=1,NDE
    9 VAR(J)=H(K,J)
      GO TO 12
   10 DO 11 J=1,NDE
   11 VAR(J)=H(K,J)+DH(K,J)
   12 PL=VAR(1)
      SL=VAR(2)
      RL=H(K,3)
      P=DEXP(CX*PL)
      S=DEXP(CX*SL)
      R=DEXP(CX*RL)
      RHO=7.463D-12*P**0.76D0
      DMDR=4.D0*PI*RHO*R**3/S
      DPDR=-6.6725985D-8*RHO*S/(R*P)
      IF(L) 100,60,70
   60 DLQ=DH(K,3)
   70 LHS(1)=PL
      RHS(1)=DPDR
      LHS(2)=SL
      RHS(2)=DMDR
C *** The following lines specify the boundary conditions to be used.
      IF(K.NE.1) GO TO 75
      BVF(1)=DMDR-3.D0
   75 IF(K.NE.NMESH) GO TO 100
      BVF(1)=PL-12.2D0
  100 RETURN
      END
