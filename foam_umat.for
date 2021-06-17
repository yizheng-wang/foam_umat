      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP, 
     2 PREDEF,DPRED,CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,  
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C    
      INCLUDE 'ABA_PARAM.INC' 
C   
      CHARACTER*8 CMNAME 
C     
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),    
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),     
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),     
     3 DFGRD0(3, 3), DFGRD1(3, 3),STRANPL(NTENS),S(NTENS),DSTRANPL(NTENS),
     4 TIME(2),DE(NTENS, NTENS), DA(NTENS, NTENS), SA(NTENS, NTENS),
     5 IS(NTENS), JS(NTENS), DE_INV(NTENS, NTENS), TT(NTENS),
     6  DT1DSIG(NTENS), DT2DSIG(NTENS),
     7   DJTDSIG(NTENS), S2(NTENS),
     8 F_V(NTENS), S_V(NTENS), DA_R(NTENS, NTENS)
          REAL*8 M, FAI, DFAIDZ,DA1, DA2,DFAIDJT,T1, T2, EN, JT,lambda,NU
C   
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0,FOUR=4.D0, 
     1SIX=6.D0,ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6, XIAO=1.0D-10) 

C ---------------------------------------------------------------
C    PROPS(1) - E(弹性模量)
C    PROPS(2) - NU(泊松比)
C    PROPS(3) - D0
C    PROPS(4) - m
C    PROPS(5) - Z10
C    PROPS(6) - H
C    PROPS(7) - Z2S
C    PROPS(8) - EN lambda的敏感性参数
C    PROPS(9) - ALPHA
C    PROPS(10) - Z20
C ---------------------------------------------------------------
C    STATEV(1) - Z1
C    STATEV(2) - Z2
C ---------------------------------------------------------------C

!  初始化，时间步等于0的时候，对状态变量初始化
! 时间步等于0的时候，历史变量赋予值，这里或许有初始化的子程序,
! 但是这里我们还是用简单的


        E = PROPS(1)
        NU = PROPS(2)
        D0 = PROPS(3)
        M = PROPS(4)
        Z10 = PROPS(5)
        H = PROPS(6)
        Z2S =PROPS(7)
        EN = PROPS(8)
        ALPHA = PROPS(9)
        if (time(2).eq.0) then 
            STATEV(1) = Z10
            STATEV(2) = PROPS(10) ! 论文中没有给，所以我用0作为初始值
        end if
        Z1 = STATEV(1) ! 赋予历史变量的含义
        Z2 = STATEV(2)
        Z = Z1 + Z2
! 赋予DDSDDE，弹性试应力
! 首先赋予弹性矩阵
      EMOD=PROPS(1)
      ENU=PROPS(2)      
      EBULK3=EMOD/(ONE-TWO*ENU)      
      EG2=EMOD/(ONE+ENU)      
      EG=EG2/TWO      
      EG3=THREE*EG      
      ELAM=(EBULK3-EG2)/THREE            
	  
C---------------------------------------------------------------C 
C    刚度
C---------------------------------------------------------------C  

      DE = ZEROS
      DO K1=1, NDI
        DO K2=1, NDI
          DE(K2, K1)=ELAM
        END DO
        DE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
        DE(K1 ,K1)=EG
      END DO
! 利用弹性矩阵，得到弹性试应力
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K1)=STRESS(K1)+DE(K1, K2)*DSTRAN(K2)
        END DO
      END DO
! 计算得lambda，这里lambda是率形式，dlambda是增量形式。
! 要计算lambda首先要计算JT以及Z
      SHYDRO=ZERO
      DO K1=1, NDI
      SHYDRO=SHYDRO+STRESS(K1)  
      END DO
      SHYDRO=SHYDRO/THREE
      DO K1=1, NDI
          S(K1)=STRESS(K1)-SHYDRO
      END DO
      DO K1=NDI+1, NTENS
        S(K1)=STRESS(K1)
      END DO
      S2 = S
      DO K1 = 4,6 ! 定义工程应力，为了后面方便
        S2(K1) = TWO * S(K1)
      END DO
      JT=ZERO
      DO K1=1, NDI
          JT=JT+(ONE/TWO)*S(K1)*S(K1) 
      END DO
      DO K1=NDI+1, NTENS
          JT=JT+S(K1)*S(K1)
      END DO
      IF (JT.LT.XIAO) THEN
        lambda = ZERO
      else
        lambda = D0/SQRT(JT)*(THREE*JT/(Z)**TWO)**EN ! 获得lambda的率形式
      ENDIF
      
      dlambda = lambda*DTIME ! 乘以时间获得lambda的增量
! 将弹性试应力拉回
      ! 这里dlambda乘以S得到的是塑性真应变，
      !但是与弹性矩阵相乘是工程应变，所以这里分为主与切部分
      DO K1=1, NDI  
        DO K2=1, NTENS
          STRESS(K1)=STRESS(K1)-dlambda*DE(K1, K2)*S(K2)
       END DO
      END DO

      DO K1=NDI+1, NTENS  
        DO K2=1, NTENS
          STRESS(K1)=STRESS(K1)-TWO*dlambda*DE(K1, K2)*S(K2)
       END DO
      END DO
! 应力更新完成，最后更新历史变量以及输出DDSDDE
      Z1 = Z1 + dlambda*M*((Z1-(ONE-ALPHA)*Z10)/Z10)*TWO*JT
      Z2 = Z2 - dlambda*H*(ONE - Z2/Z2S)*TWO*JT
      ! 更新完成后，将历史变量存起来
      STATEV(1) = Z1
      STATEV(2) = Z2
! 程序主体部分完成，下面输出DDSDDE
! 第一步需要求出DA
! 求DA的关键是对弹性矩阵求逆以及得到SA,SA是偏应力对应力求导，是一个四阶张量
! FORTRAN DATA是从列开始赋值

        SA(1,1) = TWO/THREE
        SA(2,1) = -ONE/THREE
        SA(3,1) = -ONE/THREE
        SA(4,1) = ZERO
        SA(5,1) = ZERO
        SA(6,1) = ZERO

        SA(1,2) =  -ONE/THREE
        SA(2,2) =  TWO/THREE
        SA(3,2) = -ONE/THREE
        SA(4,2) = ZERO
        SA(5,2) = ZERO
        SA(6,2) = ZERO

        SA(1,3) = -ONE/THREE
        SA(2,3) = -ONE/THREE
        SA(3,3) = TWO/THREE
        SA(4,3) = ZERO
        SA(5,3) = ZERO
        SA(6,3) = ZERO

        SA(1,4) = ZERO
        SA(2,4) = ZERO
        SA(3,4) = ZERO
        SA(4,4) = ONE
        SA(5,4) = ZERO
        SA(6,4) = ZERO

        SA(1,5) = ZERO
        SA(2,5) = -ZERO
        SA(3,5) = ZERO
        SA(4,5) = ZERO
        SA(5,5) = ONE
        SA(6,5) = ZERO

        SA(1,6) = ZERO
        SA(2,6) = ZERO
        SA(3,6) = ZERO
        SA(4,6) = ZERO
        SA(5,6) = ZERO
        SA(6,6) = ONE

        DE_INV = DE
        N = 6
        CALL BRINV(DE_INV,N,L,IS,JS) ! 对弹性矩阵求逆
        DA = DE_INV + dlambda*SA
        CALL BRINV(DA,N,L,IS,JS) ! 得到DA
! 接下来求TT， 用于求解DDSDDE，是求解DDSDDE的一部分,TT是三叉树符号
! 下面对求TT做好预备工作
        FAI = lambda ! 文中的fai就是lambda
        DFAIDZ = -EN/Z*FAI
        SNORM = TWO*JT
        DA1 = ONE/((ONE/M)+dlambda*SNORM)
        DT1DSIG = TWO*(Z1-(1-ALPHA)*Z10)/Z10*S
        DA2 = ONE/((ONE/H)+dlambda*SNORM)
        DT2DSIG = -TWO*(ONE-Z2/Z2S)*S
        IF (JT.EQ.0) THEN
          DFAIDJT = 0D0
        ELSE
          DFAIDJT = ONE/JT*(EN-ONE/TWO)*FAI
        endif
        DJTDSIG = S
        T1 = -(Z1-(ONE-ALPHA)*Z10)/Z10*2*JT
        T2 = -(ONE-Z2/Z2S)*TWO*JT
! 求TT,由于TT是对称张量，为了方便，我们将TT储存为类似工程应变的数据结构
        TT = -dlambda*DFAIDZ*DA1*DT1DSIG-dlambda*DFAIDZ*DA2*DT2DSIG
     1 +DFAIDJT*DJTDSIG
        ! 将TT的4到6，都包含，放大一倍
        DO K1 = 4,6
            TT(K1) = TWO * TT(K1)
        END DO
        THETA = ZERO
        DO I = 1, NTENS
            DO J = 1, NTENS
                THETA = THETA + TT(I)*DA(I, J)*S2(J)
            ENDDO
        ENDDO 
        THETA = THETA + DFAIDZ*DA1*T1 + DFAIDZ*DA2*T2
        OMIGA = ONE/DTIME+THETA
! 接下来求DDSDDE的后半部分
        F_V = ZERO
        DO K1=1,6
            DO K2=1,6
                F_V(K1) = F_V(K1)+DA(K1, K2) * S2(K2)
            ENDDO
        ENDDO

       S_V = ZERO
        DO K1=1,6
            DO K2=1,6
                S_V(K1) = S_V(K1)+DA(K1, K2) * TT(K2)
            ENDDO
        ENDDO

        DO I=1,6
            DO J=1,6
                DA_R(I, J) = ONE/OMIGA*F_V(I)*S_V(J)
            ENDDO
        ENDDO
        DDSDDE = DA-DA_R
c        DDSDDE = DE
c        if (MOD(TIME(2), 1.0).EQ.0) THEN
        if (NOEL.EQ.455) THEN
          if (NPT.EQ.1) THEN
            WRITE(*, *) 'TIME'
            WRITE(*, *) TIME(2)
            WRITE(*, *) 'STRESS1'
            WRITE(*, *) STRESS(1)
            WRITE(*, *) 'STRESS2'
            WRITE(*, *) STRESS(2)
            WRITE(*, *) 'STRESS3'
            WRITE(*, *) STRESS(3)
            WRITE(*, *) 'LAMBDA'
            WRITE(*, *) lambda
            WRITE(*, *) 'JT'
            WRITE(*, *) JT
            WRITE(*, *) 'Z1'
            WRITE(*, *) Z1
            WRITE(*, *) 'Z2'
            WRITE(*, *) Z2
            WRITE(*, *) 'Z'
            WRITE(*, *) Z
          ENDIF
        END IF
      RETURN
      END



c   这是求逆的子程序，A是需要求逆的矩阵，N是矩阵的阶数，L不用管，
c   IS与JS都是与矩阵阶数相同的一维数组
c   程序引用自《fortran常用算法程序集》徐士良，第二章，第三节，实矩阵求逆
        SUBROUTINE BRINV(A,N,L,IS,JS)
        DIMENSION A(N,N),IS(N),JS(N)
        DOUBLE PRECISION A,T,D
        L=1
        DO 100 K=1,N
        D=0.0
        DO 10 I=K,N
        DO 10 J=K,N
            IF (ABS(A(I,J)).GT.D) THEN
            D=ABS(A(I,J))
            IS(K)=I
            JS(K)=J
            END IF
10	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
c	    WRITE(*,20)
	    RETURN
	  END IF
20	  FORMAT(1X,'ERR**NOT INV')
	  DO 30 J=1,N
	    T=A(K,J)
	    A(K,J)=A(IS(K),J)
	    A(IS(K),J)=T
30	  CONTINUE
	  DO 40 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,JS(K))
	    A(I,JS(K))=T
40	  CONTINUE
	  A(K,K)=1/A(K,K)
	  DO 50 J=1,N
	    IF (J.NE.K) THEN
	      A(K,J)=A(K,J)*A(K,K)
	    END IF
50	  CONTINUE
	  DO 70 I=1,N
	    IF (I.NE.K) THEN
	      DO 60 J=1,N
	        IF (J.NE.K) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
	        END IF
60	      CONTINUE
	    END IF
70	  CONTINUE
	  DO 80 I=1,N
	    IF (I.NE.K) THEN
	      A(I,K)=-A(I,K)*A(K,K)
	    END IF
80	  CONTINUE
100	    CONTINUE
        DO 130 K=N,1,-1
        DO 110 J=1,N
            T=A(K,J)
            A(K,J)=A(JS(K),J)
            A(JS(K),J)=T
110	  CONTINUE
	  DO 120 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,IS(K))
	    A(I,IS(K))=T
120	  CONTINUE
130	    CONTINUE
	    RETURN
	    END      