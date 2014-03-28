C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  N                 the number of points being clustered    C
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C                    MEMBR must be initialized by R to the   C
C                    default of  rep(1, N)                   C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C  Modifications for R: Ross Ihaka, Dec 1996                 C
C                       Fritz Leisch, Jun 2000               C
C  all vars declared:   Martin Maechler, Apr 2001            C
C  PR#4195 fixed by BDR Nov 2003                             C
C  The code tried to update DISTNN(I2) on the fly,           C
C  but since it recalculated all the others it might as      C
C  well recalculate that one too.                            C
C------------------------------------------------------------C
      SUBROUTINE HCLUST2(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,
     X                  FLAG,DISS,ALPHA,BETA)
c Args
      INTEGER N, LEN, IOPT
      INTEGER IA(N),IB(N), NN(N)
      LOGICAL FLAG(N)
      DOUBLE PRECISION CRIT(N), MEMBR(N),DISS(LEN), DISNN(N)
c Var
      INTEGER IM, JJ, JM, I, NCL, J, IND, I2, J2, K, IND1, IND2, IND3
      DOUBLE PRECISION INF, DMIN, D12
c External function
      INTEGER IOFFST
c
c     was 1D+20
      DATA INF/1.D+300/
c
c     unnecessary initialization of im jj jm to keep g77 -Wall happy
c
      IM = 0
      JJ = 0
      JM = 0
C
C  Initializations
C
      DO 10 I=1,N
C        We do not initialize MEMBR in order to be able to restart the
C        algorithm from a cut.
C        MEMBR(I)=1.
 10      FLAG(I)=.TRUE.
      NCL=N
C
C  Carry out an agglomeration - first create list of NNs
C  Note NN and DISNN are the nearest neighbour and its distance
C  TO THE RIGHT of I.
C
      DO 30 I=1,N-1
         DMIN=INF
         DO 20 J=I+1,N
            IND=IOFFST(N,I,J)
            IF (DISS(IND).GE.DMIN) GOTO 20
               DMIN=DISS(IND)
               JM=J
 20         CONTINUE
         NN(I)=JM
         DISNN(I)=DMIN
 30      CONTINUE
C
  400 CONTINUE
C     Next, determine least diss. using list of NNs
      DMIN=INF
      DO 600 I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 600
         IF (DISNN(I).GE.DMIN) GOTO 600
            DMIN=DISNN(I)
            IM=I
            JM=NN(I)
  600    CONTINUE
      NCL=NCL-1
C
C  This allows an agglomeration to be carried out.
C
      I2=MIN0(IM,JM)
      J2=MAX0(IM,JM)
      IA(N-NCL)=I2
      IB(N-NCL)=J2
      CRIT(N-NCL)=DMIN
      FLAG(J2)=.FALSE.
C
C  Update dissimilarities from new cluster.
C
      DMIN=INF
      DO 50 K=1,N
         IF (.NOT.FLAG(K)) GOTO 50
         IF (K.EQ.I2) GOTO 50
         IF (I2.LT.K) THEN
                           IND1=IOFFST(N,I2,K)
                      ELSE
                           IND1=IOFFST(N,K,I2)
         ENDIF
         IF (J2.LT.K) THEN
                           IND2=IOFFST(N,J2,K)
                      ELSE
                           IND2=IOFFST(N,K,J2)
         ENDIF
         IND3=IOFFST(N,I2,J2)
         D12=DISS(IND3)
C
C  WARD.D1 MINIMUM VARIANCE METHOD - IOPT=1 and 2.
C  ### Previous and next lines modified by PL. Same below.
         IF (IOPT.LE.2) THEN
            DISS(IND1)=(MEMBR(I2)+MEMBR(K))*DISS(IND1)+
     X                 (MEMBR(J2)+MEMBR(K))*DISS(IND2) - MEMBR(K)*D12
            DISS(IND1)=DISS(IND1) / (MEMBR(I2)+MEMBR(J2)+MEMBR(K))
         ENDIF
C
C  SINGLE LINK METHOD - IOPT=3.
C
         IF (IOPT.EQ.2) THEN
            DISS(IND1)=MIN(DISS(IND1),DISS(IND2))
         ENDIF
C
C  COMPLETE LINK METHOD - IOPT=4.
C
         IF (IOPT.EQ.3) THEN
            DISS(IND1)=MAX(DISS(IND1),DISS(IND2))
         ENDIF
C
C  AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=5.
C
         IF (IOPT.EQ.4) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2))/
     X                 (MEMBR(I2)+MEMBR(J2))
         ENDIF
C
C  MCQUITTY'S METHOD - IOPT=6.
C
         IF (IOPT.EQ.5) THEN
            DISS(IND1)=0.5D0*DISS(IND1)+0.5D0*DISS(IND2)
         ENDIF
C
C  CENTROID METHOD - IOPT=7.
C
         IF (IOPT.EQ.7) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2)-
     X                  MEMBR(I2)*MEMBR(J2)*D12/(MEMBR(I2)+MEMBR(J2)))/
     X          (MEMBR(I2)+MEMBR(J2))
            ENDIF
C
C  MEDIAN (GOWER'S) METHOD - IOPT=8.
C
         IF (IOPT.EQ.8) THEN
            DISS(IND1)=0.5D0*DISS(IND1)+0.5D0*DISS(IND2)-0.25D0*D12
         ENDIF
C
C  Beta-flexible method (Lance & Williams l966a, 1967c) - IOPT=9.
C
         IF (IOPT.EQ.9) THEN
            DISS(IND1)=alpha*DISS(IND1)+alpha*DISS(IND2)-beta*D12
         ENDIF
C
 50      CONTINUE
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
C
C  Update list of nearest neighbours
C
      DO 900 I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 900
C        (Redetermine NN of I:)
         DMIN=INF
         DO 870 J=I+1,N
            IF (.NOT.FLAG(J)) GOTO 870
            IND=IOFFST(N,I,J)
            IF (DISS(IND).GE.DMIN) GOTO 870
               DMIN=DISS(IND)
               JJ=J
  870       CONTINUE
         NN(I)=JJ
         DISNN(I)=DMIN
  900    CONTINUE
C
C  Repeat previous steps until N-1 agglomerations carried out.
C
      IF (NCL.GT.1) GOTO 400
C
C  For ward.D2, take the square root of the fusion heights
      if(IOPT.EQ.2) 
         DO 910 I=1,N-1
  910    CRIT(I) <- sqrt(CRIT(I))
C
      RETURN
      END
C     of HCLUST()
C
C
      INTEGER FUNCTION IOFFST(N,I,J)
C  Map row I and column J of upper half diagonal symmetric matrix
C  onto vector.
      INTEGER N,I,J
      IOFFST=J+(I-1)*N-(I*(I+1))/2
      RETURN
      END
