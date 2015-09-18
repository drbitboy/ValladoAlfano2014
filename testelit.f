CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC TESTELIT.F - exercises ELIT.F
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C gfortran testelit.f elit.o -o testelit
      programtestelit
      implicitnone
      doubleprecision phi/1.5d0/,mu/1d0/,fe,ee
      doubleprecision dpr
      dpr = 45d0 / atan(1d0)
      dowhile(.true.)
        callelit(mu,dpr*phi,fe,ee)
        print *,dpr*phi,phi,mu,fe,ee
C       write(*,99999,advance='no')'Phi,Mu:'
        write(*,99999)'Phi,Mu (Ctrl-D to exit):  '
        read(*,*,end=99998) phi, mu
      enddo
99999 format(a,$)
99998 print*
      end
