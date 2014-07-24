      program fortran_example
      implicit none
      integer i,nexternal,status,label
      parameter (nexternal=4)
      double precision pmass(nexternal),p(0:4,nexternal)
      character*16 filename     
      double precision couplings(2),mu,virt_wgts(4)

      data (p(i,1),i=0,4)/ 0.5000000E+02 , 0.0000000E+00 ,
     &     0.0000000E+00 , 0.5000000E+02 , 0.0000000E+00 /
      data (p(i,2),i=0,4)/ 0.5000000E+02 , 0.0000000E+00 ,
     &     0.0000000E+00 ,-0.5000000E+02 , 0.0000000E+00 /
      data (p(i,3),i=0,4)/ 0.5000000E+02 , 0.1109243E+02 ,
     &     0.4448308E+02 ,-0.1995529E+02 , 0.0000000E+00 /
      data (p(i,4),i=0,4)/ 0.5000000E+02 ,-0.1109243E+02 ,
     &     -0.4448308E+02 , 0.1995529E+02 , 0.0000000E+00 /

      filename = "contract_file.lh"
      mu =100d0
      label=1

      call OLP_Start(filename//CHAR(0),status)
      call OLP_EvalSubprocess(label,p,mu,couplings,virt_wgts)
         
      write(*,*) virt_wgts

      return
      end
