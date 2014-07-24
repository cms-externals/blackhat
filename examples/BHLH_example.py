import BHLH

orderFileText="""
# BLHA options

MatrixElementSquareType CHsummed
CorrectionType QCD
MassiveParticleScheme OnShell
IRsubtractionMethod None
SubdivideProcess No
IRregularisation CDR

# BH options

Z_mass 90.0
W_mass 80.0
Z_width 10


11 -11 -> 1 -1 
11 -11  -> 1 -1 

"""

of=open('orderFile.lh','w')
of.write(orderFileText)
of.close()


BHLH.SignContract('orderFile.lh','contractFile.lh')
BHLH.Init('contract_file.lh')


momenta=[ 
# E  x  y  z    m
  1, 0, 0, 1,   0,
  1, 0, 0,-1,   0,
  1, 0, 1, 0,   0,
  1, 0,-1, 0,   0
]

double,single,finite,tree=BHLH.EvalSubprocess(1,momenta,10,0,0)

print 'A=%s/eps^2 + %s/eps + %s' % (double,single,finite)
