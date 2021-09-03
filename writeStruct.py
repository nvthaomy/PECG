"""
print out structure of PEs compatible with main.py
input: 
	packmol mix file: mix.inp, 
	charge file in ~/bin/SCOUTff/structlib/PAA_PAH
"""

pdbPrefix = 'AA6f0.5' 
molPrefix = 'PAA'
mixFile = 'mix.inp'
chargeFile = 'PAAf0.5charge.txt'
#column index in chargeFile
chainName_col = 1
seq_col = 2

atomMap = {'h': 'A', 't': 'A', 'u': 'A', 'd': 'A', 'H':'A-', 'T':'A-', 'U':'A-', 'D':'A-'}
#atomMap = {'h': 'B', 't': 'B', 'u': 'B', 'd': 'B', 'H':'B+', 'T':'B+', 'U':'B+', 'D':'B+'}

pdbs = []
nmols = []
structLib = {}
nmolLib = {}
f = open(mixFile,'r')
lines = f.readlines()

i = 0
while i <= len(lines)-1:
    line = lines[i]
    if line.startswith('structure'):
        pdb = line.split('/')[-1]
        pdb = pdb.split('.pdb')[0]
        if pdbPrefix in pdb:
            pdbs.append(pdb)
            i += 1
            line = lines[i]
            nmols.append(int(line.split()[-1]))
    i += 1        

chargeLib = {}
f = open(chargeFile,'r')
lines = f.readlines()
for line in lines:
    if not line.startswith('#'):
        key = line.split()[chainName_col] 
        val = line.split()[seq_col]
        tmp_struc = []
        for l in val:
            tmp_struc.append(atomMap[l])
        chargeLib.update({key: tmp_struc})

# make structure libraries for sim
for i,pdb in enumerate(pdbs):
    structLib.update({molPrefix+'{}'.format(i): chargeLib[pdb]})
    nmolLib.update({molPrefix+'{}'.format(i): nmols[i]})
print('structures \n{}'.format(structLib))
print('number of molecules \n{}'.format(nmolLib))
print('molecule names \n{}'.format(structLib.keys()))
