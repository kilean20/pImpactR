import os
from distutils.core import setup
import sys

# ===== compile naff, readTBT ======
CC = 'gfortran'
err = os.system('f2py -c readTBT.f90 -m readTBT')
if err!= 0:
    print('f2py readTBT.f90 failed')
err = os.system(CC+' -c -fPIC -O3 naff.f90')
if err!= 0:
    print(CC+' naff.f90 failed')
os.system('f2py -c pynaff.f90 naff.o -m pyNaff')
if err!= 0:
    print('f2py pynaff.f90 failed')
os.system('mv *.so ./pImpactR/')
os.system('rm *.o *.mod')
#====================================


setup(
    name = "pImpactR",
    version = "0.0.1",
    author = "Kilean Hwang",
    author_email = "kilean@lbl.gov",
    description = ("python wrapper of IMPACTR."),
    license = "Lawrence Berkeley National Laboratory",
    keywords = "IMPACT",
    url = "",
    packages=['pImpactR'],
    #package_data={'pImpactR': ['xmain']},
    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Utilities",
        "License :: Free for non-commercial use",
    ],
    zip_safe=False
)
