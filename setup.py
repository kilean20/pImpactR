import os
from distutils.core import setup

# compile readTBT.f90 using f2py
os.system('f2py -c ./pImpactR/readTBT.f90 -m ./oImpactR/readTBT')



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
