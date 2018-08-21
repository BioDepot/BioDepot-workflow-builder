from setuptools import setup
setup(name="Demo_DtoxS",
      packages=["Demo_DtoxS"],
      package_data={"Demo_DtoxS": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Demo_DtoxS = Demo_DtoxS"},)
      
setup(name="Demo_kallisto",
      packages=["Demo_kallisto"],
      package_data={"Demo_kallisto": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Demo_kallisto = Demo_kallisto"},)
      
setup(name="Demo_STAR",
      packages=["Demo_STAR"],
      package_data={"Demo_STAR": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Demo_STAR = Demo_STAR"},)
