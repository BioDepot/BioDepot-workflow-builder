from setuptools import setup
setup(name="RNA-seq",
      packages=["RNA_seq"],
      package_data={"RNA_seq": ["icons/*.svg"]},
      entry_points={"orange.widgets": "RNA-seq = RNA_seq"},
      )
setup(name="Utilities",
      packages=["Utilities"],
      package_data={"Utilities": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Utilities = Utilities"},
      )
setup(name="Miscellaneous",
      packages=["Miscellaneous"],
      package_data={"Miscellaneous": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Miscellaneous = Miscellaneous"},
      )
setup(name="User",
      packages=["User"],
      package_data={"User": ["icons/*.svg"]},
      entry_points={"orange.widgets": "User = User"},
      )
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
setup(name="Jupyter",
      packages=["Jupyter"],
      package_data={"Jupyter": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Jupyter = Jupyter"},)
setup(name="Scripting",
      packages=["Scripting"],
      package_data={"Scripting": ["icons/*.svg"]},
      entry_points={"orange.widgets": "Scripting = Scripting"},)
