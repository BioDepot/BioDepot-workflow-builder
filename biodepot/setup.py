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
setup(name="DToxSDemo",
      packages=["DToxSDemo"],
      package_data={"DToxSDemo": ["icons/*.svg"]},
      entry_points={"orange.widgets": "DToxSDemo = DToxSDemo"},)
setup(name="kallistoDemo",
      packages=["kallistoDemo"],
      package_data={"kallistoDemo": ["icons/*.svg"]},
      entry_points={"orange.widgets": "kallistoDemo = kallistoDemo"},)
