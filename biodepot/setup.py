from setuptools import setup
setup(name="Bwb_core",
      packages=["Bwb_core"],
      package_data={"Bwb_core": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Bwb-core = Bwb_core"},
      )
setup(name="RNA-seq",
      packages=["RNA_seq"],
      package_data={"RNA_seq": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "RNA-seq = RNA_seq"},
      )
setup(name="Utilities",
      packages=["Utilities"],
      package_data={"Utilities": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Utilities = Utilities"},
      )
setup(name="Miscellaneous",
      packages=["Miscellaneous"],
      package_data={"Miscellaneous": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Miscellaneous = Miscellaneous"},
      )
setup(name="User",
      packages=["User"],
      package_data={"User": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "User = User"},
      )
