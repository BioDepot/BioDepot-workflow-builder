from setuptools import setup

setup(name="Bwb",
      packages=["orangebiodepot"],
      package_data={"orangebiodepot": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Bwb = orangebiodepot"},
      )
setup(name="Bwb-core",
      packages=["bwb"],
      package_data={"bwb": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Bwb-core = bwb"},
      )
setup(name="RNA-seq",
      packages=["RNAseq"],
      package_data={"RNAseq": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "RNA-seq = RNAseq"},
      )
setup(name="Utilities",
      packages=["utilities"],
      package_data={"utilities": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Utilities = utilities"},
      )
setup(name="Miscellaneous",
      packages=["other"],
      package_data={"other": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "Miscellaneous = other"},
      )
setup(name="User",
      packages=["user"],
      package_data={"user": ["icons/*.svg"]},
      classifiers=["Example :: Invalid"],
      entry_points={"orange.widgets": "User = user"},
      )
