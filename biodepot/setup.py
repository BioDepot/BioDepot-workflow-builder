from setuptools import setup, find_packages

setup(
    name="biodepot",  # Use a common distribution name
    version="0.0.1",
    packages=["RNA_seq", "Utilities", "Miscellaneous", "User", "Jupyter", "Scripting"],
    package_data={
        "RNA_seq": ["icons/*.svg"],
        "Utilities": ["icons/*.svg"],
        "Miscellaneous": ["icons/*.svg"],
        "User": ["icons/*.svg"],
	"Jupyter": ["icons/*.svg"],
	"Scripting": ["icons/*.svg"]
    },
    entry_points={
        "orange.widgets": [
            "RNA-seq = RNA_seq",
            "Utilities = Utilities",
            "Miscellaneous = Miscellaneous",
            "User = User",
            "Jupyter = Jupyter",
	    "Scripting = Scripting"
        ]
    },
)
