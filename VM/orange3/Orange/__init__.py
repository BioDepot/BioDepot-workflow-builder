from .misc.lazy_module import _LazyModule
from .misc.datasets import _DatasetInfo
from .version import short_version as __version__, git_revision as __git_version__

ADDONS_ENTRY_POINT = "orange.addons"

from Orange import data

for mod_name in [
    "classification",
    "clustering",
    "distance",
    "ensembles",
    "evaluation",
    "misc",
    "modelling",
    "preprocess",
    "projection",
    "regression",
    "statistics",
    "version",
    "widgets",
]:
    globals()[mod_name] = _LazyModule(mod_name)

datasets = _DatasetInfo()

del mod_name

# If Qt is available (GUI) and Qt5, install backport for PyQt4 imports
try:
    import AnyQt.importhooks
except ImportError:
    pass
else:
    if AnyQt.USED_API == "pyqt5":
        # Make the chosen PyQt version pinned
        from AnyQt.QtCore import QObject

        del QObject

        import pyqtgraph  # import pyqtgraph first so that it can detect Qt5

        del pyqtgraph

        AnyQt.importhooks.install_backport_hook("pyqt4")
    del AnyQt
