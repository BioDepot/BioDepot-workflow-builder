foldable pip install -U setuptools pip codecov

# Don't install PyQt5 if PyQt4 is requested
[ "$PYQT4" ] && sed -i '/pyqt5/Id' requirements-doc.txt

if [[ $TRAVIS_PYTHON_VERSION == 3.4 ]]; then pip install pandas==0.20.3; fi

# Install dependencies sequentially
cat requirements-core.txt \
    requirements-gui.txt \
    requirements-dev.txt \
    requirements-opt.txt \
    requirements-doc.txt |
    while read dep; do
        dep="${dep%%#*}"  # Strip the comment
        [ "$dep" ] &&
            foldable pip install $dep
    done

# Create a source tarball from the git checkout
foldable python setup.py sdist
# Create a binary wheel from the packed source
foldable pip wheel --no-deps -w dist dist/Orange3-*.tar.gz
# Install into a testing folder
ORANGE_DIR="$(pwd)"/build/travis-test
mkdir -p "$ORANGE_DIR"
pip install --no-deps --target "$ORANGE_DIR"  dist/Orange3-*.whl

cd $TRAVIS_BUILD_DIR
