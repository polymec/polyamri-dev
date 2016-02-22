# Install required software
brew update
brew install gcc openmpi wget

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

