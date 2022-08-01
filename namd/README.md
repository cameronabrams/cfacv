The file ``Linux-x86_64.tcl`` in ``arch`` under the NAMD source tree sets environment variables that instruct the compilation where to find TcL headers and libraries.  The version here is modified to use a system-installation of TcL in which ``tcl.h`` is found in ``/usr/include`` and ``libtcl8.6.so`` is found in ``/usr/lib64``.  

If this file is used in building NAMD from source, you will NOT need to download any TcL tarballs from UIUC.  It is recommended that both ``namd2`` and ``cfacv.so`` be built using the same TcL instance, and this config file will ensure this.
