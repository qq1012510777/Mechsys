#!/bin/bash

set -e

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo
echo "****************************************************************************"
echo "* You can call this script with an option to force recompiling everything  *"
echo "* and/or an option to also download packages                               *"
echo "*                                                  recompile     download  *"
echo "*                                                          |     |         *"
echo "* Example:                                                 V     V         *"
echo "*   sh $MECHSYS_ROOT/mechsys/scripts/do_compile_deps.sh {0,1} {0,1}        *"
echo "*                                                                          *"
echo "* By default, the code will not be re-compiled if this was already done.   *"
echo "*                                                                          *"
echo "* By default, -j4 (4 processors) is given to make. This can be changed by  *"
echo "* modifying the variable NPROCS in this script.                            *"
echo "*                                                                          *"
echo "****************************************************************************"

NPROCS=4

if [ -d "$MECHSYS_ROOT/mechsys" ]; then
    echo
    echo "Found: $MECHSYS_ROOT/mechsys ==> OK"
else
    echo
    echo "Directory named 'mechsys' does not exist"
    echo "Please, download 'mechsys' first:"
    echo
    echo "   hg clone http://hg.savannah.nongnu.org/hgweb/mechsys/"
    echo
    exit 1
fi

RECOMPILE=0
if [ "$#" -gt 0 ]; then
    RECOMPILE=$1
    if [ "$RECOMPILE" -lt 0 -o "$RECOMPILE" -gt 1 ]; then
        echo
        echo "The option for re-compilation must be either 1 or 0. ($1 is invalid)"
        echo
        exit 1
    fi
fi

FORCEDOWNLOAD=0
if [ "$#" -gt 1 ]; then
    FORCEDOWNLOAD=$2
    if [ "$FORCEDOWNLOAD" -lt 0 -o "$FORCEDOWNLOAD" -gt 1 ]; then
        echo
        echo "The option for downloading and compilation of additional packages must be either 1 or 0. ($2 is invalid)"
        echo
        exit 1
    fi
fi

test -d $MECHSYS_ROOT/pkg || mkdir $MECHSYS_ROOT/pkg

VER_TRIANGLE=1.6
VER_TETGEN=1.4.3
VER_VORO=0.3.1
VER_OPENMPI=1.5.4
VER_PARMETIS=3.2.0 # 4.0.2 # MUMPS does not work with this one
VER_MUMPS=4.10.0
VER_IGRAPH=0.5.4
VER_SOPLEX=1.5.0
VER_VTK=5.8.0
VER_VTK_MAJOR=5.8

compile_scalapack() {
    INC_OPENMPI=/usr/local/lib/openmpi/include
    LIB_LAPACK=/usr/lib/liblapack.so
    LIB_BLAS=/usr/lib/libblas.so
    LDIR=$MECHSYS_ROOT/pkg/scalapack_installer/lib
    python setup.py --notesting --mpiincdir=$INC_OPENMPI --lapacklib=$LIB_LAPACK --blaslib=$LIB_BLAS
}

compile_mumps() {
    cp $MECHSYS_ROOT/mechsys/patches/mumps/Makefile.inc .
    make clean
    make -j$NPROCS
}

parmetis_links() {
    LDIR=$MECHSYS_ROOT/pkg/parmetis-$VER_PARMETIS
    MACH=`uname -m`
    ln -s $LDIR/build/Linux-$MACH/libmetis/libmetis.a $LDIR/
    ln -s $LDIR/build/Linux-$MACH/libparmetis/libparmetis.a $LDIR/
}

error_message() {
    echo
    echo
    echo "    [1;31m Error: $1 [0m"
    echo
    echo
}

download_and_compile() {
    PKG=""
    PKG_DIR=""
    EXT=tar.gz
    LOCATION=""
    EXTRA_CONF=""
    EXTRA_CMD=""
    CONF_PRMS=""
    IS_SVN=0
    DO_PATCH=0
    DO_CONF=0
    DO_CMAKECONF=0
    DO_MAKE=1
    DO_MAKE_INST=0
    case "$1" in
        triangle)
            PKG=triangle$VER_TRIANGLE
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        tetgen)
            PKG=tetgen$VER_TETGEN
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        voro)
            PKG=voro++$VER_VORO
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            DO_MAKE=0
            ;;
        openmpi)
            PKG=openmpi-$VER_OPENMPI
            EXT=tar.bz2
            LOCATION=http://www.open-mpi.org/software/ompi/v1.5/downloads/$PKG.$EXT
            DO_CONF=1
            DO_MAKE_INST=1
            ;;
        parmetis)
            #PKG=parmetis-$VER_PARMETIS  # for 4.0.2
            #EXTRA_CONF="make config"    # for 4.0.2
            #EXTRA_CMD=parmetis_links    # for 4.0.2
            PKG=ParMetis-$VER_PARMETIS
            LOCATION=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$PKG.$EXT
            LOCATION=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD/$PKG.$EXT
            ;;
        scalapack)
            PKG=scalapack_installer
            PKG_DIR=scalapack_installer
            EXT=tgz
            LOCATION=http://www.netlib.org/scalapack/$PKG.$EXT
            DO_MAKE=0
            EXTRA_CMD=compile_scalapack
            ;;
        mumps)
            PKG=MUMPS_$VER_MUMPS
            #LOCATION=""
            LOCATION=http://mumps.enseeiht.fr/$PKG.$EXT
            DO_MAKE=0
            EXTRA_CMD=compile_mumps
            ;;
        igraph)
            PKG=igraph-$VER_IGRAPH
            LOCATION=http://sourceforge.net/projects/igraph/files/C%20library/$VER_IGRAPH/$PKG.$EXT
            DO_CONF=1
            ;;
        soplex)
            PKG=soplex-$VER_SOPLEX
            EXT=tgz
            LOCATION=http://soplex.zib.de/download/$PKG.$EXT
            ;;
        vtk)
            PKG=vtk-$VER_VTK
            PKG_DIR=VTK
            LOCATION=http://www.vtk.org/files/release/$VER_VTK_MAJOR/$PKG.$EXT
            DO_CMAKECONF=1
            DO_MAKE_INST=1
            ;;
        *)
            error_message "download_and_compile: __Internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"

    # change into the packages directory
    cd $MECHSYS_ROOT/pkg

    # package filename and directory
    PKG_FILENAME=$PKG.$EXT
    if [ -z "$PKG_DIR" ]; then PKG_DIR=$PKG; fi

    # check for package that must be existing (cannot be downloaded)
    if [ -z "$LOCATION" ]; then
        if [ ! -e "$PKG_FILENAME" ]; then
            error_message "Please download <$PKG_FILENAME> first"
            return
        fi
    fi

    # (re)compile or return (erasing existing package) ?
    if [ "$IS_SVN" -eq 0 ]; then
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Erasing existing $PKG_DIR"
                rm -rf $MECHSYS_ROOT/pkg/$PKG_DIR
            else
                echo "    Using existing $PKG_DIR"
                return
            fi
        fi
    else
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Updating existing $PKG SVN repository"
                cd $PKG_DIR
                svn up
                cd $MECHSYS_ROOT/pkg
            else
                echo "    Using existing $PKG SVN repository in $PKG_DIR"
                return
            fi
        fi
    fi

    # download package
    if [ "$IS_SVN" -eq 0 ]; then
        if [ "$FORCEDOWNLOAD" -eq 1   -o   ! -e "$PKG_FILENAME" ]; then
            if [ -e "$PKG_FILENAME" ]; then
                echo "    Removing existing <$PKG_FILENAME>"
                rm $PKG_FILENAME
            fi
            echo "    Downloading <$PKG_FILENAME>"
            wget $LOCATION
        fi
    else
        if [ ! -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            echo "    Downloading new $PKG SVN repository"
            svn co $LOCATION $PKG
        fi
    fi

    # uncompress package
    if [ "$IS_SVN" -eq 0 ]; then
        echo "        . . . uncompressing . . ."
        if [ "$EXT" = "tar.bz2" ]; then
            tar xjf $PKG_FILENAME
        else
            tar xzf $PKG_FILENAME
        fi
    fi

    # change into the package directory
    cd $PKG_DIR

    # patch
    if [ "$DO_PATCH" -eq 1 ]; then
        echo "        . . . patching . . ."
        sh $MECHSYS_ROOT/mechsys/patches/${1}/do_patch.sh
    fi

    # configure
    if [ "$DO_CONF" -eq 1 ]; then
        echo "        . . . configuring . . ."
        ./configure $CONF_PRMS 2> /dev/null
    fi

    # extra configuration
    if [ ! -z "$EXTRA_CONF" ]; then
        echo "        . . . extra configuration . . . . . ."
        $EXTRA_CONF
    fi

    # cmake configuration
    if [ "$DO_CMAKECONF" -eq 1 ]; then
        echo "        . . . configuring using cmake . . ."
        cmake . 2> /dev/null
    fi

    # compilation
    if [ "$DO_MAKE" -eq 1 ]; then
        echo "        . . . compiling . . ."
        make -j$NPROCS > /dev/null 2> /dev/null
    fi

    # make install
    if [ "$DO_MAKE_INST" -eq 1 ]; then
        echo "        . . . installing . . ."
        sudo make install > /dev/null 2> /dev/null
        echo "           . . ldconfig . ."
        sudo ldconfig
    fi

    # execute specific command
    if [ ! -z "$EXTRA_CMD" ]; then
        echo "        . . . command . . . . . ."
        $EXTRA_CMD
    fi

    # finished
    echo "        . . . finished . . . . . "
}

download_and_compile triangle
download_and_compile tetgen
download_and_compile voro
download_and_compile openmpi
download_and_compile parmetis
download_and_compile scalapack
download_and_compile mumps
download_and_compile igraph
download_and_compile soplex
download_and_compile vtk

echo
echo "Finished ###################################################################"
echo
