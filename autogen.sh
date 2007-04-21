#!/bin/sh
test -d autoconf || mkdir autoconf
libtoolize -c --force
aclocal
automake -a -c --foreign
autoconf
