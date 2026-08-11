#!/bin/sh
# Generate configure script from configure.ac
echo "Generating configure script using autoreconf..."
autoreconf -fi
echo "Done. You may now run ./configure and make."
