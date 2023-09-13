#!/bin/bash
set -e

# Set variables
VERSION=$1

# Update Cargo version
sed -i "s/^__version__ = .*$/version = \"$VERSION\"/" magphylogeny/__init__.py
