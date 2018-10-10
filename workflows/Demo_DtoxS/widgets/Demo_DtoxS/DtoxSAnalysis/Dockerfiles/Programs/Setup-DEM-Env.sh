#!/bin/bash
# This script should be at the top level of a directory tree for DEG/DEP analysis:
# DEG/DEP
# ├── Configs
# ├── Programs
# ├── Counts
# ├── Params
# ├── Scripts
# └── Results


help_msg ()
{
	echo "Set up a directory tree for DEG/DEP analysis." 1>&2
	echo "Usage: "$1" [dest dir] [src dir]" 1>&2
	echo "       [dest dir] is the destination top directory for DEG/DEP analysis." 1>&2
	echo "       [src dir] is the source top directory for DEG/DEP analysis [Default: ..] (Optional)." 1>&2
	exit 1
}

check_dir ()
{
	if [ ! -d "$1" ]; then
		echo "ERROR: "$1" is not found!" 1>&2
		exit 1
	fi
}

# The main program begins here

PROG_DIR="$(dirname "$0")"
PROG_NAME="$(basename "$0")"
if [ $# -lt 1 ]; then
	help_msg "${PROG_NAME}"
else
	if [ $# -lt 2 ]; then
		# Destination top DEG/DEP directory.
		DEM_DIR="$1"
		# Source top DEG/DEP directory.
		pushd "${PROG_DIR}" &> /dev/null
		SRC_DIR="$(dirname "$(pwd -L)")"
		popd &> /dev/null
	elif [ $# -lt 3 ]; then
		# Destination top DEG/DEP directory.
		DEM_DIR="$1"
		# Source top DEG/DEP directory.
		SRC_DIR="$2"
	else
		help_msg "${PROG_NAME}"
	fi
fi

# Check source directory tree.
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Programs")
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Configs")
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Counts")
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Params")
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Scripts")
SRC_DEM_DIRS=("${SRC_DEM_DIRS[@]}" "${SRC_DIR}/Results")
for DIR in "${SRC_DEM_DIRS[@]}"; do
	check_dir "${DIR}"
done

# Create the destination top DEG/DEP directory if necessary.
if [ ! -d "${DEM_DIR}" ]; then
	mkdir -p "${DEM_DIR}"
fi

# Copy the source DEG/DEP directories to the destination top DEG/DEP directory.
echo "Copy DEG/DEP directory tree to" "${DEM_DIR}..."
for DIR in "${SRC_DEM_DIRS[@]}"; do
	if [ "$(basename "${DIR}")" == "Programs" ]; then
		COPY_OPTS="-rL"
	else
		COPY_OPTS="-rP"
	fi
	echo cp "${COPY_OPTS}" "${DIR}" "${DEM_DIR}/"
	cp "${COPY_OPTS}" "${DIR}" "${DEM_DIR}/"
done

exit 0
