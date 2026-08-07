#!/bin/bash

working_dir=@CMAKE_CURRENT_BINARY_DIR@
source_dir=@CMAKE_CURRENT_SOURCE_DIR@
crpropa_binary_dir=@CRPropa_BINARY_DIR@
crpropa_source_dir=@CRPropa_SOURCE_DIR@

rm -rf $working_dir/bindings && mkdir $working_dir/bindings
grep -rh "#include" $crpropa_source_dir/include/* \
	| sed -E "s/#include\s+\"/#include </g" \
	| sed -E s/\"/\>/g \
	| sed -E "s/\s+$//g" \
	| sort -u \
	| grep crpropa \
	> $working_dir/includes_for_bindings.h
binder \
	--root-module crpropa \
	--prefix $working_dir/bindings/ \
	--bind crpropa \
	--config $working_dir/binder.cfg \
	--include-pybind11-stl \
	--annotate-includes \
	--annotate-functions \
	$working_dir/includes_for_bindings.h \
	-p $crpropa_binary_dir/compile_commands.json