#!/bin/bash

cache_location=@CMAKE_CACHEFILE_DIR@/CMakeCache.txt
binary_dir=@CMAKE_CURRENT_BINARY_DIR@

cat $cache_location | grep -i crpropa_ | sed -E "s/:\w+=/ \"/g" | sed -E "s/-ADVANCED//g" | sed -E "/^\/\//d" | sed -E "s/^/set( /g" | sed -E "s/$/\" CACHE STRING \"\" )/g" > $binary_dir/binder_pre_chache_script.cmake