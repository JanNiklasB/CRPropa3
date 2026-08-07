#!/bin/bash

cat @CRPropa_BINARY_DIR@/compile_commands.json \
	| grep -e command -m 1 \
	| sed -E "s/ /\n/g" \
	| grep -e "-DCRPROPA" -e "-DFAST_WAVES" -e "-DHAVE_SIMD" \
	| sed -E "s/\"//g" \
	| sed -E 's/\\//g' \