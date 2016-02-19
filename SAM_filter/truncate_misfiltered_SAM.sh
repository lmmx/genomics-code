#!/bin/bash

LAST="$(tail -1 $1)"

let TRUNCATE_SIZE="${#LAST} + 1"
truncate -s -"$TRUNCATE_SIZE" "$1"
