#!/bin/bash

abs_pathtopng="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"

pathtohtml=$2

sed -e "s|@PATHTOLOGO|${abs_pathtopng}|" ${pathtohtml}