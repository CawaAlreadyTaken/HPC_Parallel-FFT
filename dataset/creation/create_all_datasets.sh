#!/bin/bash

for file in *.py; do
  python "$file" || break
done
