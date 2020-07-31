#!/bin/bash

jupyter nbconvert --to script example.ipynb
sed -i -e 's/plt.*how(.*)//g' example.py
python3 example.py
rm example.py
