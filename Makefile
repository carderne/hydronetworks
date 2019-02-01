#!/bin/bash
SHELL := /bin/bash

test:
	# convert the example notebook to a script
	jupyter nbconvert --to script example.ipynb

	# run script
	python3 example.py

	# clean up
	rm example.py