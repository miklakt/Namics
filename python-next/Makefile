.PHONY: all test clean
all:
	python3 setup.py build_ext --inplace

test: all
	python3 -m unittest discover -s tests -v

clean:
	rm -rf build namics/_namics_kernels*.so
