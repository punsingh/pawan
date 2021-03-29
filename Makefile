all: lib

lib:
	(mkdir build) || continue
	(cd build && cmake .. && make --no-print-directory && cd .. && cp build/pawan . ) || exit 1

clean: clear
	(rm -rf build) || continue
	(rm pawan) || continue

clear:
	(rm -rf data/*) || continue
	(rm -rf scratch/*) || exit 1
