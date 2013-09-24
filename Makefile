.PHONY: install all clean src

all: install
install:  
	cp -R src/* install/
src:
clean:
	rm -Rf install
	mkdir install
