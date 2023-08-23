DOCS_DIR=docs/
SOURCE_FOLDER=step2-fragment-level-analysis/

doc:
	export SPHINX_APIDOC_OPTIONS=members,private-members,show-inheritance,ignore-module-all && \
	sphinx-apidoc -Mef -o $(DOCS_DIR)/source/modalib/ $(SOURCE_FOLDER)/modalib $(SOURCE_FOLDER)/setup.py
	make -C $(DOCS_DIR) html

clean:
	make -C doc/ clean
	make -C src/ clean

.PHONY: doc clean