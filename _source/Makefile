html:
	git show master:README.md > index.md
	sphinx-build -a -b html . _build
	touch ../.nojekyll
	cp _build/index.html ../index.html
	cp _build/_static/* ../_static/
