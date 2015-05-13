from setuptools import setup, find_packages

NAME = "AutoTransGeno"

setup(
	name = NAME,
	version="",
	description = "Automate your transgenic and genomic wet work design.",
	author = "Horea Christian",
	author_email = "h.chr@mail.ru",
	url = "https://github.com/TheChymera/"+NAME,
	keywords = ["biology", "genomics", "science"],
	packages = find_packages("src"),
	package_dir = {"":"src"},
	classifiers = [],
	install_requires = [],
	provides = [NAME]
	entry_points = {'console_scripts' : \
			['organamer_reposit = organamer.cli:reposit']
		}
    )
