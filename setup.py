import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gene_graph_lib",
    version="0.1.9",
    author="DNKonanov",
    author_email="konanovdmitriy@gmail.com",
    description="mini-library for GeneGraph tool and GCB project",
    long_description='mini-library for GeneGraph tool and GCB project',
    long_description_content_type="text/markdown",
    url="https://github.com/DNKonanov/gene_graph_lib",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
