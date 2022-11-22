from setuptools import setup

setup(
    name="csp_pyxtal_dftb",
    version="0.1.0",
    install_requires=[
        "pyyaml", "pyxtal",
    ],
    entry_points={
        'console_scripts': [
            'csp_pyxtal=csp_pyxtal_dftb.csp_pyxtal:main',
            'compare_crys=csp_pyxtal_dftb.critic2:main',
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="Crytal structure prediction using PyXtal and DFTB+",
    url="https://github.com/mkatouda/csp_pyxtal_dftb",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
)
