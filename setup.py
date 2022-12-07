from setuptools import setup, find_packages

install_requires = [
    'pandas', 
    'streamlit==1.11.1',
    'molmass',
    'ms-mint'
]

config = {
    'description': 'qc_pipeline',
    'author': 'luis-ponce and soren-wacker',
    'url': 'https://github.com/luis-ponce',
    'download_url': 'https://github.com/luis-ponce/ms_mint_conc',
    'author_email': 'luisfponcinho@gmail.com',
    'version': '0.0.1',
    'install_requires': install_requires,
    'packages': find_packages(),
    'scripts': [],
    'name': 'QC_metabolomics'
}

setup(**config)
