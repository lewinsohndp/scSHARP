from setuptools import setup, find_packages

setup(
    name='scSHARP',
    version='0.1.1',
    packages=find_packages(where='src', 
                           include=['scSHARP_tool', 'scSHARP_tool.*'],
                           exclude=['test.*'])
                        #    install_requires=[
                        #        'pandas',
                        #        'scikit-learn',
                        #        'torch-cluster',
                        #        'json',
                        #        'pyg',
                        #        'anndata',
                        #        'scanpy',
                        #        'rpy2',
                        #        'captum'
                        #    ]
)