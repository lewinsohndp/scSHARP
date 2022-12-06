from setuptools import find_packages, find_namespace_packages
from distutils.core import setup

if __name__== '__main__':
    setup(include_package_data=True,
          description='test',
          long_description="""test""",
          url='test',
          version='0.1.2',
          packages=find_namespace_packages(),
          zip_safe=False,
          setup_requires=[],
          install_requires=['pandas>=1.3.4',
                            'scikit-learn>=1.0.1',
                            'torch>=1.10.0',
                            'torch-cluster>=1.5.9',
                            #'json',
                            'torch-geometric>=2.0.2',
                            'anndata>=0.8.0',
                            'scanpy>=1.9.1',
                            'rpy2>=3.5.3',
                            'captum>=0.5.0',
                            'torch_sparse==0.6.12',
                            'torch_scatter==2.0.9',
                            'torch_spline-conv==1.2.1'],
          package_data={
                '': ['scSHARP/configs/*']
          },
          name='scSHARP')