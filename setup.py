from setuptools import setup

setup(
    name='xrayflux',
    version='0.1.0',
    packages=['xrayflux'],
    license='MIT',
    url='https://github.com/nicolas-cerardi/xrayflux',
    description="A pyatomdb wrapper that performs spectrum redshifting.",
    long_description="This package uses the APEC model from `pyatomdb` to compute X-ray fluxes of distant sources. If given instrumental response files, it outputs expected count-rates in the requested energy bands.",
    long_description_content_type='markdown',
    author='Nicolas Cerardi',
    author_email='nicolas.cerardi@gmail.com',
    install_requires=[
        'numpy',
        'scipy',
        'pyatomdb',
        'astropy',
        'tqdm'
    ]
)
