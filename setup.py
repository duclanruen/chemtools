from setuptools import setup

setup(
    name='chemtools',
    version='0.0.1',    
    description='Some useful tools for IUT Students in chemistry',
    url='https://github.com/duclanruen/chemtools',
    author='Olivier Carpentier',
    author_email='duclanruen@gmail.com',
    license='MIT',
    packages=['chemtools'],
    install_requires=['numpy',
                      'scipy',
                      'matplotlib',  
                      'pandas',                   
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3',
    ],
)
