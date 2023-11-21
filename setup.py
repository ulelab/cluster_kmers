from setuptools import setup

requirements = [
  'textdistance',
  'seaborn',
  'scipy',
  'matplotlib',
  'pandas',
  'scikit-learn',
  'numpy',
  'scikit-bio',
  'logomaker',
  'seqlogo'
]

setup(name='cluster_kmers',
      version='0.0.1',
      description='cluster_kmers',
      url='http://github.com/ulelab/cluster_kmers',
      author='Klara Kuret',
      author_email='klara.kuret@ki.si',
      license='MIT',
      install_requires=requirements,
      packages=['cluster_kmers'],
      include_package_data=True,
      package_data={'cluster_kmers': ['res/*']},
      python_requires=">=3.8, <3.12",
      entry_points={"console_scripts": ["cluster_kmers=cluster_kmers.run:main",],},
      zip_safe=False)
