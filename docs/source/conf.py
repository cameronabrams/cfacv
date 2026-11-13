# Sphinx configuration for the cfacv documentation.

import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

project = 'cfacv'
copyright = '2009-2026, Cameron F Abrams'
author = 'Cameron F Abrams, Jasmine Gardner, Gourav Shrivastav, Prince Gyasi'

# Read, never copied: VERSION is the one place the version is written, and
# scripts/release.sh rewrites it.
with open(os.path.join(ROOT, 'VERSION')) as f:
    version = release = f.read().strip()

extensions = [
    'myst_parser',
    'sphinx_copybutton',
]

myst_enable_extensions = ['dollarmath', 'colon_fence']
myst_heading_anchors = 3

exclude_patterns = ['build']

html_theme = 'furo'
html_title = f'cfacv {release}'
