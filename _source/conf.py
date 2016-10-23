import os, sys
from recommonmark.parser import CommonMarkParser

extensions = []
source_parsers = {'.md': CommonMarkParser}
source_suffix = ['.md']
project = 'Kronecker'
master_doc = 'index'
exclude_patterns = ['_build']
html_title = ''
html_theme = 'alabaster'
html_theme_options = {
    'nosidebar': True,
    'github_user': 'catch22',
    'github_repo': 'kronecker',
    'github_banner': True,
    'show_powered_by': False,
    'page_width': '940px'
}
html_domain_indices = False
html_use_index = False
html_copy_source = False
html_show_copyright = False
html_show_sphinx = False
