# conf.py

# 1. Base URL and Sitemap Scheme (Already correct)
html_baseurl = 'https://th0maslavigne.github.io/FEniCSx_GMSH_tutorials/'
sitemap_url_scheme = "{link}"

# # 2. **CRITICAL FIX**: Re-declare extensions list, including sitemap
# # Note: You must include the standard extensions used by Jupyter Book/Sphinx Book Theme
# extensions = [
#     'sphinx.ext.intersphinx',
#     'sphinx_last_updated_by_git',
#     'sphinx_sitemap'  # <-- The required extension
# ]