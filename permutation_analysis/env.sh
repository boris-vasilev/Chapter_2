# Source before any `uv` command: /home is at quota, so keep uv's cache and
# managed Pythons on RDS.
export UV_CACHE_DIR=/rds/project/rds-csoP2nj6Y6Y/biv22/.cache/uv
export UV_PYTHON_INSTALL_DIR=/rds/project/rds-csoP2nj6Y6Y/biv22/.cache/uv-python
