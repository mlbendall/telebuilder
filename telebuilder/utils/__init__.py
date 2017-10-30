import os

def _get_build_from_file():
    if os.path.exists('build.txt'):
        return open('build.txt', 'rU').read().strip()
