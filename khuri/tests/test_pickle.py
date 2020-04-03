import contextlib
import os
import shutil


@contextlib.contextmanager
def environment():
    directory = './TEST_PICKLE_DIRECTORY'
    os.mkdir(directory)
    try:
        yield directory
    finally:
        shutil.rmtree(directory)


def test():
    with environment() as env:
        print(env)
