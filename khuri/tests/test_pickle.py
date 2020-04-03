import contextlib
import pickle
import os
import shutil

import pytest

from khuri.tests.test_omnes import OMNES_LIST


@pytest.fixture
def omnes_function():
    return OMNES_LIST[0]


@contextlib.contextmanager
def environment():
    directory = './TEST_PICKLE_DIRECTORY'
    os.mkdir(directory)
    try:
        yield directory
    finally:
        shutil.rmtree(directory)


def test_omnes(omnes_function):
    with environment() as env:
        path = os.path.join(env, 'omnes.pickle')
        with open(path, 'wb') as pickle_file:
            pickle.dump(omnes_function, pickle_file)
        with open(path, 'rb') as pickle_file:
            func = pickle.load(pickle_file)
            assert func(0.0) == pytest.approx(0.0)
