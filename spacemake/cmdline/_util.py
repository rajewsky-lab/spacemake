import errno
import os

from ..errors import FileWrongExtensionError

LINE_SEPARATOR = '-'*50+'\n'

bool_in_str = ['True', 'true', 'False', 'false']

def assert_file(file_path, default_value='none', extension='all'):
    if file_path == default_value:
        # file doesn't exist but has the default value,
        # so we do not need to assert anything
        return False

    # check if file exists, raise error if not
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), file_path)


    if not file_path.endswith(extension) and extension != 'all':
        raise FileWrongExtensionError(file_path, extension)

    # return true if file exists and every test was good
    return True

def str2bool(var):
    if isinstance(var, bool):
        return var

    if var in ['True', 'true']:
        return True
    elif var in ['False', 'false']:
        return False
    else:
        raise ValueError(f'variable should be boolean, or one of: {bool_in_str}')
