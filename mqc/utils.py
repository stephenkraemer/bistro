import magic
import os
import gzip


def open_gzip_or_plain_file(filepath, mode='rt'):
    """Recognizes and opens plain and gzip text files

    Throws informative exceptions if
    - file does not exist
    - file is empty

    Uses python-magic to identify the file type; python-magic is a wrapper
    for libmagic, so this should be quite robust.

    Args:
        filepath (str): relative or absolute path to file
        mode (str): Defaults to 'rt'

    Returns:
        file-like object (text file)
    """


    # Note: There is a race condition here, as the file may
    # be deleted between those checks and the opening of the file
    # Improve or don't use where such cases are a problem

    # These checks work even if there are no access permissions to the file


    #TODO-low-prio: What happens if there are no access permissions to the dir
    #      containing the file? Will these checks fail with cryptic errors?
    if not os.path.exists(filepath):
        raise FileNotFoundError(f'File {filepath} not found')
    if os.path.getsize(filepath) == 0:
        raise IOError(f'File {filepath} is empty')

    try:
        filetype_mime = magic.from_file(filepath, mime=True)
    except OSError:
        raise OSError(f"Can't open {filepath}")

    # Need to account for different gzip media types
    # see: https://tools.ietf.org/id/draft-levine-application-gzip-03.html
    gzip_media_types = ['application/gzip', 'application/gzip-compressed',
                        'application/gzipped', 'application/x-gzip',
                        'application/x-gzip-compressed', 'gzip/document']
    if filetype_mime in gzip_media_types:
        try:
            fobj = gzip.open(filepath, mode)
        except OSError:
            raise OSError(f"Can't open {filepath}")
    else:
        try:
            fobj = open(filepath, mode)
        except OSError:
            raise OSError(f"Can't open {filepath}")

    return fobj
