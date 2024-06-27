# -*- coding: utf-8 -*-
""" Utility Functions

"""
from typing import Union, Optional

import sys
import os
import re
from itertools import groupby

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2023 Matthew L. Bendall"

def numstr(x: str) -> Union[int, float, None]:
    """ Convert argument to numeric type.

    Attempts to convert argument to an integer. If this fails, attempts to
    convert to float. If both fail, return as string.

    Args:
        x: Argument to convert.

    Returns:
        The argument as int, float, or string.

    """
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return str(x)


def attrstr(x):
    """ Convert argument to dictionary.

    Converts a GTF attribute string to a dictionary. The attribute string
    should be a semicolon-separated list of tag-value pairs. See
    http://www.ensembl.org/info/website/upload/gff.html.

    Args:
        x: Argument to convert.

    Returns:
        A dictionary containing the tag-value pairs.

    """
    if type(x) is dict:
        return x
    ret = {}
    for t in re.findall('(\S+)\s+"([\s\S]*?)";', x):
        ret[t[0]] = numstr(t[1])
    return ret

def nonedot(v: Union[str, None]) -> str:
    """ Return a dot if v is None"""
    return '.' if v is None else str(v)

def dotnone(s: str) -> Union[str, None]:
    """ Return None if v is a dot """
    return None if s.strip() == '.' else s



def fileopen(fun):
    """ Open read-only filehandle if first argument is string

    This function can be used to decorate functions that expect a filehandle as
    the first argument. If the first argument is a string, it is assumed to be
    the path to the input file, and the file is opened in read mode.

    Args:
        fun: A function that expects a read-only filehandle as the first arg.

    Returns:
        function: A function that expects a path (str) OR filehandle as the
            first argument.

    """
    def wrapper(*args, **kwargs):
        if isinstance(args[0], str):
            try:
                fh = open(args[0], 'r')
                return fun(fh, *args[1:], **kwargs)
            except IOError as e:
                raise e
        elif hasattr(args[0], 'read'):
                return fun(*args, **kwargs)
        else:
            raise IOError("Unknown argument for fileopen '%s'" % str(args[0]))
    return wrapper


@fileopen
def tsv(infile, comment=None):
    """ Returns a generator for tab-delmited file.

    Args:
        infile: Input file as a file-like object
        comment (str): Rows beginning with this string will be ignored.

    Returns:
        generator: A generator that yields each row in the file as a list.

    """
    if comment is None:
        return (l.strip('\n').split('\t') for l in infile)
    else:
        return (l.strip('\n').split('\t') for l in infile if not l.startswith(comment))



def collapse_list(l):
    """ Collapse identical adjacent elements in list

    Identical adjacent elements are collapsed into a single element. Similar to
    the `uniq` linux utility.

    Args:
        l: List of objects

    Returns:
        list: List of objects with adjacent identical elements collapsed.

    Examples:
        >>> collapse_list([1, 1, 1, 2, 2, 3, 3, 1, 1,])
        [1, 2, 3, 1]
    """
    return [k for k, g in groupby(l)]


def overlap_length(a, b):
    """ Calculate the length of overlap for two intervals.

    Args:
        a (int, int): First interval (as 2-tuple of integers).
        b (int, int): Second interval (as 2-tuple of integers).

    Returns:
        int: The length of overlap for a and b, or 0 if there is no overlap.

    Examples:
        >>>
        10
        >>> overlap_length((60, 90), (30, 45))
        0
        >>> overlap_length((20, 40), (36, 86))
        4
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def raw_input_stderr(*args):
    """ Prompt user for input using standard error.

    Same as built-in function `raw_input`, but prompt is written to standard
    error instead of standard output.

    Args:
        *args: Arguments passed to raw_input. If `prompt` is present it is
            written to standard error.

    Returns:
        str: User input (with no trailing newline)

    """
    sys.stdout = sys.stderr
    x = input(*args)
    sys.stdout = sys.__stdout__
    print(x)
    return x    


def wraplines(s, wrap=60):
    """ Wrap lines

    Args:
        s: String
        wrap: Number of characters per line

    Returns:
        str: String with newlines.

    """
    return '\n'.join(s[i:(i+wrap)] for i in range(0, len(s), wrap))
