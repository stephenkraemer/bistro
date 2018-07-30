import re
from textwrap import dedent

def to_tsv(multi_line_string: str) -> str:
    """Dedent multi-line string and convert whitespace stretches to tab"""
    return dedent(re.sub(' +', '\t', multi_line_string))
