
"""
A higher level HMMER tlbout representation in Python

=======
License
=======

Copyright (c) 2022 Mateo Torres <mateo.torres@fgv.br>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""

__author__ = "Mateo Torres"
__email__ = "mateo.torres@fgv.br"
__copyright__ = "Copyright (c) 2022, Tamas Nepusz"
__license__ = "MIT"
__version__ = "0.1"

__all__ = ["HmmerTbloutLine", "HmmerTbloutFile"]

from dataclasses import dataclass
from gzip import GzipFile

@dataclass
class HmmerTbloutLine(object):
    
    target_name: str
    target_accession: str
    query_name: str
    query_accession: str
    e_value: float
    score: float
    bias: float
    e_value_domain: float
    score_domain: float
    bias_domain: float
    exp: float
    reg: float
    clu: float
    ov: float
    env: float
    dom: float
    rep: float
    inc: float
    description_of_target: str
    
    @staticmethod
    def from_line(line):
        """Constructs an annotation from a line"""
        args = line.strip().split(maxsplit=18)
        for i in range(4, 18):
            args[i] = float(args[i])
        return HmmerTbloutLine(*args)

                
class HmmerTbloutFile(object):
    def __init__(self, fp):
        if isinstance(fp, str):
            if fp[-3:] == ".gz":
                self.fp = GzipFile(fp)
            else:
                self.fp = open(fp)
        else:
            self.fp = fp
        
    def annotations(self):
        for ln, line in enumerate(self.fp):
            if not line or line[0] == '#':
                continue
            try:
                yield HmmerTbloutLine.from_line(line)
            except TypeError as ex:
                raise SyntaxError("cannot parse HMMER line", ln+1)
        
    def __iter__(self):
        return self.annotations()
