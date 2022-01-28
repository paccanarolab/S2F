
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

class HmmerTbloutLine(object):
    
    __slots__ = ["target_name", "target_accession", "query_name",
                 "query_accession", "e_value", "score", "bias",
                 "e_value_domain", "score_domain", "bias_domain",
                 "exp", "reg", "clu", "ov", "env", "dom", "rep",
                 "inc", "description_of_target"]
    
    def __init__(self, *args, **kwds):
        """Constructs an annotation. Use keyword arguments to specify the
        values of the different attributes. If you use positional arguments,
        the order of the arguments must be the same as they are in the HMMER 
        tblout file. No syntax checking is done on the values entered"""
        if len(args) == 1 and not kwds:
            args = args[0].strip().split(maxsplit=len(self.__slots__))
        for (name, value) in zip(self.__slots__, args):
            setattr(self, name, value)
        for name, value in kwds.items():
            setattr(self, name, kwds[value])
        for name in self.__slots__:
            if not hasattr(self, name):
                setattr(self, name, "")
        self._polish_attributes()
        
    def _polish_attributes(self):
        """Ensures that the atributes are of the right type"""
        self._ensure_float("e_value")
        self._ensure_float("score")
        self._ensure_float("bias")
        self._ensure_float("e_value_domain")
        self._ensure_float("score_domain")
        self._ensure_float("bias_domain")

    def _ensure_float(self, attr):
        """Ensures that a given attribute is a float and not a string"""
        value = getattr(self, attr)
        if not isinstance(value, float):
            if value == "":
                setattr(self, attr, 0)
            else:
                setattr(self, attr, float(value))
                
class HmmerTbloutFile(object):
    def __init__(self, fp):
        if isinstance(fp, str):
            if fp[-3:] == ".gz":
                from gzip import GzipFile
                self.fp = GzipFile(fp)
            else:
                self.fp = open(fp)
        else:
            self.fp = fp
        self.lineno = 0

        
    def annotations(self):
        for line in self.fp:
            self.lineno += 1
            if not line or line[0] == '#':
                continue
            try:
                line = line.strip()
                yield HmmerTbloutLine(line)
            except TypeError as ex:
                raise SyntaxError("cannot parse HMMER line", self.lineno)
        
    def __iter__(self):
        return self.annotations()