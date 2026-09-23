# Third-party notices

The `duckhts.duckdb_extension.wasm` files in this package are the DuckHTS DuckDB extension compiled for duckdb-wasm. DuckHTS itself is distributed here under the GNU General Public License, version 2 or (at your option) any later version; see `LICENSE`. The binaries also contain the components below, each under its own licence, reproduced verbatim from the source named for it. The component list comes from the wasm build's link inputs (`libhts.a`, `libz.a`, `libbz2.a`, `liblzma.a`, and DuckHTS's own sources including the vendored cgranges, libBigWig and VariantKey code).

DuckHTS also contains a computational subset ported from Somalier v0.3.4 (Brent Pedersen and contributors, MIT), and is compiled against the DuckDB C API headers (DuckDB contributors, MIT); no Somalier or DuckDB source files are bundled.

## HTSlib

- Copyright: Genome Research Ltd and contributors
- Licence: MIT/Expat outside `cram/`; Modified BSD 3-Clause for `cram/`
- Text from: third_party/htslib/LICENSE

```text
[Files in this distribution outwith the cram/ subdirectory are distributed
according to the terms of the following MIT/Expat license.]

The MIT/Expat License

Copyright (C) 2012-2026 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.


[Files within the cram/ subdirectory in this distribution are distributed
according to the terms of the following Modified 3-Clause BSD license.]

The Modified-BSD License

Copyright (C) 2012-2026 Genome Research Ltd.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
   nor the names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


[The use of a range of years within a copyright notice in this distribution
should be interpreted as being equivalent to a list of years including the
first and last year specified and all consecutive years between them.

For example, a copyright notice that reads "Copyright (C) 2005, 2007-2009,
2011-2012" should be interpreted as being identical to a notice that reads
"Copyright (C) 2005, 2007, 2008, 2009, 2011, 2012" and a copyright notice
that reads "Copyright (C) 2005-2012" should be interpreted as being identical
to a notice that reads "Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010,
2011, 2012".]
```

## htscodecs (bundled with HTSlib)

- Copyright: Genome Research Ltd and contributors
- Licence: BSD 3-Clause (see text for exceptions)
- Text from: third_party/htslib/htscodecs/LICENSE.md

```text
All files except those explicitly listed below are copyright Genome
Research Limited and are made available under the BSD license.

> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions
> are met:
> 
>     (1) Redistributions of source code must retain the above copyright
>     notice, this list of conditions and the following disclaimer. 
> 
>     (2) Redistributions in binary form must reproduce the above copyright
>     notice, this list of conditions and the following disclaimer in
>     the documentation and/or other materials provided with the distribution.  
>     
>     (3)The name of the author may not be used to endorse or promote
>     products derived from this software without specific prior written
>     permission.
> 
> THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
> IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
> WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
> INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
> (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
> SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
> HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
> STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
> IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
> POSSIBILITY OF SUCH DAMAGE. 

c_range_coder.h is Public Domain, derived from work by Eugene
Shelwien.

rANS_byte.h and rANS_word.h are derived from Fabien Giesen's work and
is Public Domain.  https://github.com/rygorous/ryg_rans This work was
in turn based on the ANS family of entropy encoders as described by
Jarek Duda's paper: http://arxiv.org/abs/1311.2540

> To the extent possible under law, Fabian Giesen has waived all
> copyright and related or neighboring rights to ryg_rans, as
> per the terms of the CC0 license:
> 
>   https://creativecommons.org/publicdomain/zero/1.0
> 
> This work is published from the United States.
```

## libBigWig

- Copyright: Devon Ryan
- Licence: MIT
- Text from: third_party/libBigWig/LICENSE

```text
The MIT License (MIT)

Copyright (c) 2015 Devon Ryan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## cgranges

- Copyright: Dana-Farber Cancer Institute (Heng Li)
- Licence: MIT
- Text from: third_party/cgranges/cgranges.h (header)

```text
The MIT License

   Copyright (c) 2019 Dana-Farber Cancer Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
```

## VariantKey / RegionKey

- Copyright: Nicola Asuni and contributors
- Licence: MIT
- Text from: third_party/variantkey/LICENSE

```text
MIT LICENSE

Copyright (c) 2017-2019, 2023 GENOMICS plc
Copyright (c) 2017-2018 Nicola Asuni - Tecnick.com (binsearch)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

## zlib 1.3.1

- Copyright: Jean-loup Gailly and Mark Adler
- Licence: Zlib
- Text from: vcpkg zlib (wasm32-emscripten) copyright file

```text
Copyright notice:

 (C) 1995-2022 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu
```

## bzip2 1.0.8

- Copyright: Julian Seward
- Licence: bzip2-1.0.6
- Text from: vcpkg bzip2 (wasm32-emscripten) copyright file

```text
--------------------------------------------------------------------------

This program, "bzip2", the associated library "libbzip2", and all
documentation, are copyright (C) 1996-2019 Julian R Seward.  All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. The origin of this software must not be misrepresented; you must 
   not claim that you wrote the original software.  If you use this 
   software in a product, an acknowledgment in the product 
   documentation would be appreciated but is not required.

3. Altered source versions must be plainly marked as such, and must
   not be misrepresented as being the original software.

4. The name of the author may not be used to endorse or promote 
   products derived from this software without specific prior written 
   permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Julian Seward, jseward@acm.org
bzip2/libbzip2 version 1.0.8 of 13 July 2019

--------------------------------------------------------------------------
```

## liblzma (XZ Utils) 5.8.1

- Copyright: Lasse Collin and contributors
- Licence: 0BSD for liblzma (see text)
- Text from: vcpkg liblzma (wasm32-emscripten) copyright file

```text
XZ Utils Licensing
==================

    Different licenses apply to different files in this package. Here
    is a summary of which licenses apply to which parts of this package:

      - liblzma is under the BSD Zero Clause License (0BSD).

      - The command line tools xz, xzdec, lzmadec, and lzmainfo are
        under 0BSD except that, on systems that don't have a usable
        getopt_long, GNU getopt_long is compiled and linked in from the
        'lib' directory. The getopt_long code is under GNU LGPLv2.1+.

      - The scripts to grep, diff, and view compressed files have been
        adapted from GNU gzip. These scripts (xzgrep, xzdiff, xzless,
        and xzmore) are under GNU GPLv2+. The man pages of the scripts
        are under 0BSD; they aren't based on the man pages of GNU gzip.

      - Most of the XZ Utils specific documentation that is in
        plain text files (like README, INSTALL, PACKAGERS, NEWS,
        and ChangeLog) are under 0BSD unless stated otherwise in
        the file itself. The files xz-file-format.txt and
        lzma-file-format.xt are in the public domain but may
        be distributed under the terms of 0BSD too.

      - Translated messages and man pages are under 0BSD except that
        some old translations are in the public domain.

      - Test files and test code in the 'tests' directory, and
        debugging utilities in the 'debug' directory are under
        the BSD Zero Clause License (0BSD).

      - The GNU Autotools based build system contains files that are
        under GNU GPLv2+, GNU GPLv3+, and a few permissive licenses.
        These files don't affect the licensing of the binaries being
        built.

      - The 'extra' directory contains files that are under various
        free software licenses. These aren't built or installed as
        part of XZ Utils.

    The following command may be helpful in finding per-file license
    information. It works on xz.git and on a clean file tree extracted
    from a release tarball.

        sh build-aux/license-check.sh -v

    For the files under the BSD Zero Clause License (0BSD), if
    a copyright notice is needed, the following is sufficient:

        Copyright (C) The XZ Utils authors and contributors

    If you copy significant amounts of 0BSD-licensed code from XZ Utils
    into your project, acknowledging this somewhere in your software is
    polite (especially if it is proprietary, non-free software), but
    it is not legally required by the license terms. Here is an example
    of a good notice to put into "about box" or into documentation:

        This software includes code from XZ Utils <https://tukaani.org/xz/>.

    The following license texts are included in the following files:
      - COPYING.0BSD: BSD Zero Clause License
      - COPYING.LGPLv2.1: GNU Lesser General Public License version 2.1
      - COPYING.GPLv2: GNU General Public License version 2
      - COPYING.GPLv3: GNU General Public License version 3

    If you have questions, don't hesitate to ask for more information.
    The contact information is in the README file.
```
