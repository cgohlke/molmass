#! /usr/bin/env python
# -*- coding: utf-8 -*-
# molmass_web.py

# Copyright (c) 2005-2018, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Molecular Mass Calculator - A simple web interface for molmass.py.

Run ``python molmass_web.py`` to execute the script in a local web server.

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:Version: 2018.8.15

Requirements
------------
* `CPython 2.7 or 3.5+ <https://www.python.org>`_
* `Molmass.py 2018.8.15 <https://www.lfd.uci.edu/~gohlke/>`_

Revisions
---------
2018.8.15
    Move module into molmass package.
2018.5.29
    Use CSS flex layout.
    Separate styles from content.
    Run http server from this file's directory.
2018.5.25
    Backwards incompatible changes:
    Rename to molmass_web.py.
    Accept Flask request.args.
    Style and template changes.
    Fix regular expressions.
2018.2.6
    Style and docstring fixes.
    Escape error messages.
2005.x.x
    Initial release.

"""

from __future__ import division, print_function

import os
import re
import sys
from html import escape

if __package__:
    from . import molmass
else:
    import molmass

__version__ = '2018.8.15'
__docformat__ = 'restructuredtext en'


PAGE = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<style type="text/css">
html {{font-size:16px}}
body {{min-width:576px}}
div.header {{display:flex;flex-wrap:wrap;align-items:baseline}}
div.header a {{font-weight:bold;text-decoration:none}}
h1 {{margin:0;padding:0 0.5em 0 0;}}
h1 a {{text-decoration:none;color:#000000}}
form {{display:flex;flex-wrap:nowrap;background-color:#eeeeee;
border:1px solid #aaaaaa;padding:1em;}}
input {{margin-left:0.5em;min-width:4.5em;}}
input[type=text] {{width:100%}}
div.formula {{flex:1 1 100%;min-width:2em;margin-right:1em}}
div.buttons {{white-space:nowrap}}
</style>
{heads}
<meta name="generator" content="molmass_web.py" />
<meta name="robots" content="noarchive" />
<meta name="format-detection" content="telephone=no" />
<meta name="viewport" content="width=608px" />
<title>Molecular Mass Calculator - Christoph Gohlke</title>
</head>
<body>
<div class="header">
<h1><a href="{url}"
title="Molecular Mass Calculator version {version} by Christoph Gohlke"
>Molecular Mass Calculator</a></h1>
<p>by <a href="https://www.lfd.uci.edu/~gohlke/">Christoph Gohlke</a></p>
</div>
<form id="molmass" method="get" action="">
<div><label for="q"><strong>Formula:</strong></label></div>
<div class="formula">
<input type="text" name="q" id="q" value="{formula}" />
</div>
<div class="buttons">
<input type="submit" id="a" value="Submit" />
<input type="reset" value="Reset" onclick="window.location='{url}'"/>
</div>
</form>
<div class="content">
{content}
</div>
</body>
</html>"""

HELP = """<p>This form calculates the molecular mass (average, monoisotopic,
and nominal), the elemental composition, and the mass distribution spectrum
of a molecule given by its chemical formula, relative element weights,
or sequence.</p>
<p>Calculations are based on the
<a href="?q=isotopes" rel="nofollow">isotopic composition of the elements</a>.
Mass deficiency due to chemical bonding is not considered.
</p>
<h3>Examples</h3>
<ul>
<li>Simple chemical formulas:
<a href="?q=H2O" rel="nofollow">H2O</a> or
<a href="?q=CH3COOH" rel="nofollow">CH3COOH</a>
</li>
<li>Specific isotopes:
<a href="?q=D2O" rel="nofollow">D2O</a> or
<a href="?q=[30Si]3O2" rel="nofollow">[30Si]3O2</a>
</li>
<li><a href="?q=groups" rel="nofollow">Abbreviations</a> of chemical groups:
<a href="?q=EtOH" rel="nofollow">EtOH</a> or
<a href="?q=PhOH" rel="nofollow">PhOH</a>
</li>
<li>Simple arithmetic:
<a href="?q=(COOH)2" rel="nofollow">(COOH)2</a> or
<a href="?q=CuSO4.5H2O" rel="nofollow">CuSO4.5H2O</a>
</li>
<li>Relative element weights:
<a href="?q=O: 0.26, 30Si: 0.74" rel="nofollow">O: 0.26, 30Si: 0.74</a>
</li>
<li>Nucleotide sequences:
<a href="?q=CGCGAATTCGCG" rel="nofollow">CGCGAATTCGCG</a> or
<a href="?q=dsrna(CCUU)" rel="nofollow">dsrna(CCUU)</a>
</li>
<li>Peptide sequences:
<a href="?q=MDRGEQGLLK" rel="nofollow">MDRGEQGLLK</a> or
<a href="?q=peptide(CPK)" rel="nofollow">peptide(CPK)</a>
</li>
</ul>
<p>Formulas are case sensitive and &#8217;+&#8217; denotes
the arithmetic operator, not an ion charge.</p>
<div class="disclaimer">
<h3>Disclaimer</h3>
<p>Because this service is provided free of charge, there is no
warranty for the service, to the extent permitted by applicable law.
The service is provided &quot;as is&quot; without warranty of any kind,
either expressed or implied, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose.
The entire risk as to the quality and performance is with you.</p>
</div>
<div class="about">
<h3>About</h3>
<p>Molecular Mass Calculator version {version} by
<a href="https://www.lfd.uci.edu/~gohlke/">Christoph Gohlke</a>.
Source code is available at
<a href="https://www.lfd.uci.edu/~gohlke/code/molmass.py.html"
>molmass.py</a>,
<a href="https://www.lfd.uci.edu/~gohlke/code/molmass_web.py.html"
>molmass_web.py</a>, and
<a href="https://www.lfd.uci.edu/~gohlke/code/elements.py.html"
>elements.py</a>.</p>
</div>"""


def response(form, url, template=PAGE, help=HELP, heads=''):
    """Return HTML document from submitted form data."""
    formula = form.get('q', '')
    if not formula:
        formula = ''
    if formula == 'groups':
        formula = ''
        content = groups()
    elif formula == 'isotopes':
        formula = ''
        content = isotopes()
    elif formula:
        content = analyze(formula[:100])
    else:
        content = help.format(version=__version__)
    if formula:
        formula = escape(formula, True)
    return template.format(formula=formula, url=url, version=__version__,
                           content=content, heads=heads.strip())


def analyze(formula, maxatoms=250):
    """Return analysis of formula as HTML string."""

    def html(formula):
        """Return formula as HTML string."""
        formula = re.sub(
            r'\[(\d+)([A-Za-z]{1,2})\]', r'<sup>\1</sup>\2', formula)
        formula = re.sub(
            r'([A-Za-z]{1,2})(\d+)', r'\1<sub>\2</sub>', formula)
        return formula

    result = []
    try:
        f = molmass.Formula(formula)
        result.append(
            '<p><strong>Hill notation</strong>: %s</p>' % html(f.formula))
        if f.formula != f.empirical:
            result.append(
                '<p><strong>Empirical formula</strong>: %s</p>' % html(
                    f.empirical))

        prec = molmass.precision_digits(f.mass, 8)
        if f.mass != f.isotope.mass:
            result.append(
                '<p><strong>Average mass</strong>: %.*f</p>' % (prec, f.mass))
        result.extend((
            '<p><strong>Monoisotopic mass</strong>: %.*f</p>' % (
                prec, f.isotope.mass),
            '<p><strong>Nominal mass</strong>: %i</p>' % f.isotope.massnumber))

        c = f.composition()
        if len(c) > 1:
            result.extend((
                '<h3>Elemental Composition</h3>'
                '<table border="0" cellpadding="2">',
                '<tr>',
                '<th scope="col" align="left">Element</th>',
                '<th scope="col" align="right">Number</th>',
                '<th scope="col" align="right">Relative mass</th>',
                '<th scope="col" align="right">Fraction %</th>',
                '</tr>'))

            for symbol, count, mass, fraction in c:
                symbol = re.sub(r'^(\d+)(.+)', r'<sup>\1</sup>\2', symbol)
                result.extend((
                    '<tr><td>%s</td>' % symbol,
                    '<td align="center">%i</td>' % count,
                    '<td align="right">%.*f</td>' % (prec, mass),
                    '<td align="right">%.4f</td></tr>' % (fraction * 100)))

            count, mass, fraction = c.total
            result.extend((
                '<tr><td><em>Total</em></td>',
                '<td align="center"><em>%i</em></td>' % count,
                '<td align="right"><em>%.*f</em></td>' % (prec, mass),
                '<td align="right"><em>%.4f</em></td></tr>' % (fraction * 100),
                '</table>'))

        if f.atoms < maxatoms:
            s = f.spectrum()
            if len(s) > 1:
                norm = 100.0 / s.peak[1]
                result.extend((
                    '<h3>Mass Distribution</h3>',
                    '<p><strong>Most abundant mass</strong>: ',
                    '%.*f (%.3f%%)</p>' % (prec, s.peak[0], s.peak[1] * 100),
                    '<p><strong>Mean mass</strong>: %.*f</p>' % (prec, s.mean),
                    '<table border="0" cellpadding="2">',
                    '<tr>',
                    '<th scope="col" align="left">Relative mass</th>',
                    '<th scope="col" align="right">Fraction %</th>',
                    '<th scope="col" align="right">Intensity</th>',
                    '</tr>'))
                for mass, fraction in s.values():
                    result.extend((
                        '<tr><td>%.*f</td>' % (prec, mass),
                        '<td align="right">%.6f</td>' % (fraction * 100.0),
                        '<td align="right">%.6f</td></tr>' % (fraction*norm)))
                result.append('</table>')

    except Exception as e:
        e = str(e).splitlines()
        text = e[0][0].upper() + e[0][1:]
        details = '\n'.join(e[1:])
        result.append('<h2>Error</h2><p>%s</p><pre>%s</pre>' % (
            escape(text, True), escape(details, True)))

    return '\n'.join(result)


def isotopes():
    """Return HTML table of isotope masses and abundances."""
    template = """
    <h3>Isotopic Composition of the Elements</h3>
    <table border="0" cellpadding="2">
    <tr>
    <th align="left" scope="col">Element/Isotope</th>
    <th align="right" scope="col">Relative mass</th>
    <th align="right" scope="col">Abundance</th>
    </tr>
    {rows}
    </table>"""
    rows = []
    for ele in molmass.ELEMENTS:
        rows.extend((
            '<tr><td><em>%s</em></td>' % ele.name,
            '<td align="right">%12.8f</td></tr>' % ele.mass))
        for massnumber in sorted(ele.isotopes):
            iso = ele.isotopes[massnumber]
            rows.extend((
                '<tr><td align="right"><sup>%i</sup>%s</td>' % (
                    massnumber, ele.symbol),
                '<td align="right">%.8f</td>' % iso.mass,
                '<td align="right">%.8f</td></tr>' % (iso.abundance * 100.0)))
    return template.format(rows='\n'.join(rows))


def groups():
    """Return HTML table of chemical groups."""
    template = """
    <h3>{title}</h3>
    <table border="0" cellpadding="2">
    <tr>
    <th align="left" scope="col">Name</th>
    <th align="left" scope="col">Formula</th>
    </tr>
    {rows}
    </table>"""
    result = []
    for group, title in ((molmass.GROUPS, 'Chemical Groups'),
                         (molmass.AMINOACIDS, 'Amino Acids'),
                         (molmass.DEOXYNUCLEOTIDES, 'Deoxynucleotides'),
                         (molmass.NUCLEOTIDES, 'Nucleotides')):
        rows = []
        for key in sorted(group):
            value = group[key]
            if isinstance(value, str):
                rows.append('<tr><td>%s</td><td>%s</td></tr>' % (key, value))
        result.append(template.format(title=title, rows='\n'.join(rows)))
    return '\n'.join(result)


def main(url='http://localhost:9000/%s'):
    """Run web application in local web server."""
    import cgi
    import cgitb
    cgitb.enable()

    dirname, filename = os.path.split(__file__)
    if filename[-1] != 'y':
        filename = filename[:-1]  # don't use .pyc or .pyo
    url = url % filename
    if dirname:
        os.chdir(dirname)

    if os.getenv('SERVER_NAME'):
        print('Content-type: text/html\n\n')
        request = cgi.FieldStorage()
        request.get = request.getfirst
        print(response(request, url))
    else:
        import webbrowser
        if sys.version_info[0] == 2:
            from urlparse import urlparse
            from BaseHTTPServer import HTTPServer
            from CGIHTTPServer import CGIHTTPRequestHandler
        else:
            from urllib.parse import urlparse
            from http.server import HTTPServer, CGIHTTPRequestHandler

        def is_cgi(self):
            """Monkey patch for CGIHTTPRequestHandler.is_cgi()."""
            if filename in self.path:
                self.cgi_info = '', self.path[1:]
                return True
            return False

        CGIHTTPRequestHandler.is_cgi = is_cgi
        print('Serving CGI script at', url)
        webbrowser.open(url)
        url = urlparse(url)
        HTTPServer((url.hostname, url.port),
                   CGIHTTPRequestHandler).serve_forever()


if __name__ == '__main__':
    sys.exit(main())
