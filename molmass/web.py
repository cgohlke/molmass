#!/usr/bin/env python3
# molmass/web.py

# Copyright (c) 2005-2024, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Molecular Mass Calculator web application.

Run the web application in a local web server::

    python -m molmass.web

The application is run in a Flask built-in server, and a web browser is opened.

If Flask is not installed, the application is run in a built-in CGI server.
The cgi module is deprecated and slated for removal in Python 3.13.
To run in CGI mode, this script must be made executable on UNIX systems::

    chmod -x ./web.py

Do not run the built-in Flask or CGI servers in a production deployment.
Instead, for example, create a Flask app and serve it on a production server::

    # app.py
    from molmass.web import response
    from flask import Flask, request
    app = Flask(__name__)
    @app.route('/', methods=['GET'])
    def root():
        return response(request.args, url=request.base_url)

"""

from __future__ import annotations

__all__ = ['main', 'response']

import os
import re
import sys
from html import escape

try:
    from . import molmass
except ImportError:
    import molmass  # type: ignore

PAGE = """<!DOCTYPE html PUBLIC
"-//W3C//DTD XHTML 1.1 plus MathML 2.0 plus SVG 1.1//EN"
"http://www.w3.org/2002/04/xhtml-math-svg/xhtml-math-svg.dtd">
<html xmlns="http://www.w3.org/1999/xhtml"
xmlns:mathml="http://www.w3.org/1998/Math/MathML"
xmlns:svg="http://www.w3.org/2000/svg" xml:lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<style type="text/css">
html {{}}
body {{min-width: 576px, margin: 0.5em; padding: 0.25em 0.5em 0 0.5em}}
div.header {{display: flex; flex-wrap: wrap; align-items: baseline;
  padding-bottom: 0.4em}}
div.header p a {{text-decoration: none; color: #000000}}
div.header p a:hover {{text-decoration: underline; color: #E00000}}
h1 {{margin: 0; padding: 0 0.5em 0 0}}
h1 a {{text-decoration: none; color: #000000}}
h2.hidden {{display: none}}
form {{display: flex; flex-wrap: nowrap; background-color: #eeeeee;
  border: 1px solid #aaaaaa; padding: 1em}}
input {{margin-left: 0.5em; min-width: 4.5em}}
input[type=text] {{width: 100%}}
div.formula {{flex: 1 1 100%; min-width: 2em; margin-right: 1em}}
div.buttons {{white-space: nowrap}}
table {{border-spacing: 0.6em 0.3em}}
table caption {{padding: 1em 0 0.5em 0; text-align: left;
  font-weight: bold; font-size: larger}}
table th {{text-align: left}}
table.results th[scope=row] {{text-align: right}}
table.table th[scope=row] {{font-weight: normal}}
tr.spacer * {{padding-top: 0.6em}}
a:hover {{text-decoration: underline; color: #E00000}}
</style>{heads}
<meta name="generator" content="molmass" />
<meta name="robots" content="noarchive" />
<meta name="format-detection" content="telephone=no" />
<meta name="viewport" content="width=608px" />
<title>Molecular Mass Calculator v{version}</title>
</head>
<body>
<div class="header">
<h1>Molecular Mass Calculator</h1>
<p>by <a href="https://www.cgohlke.com">Christoph Gohlke</a></p>
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

HELP = """<h2 class="hidden">Help</h2>
<p>This web application calculates the molecular mass (average,
monoisotopic, and nominal), the elemental composition, and the mass
distribution spectrum of a molecule given by its chemical formula, relative
element weights, or sequence.</p>
<p>Calculations are based on the
<a href="?q=isotopes" rel="nofollow">isotopic composition of the elements</a>.
Mass deficiency due to chemical bonding is not considered.
</p>
<h3>Examples</h3>
<ul>
<li>Elements and counts:
<a href="?q=H2O" rel="nofollow">H2O</a> or
<a href="?q=C8H10N4O2" rel="nofollow">C8H10N4O2</a>
</li>
<li>Specific isotopes:
<a href="?q=D2O" rel="nofollow">D2O</a> or
<a href="?q=[30Si]3O2" rel="nofollow">[30Si]3O2</a>
</li>
<li>Ion charges:
<a href="?q=SO4_2-" rel="nofollow">SO4_2-</a> or
<a href="?q=[AsO4]3-" rel="nofollow">[AsO4]3-</a>
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
<a href="https://www.cgohlke.com">Christoph Gohlke</a>.
Source code is available on
<a href="https://github.com/cgohlke/molmass">GitHub</a>.</p>
</div>"""


def response(
    form,
    /,  # for compatibility
    url: str,
    template: str | None = None,
    help: str | None = None,
    heads: str = '',
) -> str:
    """Return HTML document from submitted web form.

    Parameters:
        form: Flask.request.form or cgi.FieldStorage.
        url: URL of web application.
        template: HTML page template. The default is :py:attr:`PAGE`.
        help: Help text. The default is :py:attr:`HELP`.
        heads: Additional HTML page head sections.

    """
    if template is None:
        template = PAGE
    if help is None:
        help = HELP

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
        content = help.format(version=molmass.__version__)
    if formula:
        formula = escape(formula, True)
    return template.format(
        formula=formula,
        url=url,
        version=molmass.__version__,
        content=content,
        heads=heads.strip(),
    )


def analyze(
    formula: str, /, *, maxatoms: int = 512, min_intensity: float = 1e-4
) -> str:
    """Return analysis of formula as HTML string.

    Parameters:
        formula: Chemical formula.
        maxatoms: Number of atoms below which to calculate spectrum.
        min_intensity: Minimum intensity to include in spectrum.

    """

    def html(formula: str) -> str:
        # return formula as HTML string
        formula = re.sub(r'\](\d*[+-]+)$', r']<sup>\1</sup>', formula)
        formula = re.sub(
            r'\[(\d+)([A-Za-z]{1,2})\]', r'<sup>\1</sup>\2', formula
        )
        formula = re.sub(r'([A-Za-z]{1,2})(\d+)', r'\1<sub>\2</sub>', formula)
        return formula

    spectrum: molmass.Spectrum | None
    result: list[str] = ['<h2 class="hidden">Results</h2>']
    try:
        f = molmass.Formula(formula)
        composition = f.composition()
        if f.atoms < maxatoms:
            spectrum = f.spectrum(min_intensity=min_intensity)
            if len(spectrum) < 2:
                spectrum = None
        else:
            spectrum = None
        result.append('<table class="results"><caption></caption>')
        result.append(
            '<tr><th scope="row">Hill notation</th>'
            f'<td colspan="2">{html(f.formula)}</td></tr>'
        )
        if f.formula != f.empirical:
            result.append(
                f'<tr><th scope="row">Empirical formula</th>'
                f'<td colspan="2">{html(f.empirical)}</td></tr>'
            )

        prec = max(
            (
                1,
                molmass.precision_digits(f.mass, 9),
                molmass.precision_digits(f.isotope.mass, 9),
            )
        )

        result.append(
            f'<tr class="spacer"><th scope="row">Nominal mass</th>'
            f'<td>{f.nominal_mass}</td></tr>'
        )
        if f.mass != f.isotope.mass:
            # formula is not an isotope
            result.append(
                f'<tr><th scope="row">Average mass</th>'
                f'<td>{f.mass:.{prec}f}</td></tr>'
            )
        result.append(
            f'<tr><th scope="row">Monoisotopic mass</th>'
            f'<td>{f.isotope.mass:.{prec}f}</td>'
            f'<td>{f.isotope.abundance * 100:.4f}%</td></tr>'
        )
        if spectrum:
            result.append(
                '<tr><th scope="row">Most abundant mass</th>'
                f'<td>{spectrum.peak.mass:.{prec}f}</td>'
                f'<td>{spectrum.peak.fraction * 100:.4f}%</td></tr>'
            )
            result.append(
                '<tr><th scope="row">Mean of distribution</th>'
                f'<td>{spectrum.mean:.{prec}}</td></tr>'
            )
        result.append(
            '<tr class="spacer"><th scope="row">Number of atoms</th>'
            f'<td>{f.atoms}</td></tr>'
        )
        result.append('</table>')

        if len(composition) > 1:
            result.extend(
                (
                    '<!--h3>Elemental Composition</h3-->',
                    '<table class="table">'
                    '<caption>Elemental Composition</caption>',
                    '<tr>',
                    '<th scope="col">Element</th>',
                    '<th scope="col">Count</th>',
                    '<th scope="col">Relative mass</th>',
                    '<th scope="col">Fraction %</th>',
                    '</tr>',
                )
            )
            for i in composition.values():
                symbol = re.sub(r'^(\d+)(.+)', r'<sup>\1</sup>\2', i.symbol)
                result.extend(
                    (
                        '<tr>',
                        f'<th scope="row">{symbol}</th>',
                        f'<td>{i.count}</td>',
                        f'<td>{i.mass:.{prec}f}</td>',
                        f'<td>{i.fraction * 100:.4f}</td>',
                        '</tr>',
                    )
                )
            result.append('</table>')

        if spectrum is not None and len(spectrum) > 1:
            if abs(spectrum._charge) > 1:
                mz = '\n<th scope="col">m/z</th>'
            else:
                mz = ''
            result.extend(
                (
                    '<!--h3>Mass Distribution</h3-->',
                    '<table class="table">',
                    '<caption>Mass Distribution</caption>',
                    '<tr>',
                    '<th scope="col">Mass number</th>',
                    '<th scope="col">Relative mass</th>',
                    '<th scope="col">Fraction %</th>',
                    '<th scope="col">Intensity %</th>',
                    mz,
                    '</tr>',
                )
            )
            for item in spectrum.values():
                if mz:
                    mz = f'<td>{item.mz:.{prec}f}</td>'
                result.extend(
                    (
                        '<tr>',
                        f'<th scope="row">{item.massnumber}</th>',
                        f'<td>{item.mass:.{prec}f}</td>',
                        f'<td>{item.fraction*100.:.6}</td>',
                        f'<td>{item.intensity:.6}</td>' f'{mz}',
                        '</tr>',
                    )
                )
            result.append('</table>')

    except Exception as exc:
        msg = str(exc).splitlines()
        text = msg[0][0].upper() + msg[0][1:]
        details = '\n'.join(msg[1:])
        result.append(
            f'<h3>Error: {escape(text, True)}</h3>'
            f'<pre>{escape(details, True)}</pre>'
        )

    return '\n'.join(i for i in result if i)


def isotopes() -> str:
    """Return HTML table of isotope masses and abundances."""
    template = """
    <!--h2>Isotopic Composition of the Elements</h2-->
    <table>
    <caption>Isotopic Composition of the Elements</caption>
    <tr>
    <th scope="col">Element/Isotope</th>
    <th scope="col" align="right">Relative mass</th>
    <th scope="col" align="right">Abundance</th>
    </tr>
    {rows}
    </table>""".strip().replace(
        '    ', ''
    )

    rows: list[str] = []
    for ele in molmass.ELEMENTS:
        rows.extend(
            (
                f'<tr><td><a href="?q={ele.symbol}">{ele.name}</a></td>',
                f'<td align="right">{ele.mass:12.8f}</td></tr>',
            )
        )
        for massnumber in sorted(ele.isotopes):
            iso = ele.isotopes[massnumber]
            rows.extend(
                (
                    '<tr>',
                    '<td align="right">'
                    f'<a href="?q={massnumber}{ele.symbol}">'
                    f'<sup>{massnumber}</sup>{ele.symbol}</a></td>',
                    f'<td align="right">{iso.mass:.8f}</td>',
                    f'<td align="right">{iso.abundance * 100.:.8f}</td>',
                    '</tr>',
                )
            )
    return template.format(rows='\n'.join(rows))


def groups() -> str:
    """Return HTML table of chemical groups."""
    template = """
    <!--h3>{title}</h3-->
    <table>
    <caption>{title}</caption>
    <tr>
    <th scope="col">Name</th>
    <th scope="col">Formula</th>
    </tr>
    {rows}
    </table>""".strip().replace(
        '    ', ''
    )

    result: list[str] = ['<!--h2>Abbreviations of chemical groups</h2-->']
    for group, title in (
        (molmass.GROUPS, 'Chemical Groups'),
        (molmass.AMINOACIDS, 'Amino Acids'),
        (molmass.DEOXYNUCLEOTIDES, 'Deoxynucleotides'),
        (molmass.NUCLEOTIDES, 'Nucleotides'),
    ):
        rows: list[str] = []
        for key in sorted(group):
            value = group[key]
            if isinstance(value, str):
                link = key if group is molmass.GROUPS else value
                rows.append(
                    f'<tr><td>{key}</td>'
                    f'<td><a href="?q={link}">{value}</a></td></tr>'
                )
        result.append(template.format(title=title, rows='\n'.join(rows)))
    return '\n'.join(result)


def webbrowser(url: str, /, delay: float = 1.0) -> None:
    """Open url in web browser after delay.

    Parameters:
        url: URL to open in web browser.
        delay: Delay in seconds before opening web browser.

    """
    import threading
    import webbrowser

    threading.Timer(delay, lambda: webbrowser.open(url)).start()


def cgi(url: str, *, open_browser: bool = True, debug: bool = True) -> int:
    """Run web application in local CGI server.

    Parameters:
        url: URL at which the web application is served.
        open_browser: Open `url` in web browser.
        debug: Enable debug mode.

    """
    import cgi
    import cgitb

    if debug:
        cgitb.enable()

    dirname, filename = os.path.split(__file__)
    if filename[-1] != 'y':
        filename = filename[:-1]  # don't use .pyc or .pyo
    if not url.endswith('/'):
        url += '/'
    url += filename
    if dirname:
        os.chdir(dirname)

    if os.getenv('SERVER_NAME'):
        print('Content-type: text/html\n\n')
        request = cgi.FieldStorage()
        request.get = request.getfirst  # type: ignore
        print(response(request, url))
    else:
        from http.server import CGIHTTPRequestHandler, HTTPServer
        from urllib.parse import urlparse

        def is_cgi(self) -> bool:
            # monkey patch for CGIHTTPRequestHandler.is_cgi
            if filename in self.path:
                self.cgi_info = '', self.path[1:]
                return True
            return False

        CGIHTTPRequestHandler.is_cgi = is_cgi  # type: ignore
        print('Running CGI script at', url)
        if open_browser:
            webbrowser(url)
        urlp = urlparse(url)
        if urlp.hostname is None or urlp.port is None:
            raise ValueError(f'invalid URL {url!r}')
        HTTPServer(
            (urlp.hostname, urlp.port), CGIHTTPRequestHandler
        ).serve_forever()
    return 0


def main(
    *,
    form: dict[str, str] | None = None,
    url: str | None = None,
    open_browser: bool = True,
    debug: bool = False,
) -> int:
    """Run web application in local Flask or fallback CGI server.

    Parameters:
        form:
            Form data to use on first run.
        url:
            URL at which the web application is served.
            The default is 'http://127.0.0.1:5001/'.
        open_browser:
            Open `url` in web browser.
        debug:
            Enable debug mode.

    """
    if url is None:
        url = 'http://127.0.0.1:5001/'

    try:
        from flask import Flask, request
    except ImportError:
        return cgi(url, open_browser=open_browser, debug=debug)

    from urllib.parse import urlparse

    urlp = urlparse(url)
    if urlp.hostname is None or urlp.port is None:
        raise ValueError(f'invalid URL {url!r}')

    app = Flask(__name__)

    @app.route('/', methods=['GET'])
    def root() -> str:
        nonlocal form
        r = response(
            request.args if form is None else form, url=request.base_url
        )
        form = None
        return r

    if open_browser:
        webbrowser(url)

    app.run(host=urlp.hostname, port=urlp.port, debug=debug)
    return 0


if __name__ == '__main__':
    sys.exit(main())
