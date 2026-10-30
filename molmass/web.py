#!/usr/bin/env python3
# molmass/web.py

# Copyright (c) 2005-2026, Christoph Gohlke
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

The application is run in a local server based on the http.server module,
and a web browser is opened.

Do not run the built-in server in a production deployment.
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

__all__ = ['HELP', 'analyze', 'groups', 'isotopes', 'main', 'response']

import re
import sys
from html import escape
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from http.server import ThreadingHTTPServer
    from typing import Any

try:
    from . import molmass
    from .elements import ELEMENTS
except ImportError:
    from elements import ELEMENTS  # type: ignore[no-redef]

    import molmass  # type: ignore[no-redef]

DEBUG = False

MAXLENGTH = 2000
"""Maximum length of formula string accepted by :py:func:`analyze`."""

PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="author" content="Christoph Gohlke">
<meta name="generator" content="molmass.py">
<meta name="robots" content="noarchive">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<style>
body {{min-width: 455px; margin: 0.5em; padding: 0.25em 0.5em 0 0.5em;}}
header {{display: flex; flex-wrap: wrap; align-items: baseline;
 padding-bottom: 0.4em;}}
header p a {{text-decoration: none; color: #000000;}}
header p a:hover {{text-decoration: underline; color: #E00000;}}
h1 {{margin: 0; padding: 0 0.5em 0 0;}}
h1 a {{text-decoration: none; color: #000000;}}
h2.hidden {{display: none;}}
form {{display: flex; flex-wrap: nowrap; background-color: #eeeeee;
 border: 1px solid #aaaaaa; padding: 1em;}}
input {{margin-left: 0.5em; min-width: 4.5em;}}
input[type=text] {{width: 100%;}}
div.formula {{flex: 1 1 100%; min-width: 2em; margin-right: 1em;}}
div.buttons {{white-space: nowrap;}}
table {{border-spacing: 0.6em 0.3em;}}
table caption {{padding: 1em 0 0.5em 0; text-align: left;
 font-weight: bold; font-size: larger;}}
table th {{text-align: left;}}
table.results th[scope=row] {{text-align: right;}}
table.table th[scope=row] {{font-weight: normal;}}
th.right {{text-align: right;}}
td.right {{text-align: right;}}
tr.spacer * {{padding-top: 0.6em;}}
a:hover {{text-decoration: underline; color: #E00000;}}
</style>{heads}
<title>Molecular Mass Calculator v{version}</title>
</head>
<body>
<header>
<h1>Molecular Mass Calculator</h1>
<p>by <a href="https://www.cgohlke.com">Christoph Gohlke</a></p>
</header>
<main>
<form id="molmass" method="get">
<div><label for="q"><strong>Formula:</strong></label></div>
<div class="formula">
<input type="text" name="q" id="q" value="{formula}">
</div>
<div class="buttons">
<input type="submit" id="a" value="Submit">
<input type="reset" value="Reset" data-url="{url}"
 onclick="window.location=this.dataset.url">
</div>
</form>
<div class="content">
{content}
</div>
</main>
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
<a href="?q=%5B30Si%5D3O2" rel="nofollow">[30Si]3O2</a>
</li>
<li>Ion charges:
<a href="?q=SO4_2-" rel="nofollow">SO4_2-</a> or
<a href="?q=%5BAsO4%5D3-" rel="nofollow">[AsO4]3-</a>
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
<a href="?q=O:%200.26,%2030Si:%200.74" rel="nofollow">O: 0.26, 30Si: 0.74</a>
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
    form: Any,
    /,  # for compatibility
    url: str,
    template: str | None = None,
    help: str | None = None,  # noqa: A002
    heads: str = '',
) -> str:
    """Return HTML document from submitted web form.

    Parameters:
        form:
            Mapping of form fields, for example, Flask ``request.args``.
        url:
            URL of web application.
        template:
            HTML page template. The default is :py:attr:`PAGE`.
        help:
            Help text. The default is :py:attr:`HELP`.
        heads:
            Additional HTML page head sections.

    """
    if template is None:
        template = PAGE
    if help is None:
        help = HELP  # noqa: A001

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
        content = analyze(formula)
    else:
        content = help.format(version=molmass.__version__)
    if formula:
        formula = escape(formula, quote=True)
    return template.format(
        formula=formula,
        url=escape(url, quote=True),
        version=molmass.__version__,
        content=content,
        heads=heads.strip(),
    )


def analyze(
    formula: str, /, *, maxatoms: int = 512, min_intensity: float = 1e-4
) -> str:
    """Return analysis of formula as HTML string.

    Parameters:
        formula:
            Chemical formula of at most :py:attr:`MAXLENGTH` characters.
        maxatoms:
            Number of atoms below which to calculate spectrum.
        min_intensity:
            Minimum intensity to include in spectrum.

    """

    def html(formula: str) -> str:
        # return formula as HTML string
        formula = re.sub(r'\](\d*[+-]+)$', r']<sup>\1</sup>', formula)
        formula = re.sub(
            r'\[(\d+)([A-Za-z]{1,2})\]', r'<sup>\1</sup>\2', formula
        )
        return re.sub(r'([A-Za-z]{1,2})(\d+)', r'\1<sub>\2</sub>', formula)

    spectrum: molmass.Spectrum | None
    result: list[str] = ['<h2 class="hidden">Results</h2>']
    if len(formula) > MAXLENGTH:
        result.append(
            f'<h3>Error: Formula too long, maximum is {MAXLENGTH} characters'
            '</h3>'
        )
        return '\n'.join(result)
    try:
        f = molmass.Formula(formula)
        composition = f.composition()
        spectrum = None
        if f.atoms < maxatoms:
            spectrum = f.spectrum(min_intensity=min_intensity)
            if len(spectrum) == 0:
                spectrum = None
        result.append('<table class="results"><caption></caption>')
        result.append(
            '<tr><th scope="row">Hill notation</th>'
            f'<td>{html(f.formula)}</td><td></td></tr>'
        )
        if f.formula != f.empirical:
            result.append(
                f'<tr><th scope="row">Empirical formula</th>'
                f'<td>{html(f.empirical)}</td><td></td></tr>'
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
            f'<td>{f.nominal_mass}</td><td></td></tr>'
        )
        if f.mass != f.isotope.mass:
            # formula is not an isotope
            result.append(
                f'<tr><th scope="row">Average mass</th>'
                f'<td>{f.mass:.{prec}f}</td><td></td></tr>'
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
                f'<td>{spectrum.mean:.{prec}f}</td><td></td></tr>'
            )
        result.append(
            '<tr class="spacer"><th scope="row">Number of atoms</th>'
            f'<td>{f.atoms}</td><td></td></tr>'
        )
        result.append('</table>')

        if len(composition) > 1:
            result.extend(
                (
                    '<!--h3>Elemental Composition</h3-->',
                    (
                        '<table class="table">'
                        '<caption>Elemental Composition</caption>'
                    ),
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
            mz = '\n<th scope="col">m/z</th>' if abs(f.charge) > 1 else ''
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
                        f'<td>{item.fraction * 100:.6f}</td>',
                        f'<td>{item.intensity:.6f}</td>{mz}',
                        '</tr>',
                    )
                )
            result.append('</table>')

    except Exception as exc:
        if DEBUG:
            raise
        lines = str(exc).splitlines() or [type(exc).__name__]
        text = lines[0][:1].upper() + lines[0][1:]
        details = '\n'.join(lines[1:])
        result.append(
            f'<h3>Error: {escape(text, quote=True)}</h3>'
            f'<pre>{escape(details, quote=True)}</pre>'
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
    <th scope="col" class="right">Relative mass</th>
    <th scope="col" class="right">Abundance %</th>
    </tr>
    {rows}
    </table>""".strip().replace('    ', '')

    rows: list[str] = []
    for ele in ELEMENTS:
        rows.extend(
            (
                f'<tr><td><a href="?q={ele.symbol}">{ele.name}</a></td>',
                f'<td class="right">{ele.mass:12.8f}</td><td></td></tr>',
            )
        )
        for massnumber in sorted(ele.isotopes):
            iso = ele.isotopes[massnumber]
            rows.extend(
                (
                    '<tr>',
                    (
                        '<td class="right">'
                        f'<a href="?q={massnumber}{ele.symbol}">'
                        f'<sup>{massnumber}</sup>{ele.symbol}</a></td>'
                    ),
                    f'<td class="right">{iso.mass:.8f}</td>',
                    f'<td class="right">{iso.abundance * 100.:.8f}</td>',
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
    </table>""".strip().replace('    ', '')

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
        url:
            URL to open in web browser.
        delay:
            Delay in seconds before opening web browser.

    """
    import threading
    import webbrowser

    threading.Timer(delay, lambda: webbrowser.open(url)).start()


def serve(
    url: str,
    *,
    form: dict[str, str] | None = None,
    open_browser: bool = True,
) -> int:
    """Run web application in local http.server.

    Parameters:
        url:
            URL at which the web application is served.
        form:
            Form data to use on first request.
        open_browser:
            Open `url` in web browser.

    """
    server = http_server(url, form)
    print('Running web application at', url)  # noqa: T201
    if open_browser:
        webbrowser(url)
    with server:
        server.serve_forever()
    return 0


def http_server(
    url: str, form: dict[str, str] | None = None, /
) -> ThreadingHTTPServer:
    """Return http.server serving web application at root of url."""
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    from urllib.parse import parse_qs, urlparse

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            nonlocal form
            request = urlparse(self.path)
            if request.path != '/':
                self.send_error(404)
                return
            if form is None:
                query = {k: v[0] for k, v in parse_qs(request.query).items()}
            else:
                query, form = form, None
            body = response(query, url=url).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    return ThreadingHTTPServer(host_port(url), Handler)


def host_port(url: str, /) -> tuple[str, int]:
    """Return host name and port of url."""
    from urllib.parse import urlparse

    urlp = urlparse(url)
    if urlp.hostname is None or urlp.port is None:
        msg = f'invalid URL {url!r}'
        raise ValueError(msg)
    return urlp.hostname, urlp.port


def main(
    *,
    form: dict[str, str] | None = None,
    url: str | None = None,
    open_browser: bool = True,
) -> int:
    """Run web application in local http.server.

    Parameters:
        form:
            Form data to use on first request.
        url:
            URL at which the web application is served.
            The default is 'http://127.0.0.1:5001/'.
        open_browser:
            Open `url` in web browser.

    """
    if url is None:
        url = 'http://127.0.0.1:5001/'
    return serve(url, form=form, open_browser=open_browser)


if __name__ == '__main__':
    sys.exit(main())
