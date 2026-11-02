"""Build the GitHub Pages site from Markdown in the repository, so each page and its source stay one document.

  python3 docs/build_site.py      # needs pandoc

- docs/index.html from GUIDE.md: the visual guide to all the studies.
- docs/<study>/index.html from <study>/COMPANION.md, for every study folder that has one: the study's visual
  companion, with links back to the guide and to the study's paper (<study>/paper/paper.pdf on GitHub). The page
  title is the first "# " heading of COMPANION.md.

Mermaid diagrams render in the browser. Links to files in the repository point to them on GitHub; in a companion they
are relative to the study folder. Images in a companion (Markdown or <img>) with a relative path are copied into
docs/<study>/ so the Pages site serves them: a path inside the study folder keeps its place (figures/x.png is copied
to docs/<study>/figures/x.png), and a path outside it (../other/x.png) is copied to docs/<study>/_repo/other/x.png
with the src rewritten.
"""
import posixpath
import re
import shutil
import subprocess
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/rewire-bio/amp-2026-studies"

TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>@@TITLE@@</title>
<meta name="description" content="@@DESCRIPTION@@">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root {
  --bg: #fbfaf7; --surface: #ffffff; --ink: #1f2430; --muted: #5e6874; --line: #e4e2dc;
  --accent: #df007d; --accent-soft: #fde2ef; --blue: #005eb8; --blue-soft: #e3eefa;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #14161b; --surface: #1c1f26; --ink: #e8e9ec; --muted: #a3a9b3; --line: #2e323b;
    --accent: #ff4fa8; --accent-soft: #3a1830; --blue: #6aa8ff; --blue-soft: #16263d;
  }
}
:root[data-theme="dark"] {
  --bg: #14161b; --surface: #1c1f26; --ink: #e8e9ec; --muted: #a3a9b3; --line: #2e323b;
  --accent: #ff4fa8; --accent-soft: #3a1830; --blue: #6aa8ff; --blue-soft: #16263d;
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 16px/1.65 Inter, system-ui, sans-serif; }
header.hero { background: linear-gradient(135deg, #2b2f3b 0%, #3d2240 100%); color: #fff; padding: 56px 16px 48px; }
header.hero .inner { max-width: 1080px; margin: 0 auto; }
header.hero .kicker { text-transform: uppercase; letter-spacing: .12em; font-size: 13px; color: #f29ccb; font-weight: 600; }
header.hero h1 { font-size: clamp(30px, 5vw, 46px); line-height: 1.15; margin: 10px 0 14px; }
header.hero p { max-width: 720px; color: #d9dbe2; font-size: 18px; margin: 0 0 18px; }
header.hero a.btn { display: inline-block; color: #fff; border: 1px solid #ffffff55; border-radius: 999px; padding: 6px 16px;
  text-decoration: none; font-size: 14px; margin-right: 8px; }
header.hero a.btn:hover { background: #ffffff1a; }
.layout { max-width: 1080px; margin: 0 auto; padding: 0 16px; display: grid; grid-template-columns: 220px minmax(0, 1fr); gap: 40px; }
nav.toc { position: sticky; top: 24px; align-self: start; padding-top: 40px; font-size: 14px; }
nav.toc a { display: block; color: var(--muted); text-decoration: none; padding: 4px 0 4px 12px; border-left: 2px solid var(--line); }
nav.toc a:hover { color: var(--accent); border-left-color: var(--accent); }
main { padding: 24px 0 80px; min-width: 0; }
main img { max-width: 100%; height: auto; }
main h1:first-child { display: none; }
main h2 { font-size: 28px; margin: 56px 0 12px; padding-top: 8px; border-top: 1px solid var(--line); }
main h3 { font-size: 21px; margin: 36px 0 8px; color: var(--accent); }
main a { color: var(--blue); }
main blockquote { margin: 16px 0; padding: 16px 20px; background: var(--surface); border: 1px solid var(--line);
  border-left: 4px solid var(--accent); border-radius: 10px; }
main blockquote p { margin: 6px 0; }
pre.mermaid { background: var(--surface); border: 1px solid var(--line); border-radius: 12px; padding: 18px; margin: 18px 0;
  overflow-x: auto; text-align: center; }
pre.mermaid code { font: inherit; white-space: pre; }
table { border-collapse: collapse; width: 100%; font-size: 14.5px; margin: 14px 0; display: block; overflow-x: auto; }
th, td { border-bottom: 1px solid var(--line); padding: 8px 10px; text-align: left; vertical-align: top; }
th { color: var(--muted); font-weight: 600; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .9em; background: var(--blue-soft);
  padding: 1px 5px; border-radius: 4px; }
ul li, ol li { margin: 4px 0; }
footer { border-top: 1px solid var(--line); color: var(--muted); font-size: 14px; padding: 24px 16px 48px; text-align: center; }
@media (max-width: 860px) { .layout { grid-template-columns: 1fr; gap: 0; } nav.toc { display: none; } }
</style>
</head>
<body>
<header class="hero"><div class="inner">
  <div class="kicker">@@KICKER@@</div>
  <h1>@@HEADING@@</h1>
  <p>@@LEAD@@</p>
@@BUTTONS@@
</div></header>
<div class="layout">
<nav class="toc">@@TOC@@</nav>
<main>
@@BODY@@
</main>
</div>
<footer>Tim Richardson · Built from <a href="@@REPO@@/blob/main/@@SOURCE@@">@@SOURCE@@</a> · Code under the MIT licence</footer>
<script src="https://cdn.jsdelivr.net/npm/mermaid@11.4.1/dist/mermaid.min.js"></script>
<script>
  // Draw after the web font has loaded, so labels are measured in the font they are shown in.
  const dark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  mermaid.initialize({ startOnLoad: false, theme: dark ? "dark" : "default", securityLevel: "strict",
    fontFamily: "Inter, system-ui, sans-serif",
    flowchart: { htmlLabels: true, curve: "basis", wrappingWidth: 260, subGraphTitleMargin: { top: 8, bottom: 16 } } });
  document.fonts.ready.then(() => mermaid.run());
</script>
</body>
</html>
"""


GUIDE_BUTTONS = """  <a class="btn" href="@@REPO@@">Code, data and results on GitHub</a>
  <a class="btn" href="#5-reading-list-in-order">Reading list</a>"""


def render(md: Path) -> str:
    """Markdown to an HTML fragment, with mermaid blocks left as text for the browser to draw."""
    body = subprocess.run(["pandoc", "-f", "gfm", "-t", "html", "--wrap=none", str(md)],
                          check=True, capture_output=True, text=True).stdout
    # Mermaid reads the block's text; drop pandoc's <code> wrapper.
    return re.sub(r'<pre class="mermaid"><code>(.*?)</code></pre>', r'<pre class="mermaid">\1</pre>', body, flags=re.S)


def github_links(body: str, base: str = "") -> str:
    """Point relative links at the file on GitHub. base is the folder the Markdown sits in, relative to the root."""
    def url(m: re.Match) -> str:
        href = m.group(1)
        path = posixpath.normpath(posixpath.join(base, href)) if base else href
        if path.startswith("../"):
            raise SystemExit(f"link outside the repository: {href} (in {base or '.'})")
        kind = "tree" if href.endswith("/") or path == "." else "blob"
        path = "" if path == "." else path + ("/" if href.endswith("/") and not path.endswith("/") else "")
        return f'href="{REPO}/{kind}/main/{path}"'
    return re.sub(r'href="(?!https?:|#|mailto:)([^"]+)"', url, body)


def page(**fields: str) -> str:
    out = TEMPLATE
    # Body last, so placeholder names inside the Markdown are left alone.
    for key in ("TITLE", "DESCRIPTION", "KICKER", "HEADING", "LEAD", "BUTTONS", "TOC", "SOURCE", "REPO", "BODY"):
        out = out.replace(f"@@{key}@@", fields[key])
    return out


def toc(body: str) -> str:
    return "".join(f'<a href="#{i}">{re.sub("<[^>]+>", "", t)}</a>'
                   for i, t in re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body))


def build_guide() -> None:
    body = github_links(render(ROOT / "GUIDE.md"))
    html = page(
        TITLE="AMP 2026 studies: a visual guide",
        DESCRIPTION="What each study behind the AMP 2026 talk on genomic foundation models asks, why, how it works, "
                    "what it found, and what to read.",
        KICKER="AMP 2026 · Scoring variants with genomic foundation models",
        HEADING="A visual guide to the studies",
        LEAD="What each study behind the talk asks, why it was worth asking, how it works, what it found, "
             "and what to read first.",
        BUTTONS=GUIDE_BUTTONS.replace("@@REPO@@", REPO), TOC=toc(body), SOURCE="GUIDE.md", REPO=REPO, BODY=body)
    (ROOT / "docs" / "index.html").write_text(html)
    print(f"wrote docs/index.html ({len(html):,} bytes)")


def copy_images(body: str, study: str, out: Path) -> str:
    """Copy each relative image into docs/<study>/ and point its src at the copy."""
    folder = ROOT / study

    def fix(m: re.Match) -> str:
        src = m.group(2)
        if re.match(r"[a-z][a-z0-9+.-]*:|//|#", src, re.I):   # absolute URL, data: URI or anchor
            return m.group(0)
        rel = posixpath.normpath(posixpath.join(study, src.split("#")[0].split("?")[0]))
        if rel.startswith("../"):
            raise SystemExit(f"{study}/COMPANION.md: image outside the repository: {src}")
        if not (ROOT / rel).is_file():
            raise SystemExit(f"{study}/COMPANION.md: image not found: {src} (looked for {rel})")
        inside = rel.startswith(study + "/")
        new = rel[len(study) + 1:] if inside else f"_repo/{rel}"
        dest = out / new
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, dest)
        print(f"  copied {rel} to docs/{study}/{new}")
        return f'{m.group(1)}"{new}"'

    return re.sub(r'(<img\b[^>]*?\ssrc=)"([^"]+)"', fix, body)


def build_companion(md: Path) -> None:
    study = md.parent.name
    out = ROOT / "docs" / study
    out.mkdir(parents=True, exist_ok=True)
    body = render(md)
    body = copy_images(body, study, out)
    body = github_links(body, base=study)
    h1 = re.search(r"<h1[^>]*>(.*?)</h1>", body, flags=re.S)
    heading = h1.group(1).strip() if h1 else escape(study)
    title = re.sub("<[^>]+>", "", heading)
    pdf = f"{REPO}/blob/main/{study}/paper/paper.pdf"
    raw = f"{REPO}/raw/main/{study}/paper/paper.pdf"
    buttons = "\n".join([
        '  <a class="btn" href="../">Back to the guide</a>',
        f'  <a class="btn" href="{pdf}">Paper (PDF)</a>',
        f'  <a class="btn" href="{raw}">Paper (raw PDF)</a>',
        f'  <a class="btn" href="{REPO}/tree/main/{study}">Code, data and results on GitHub</a>',
    ])
    html = page(
        TITLE=f"{title}: visual companion",
        DESCRIPTION=f"Diagrams, charts and findings of the {title} study, with a link to its paper.",
        KICKER=f"AMP 2026 · Visual companion · {escape(study)}/",
        HEADING=heading,
        LEAD="The study's diagrams, charts and findings on one page. The paper has the full method and every number.",
        BUTTONS=buttons, TOC=toc(body), SOURCE=f"{study}/COMPANION.md", REPO=REPO, BODY=body)
    (out / "index.html").write_text(html)
    print(f"wrote docs/{study}/index.html ({len(html):,} bytes)")


def main() -> None:
    build_guide()
    for md in sorted(ROOT.glob("*/COMPANION.md")):
        build_companion(md)


if __name__ == "__main__":
    main()
