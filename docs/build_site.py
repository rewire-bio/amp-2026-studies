"""Build the GitHub Pages explainer (docs/index.html) from GUIDE.md, so the page and the guide stay one document.

  python3 docs/build_site.py      # needs pandoc

Mermaid diagrams render in the browser. Links to files in the repository point to them on GitHub.
"""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/rewire-bio/amp-2026-studies"

TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AMP 2026 studies: a visual guide</title>
<meta name="description" content="What each study behind the AMP 2026 talk on genomic foundation models asks, why, how it works, what it found, and what to read.">
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
  <div class="kicker">AMP 2026 · Scoring variants with genomic foundation models</div>
  <h1>A visual guide to the studies</h1>
  <p>What each study behind the talk asks, why it was worth asking, how it works, what it found, and what to read first.</p>
  <a class="btn" href="REPO">Code, data and results on GitHub</a>
  <a class="btn" href="#5-reading-list-in-order">Reading list</a>
</div></header>
<div class="layout">
<nav class="toc">TOC</nav>
<main>
BODY
</main>
</div>
<footer>Tim Richardson · Built from <a href="REPO/blob/main/GUIDE.md">GUIDE.md</a> · Code under the MIT licence</footer>
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


def main() -> None:
    body = subprocess.run(["pandoc", "-f", "gfm", "-t", "html", "--wrap=none", str(ROOT / "GUIDE.md")],
                          check=True, capture_output=True, text=True).stdout
    # Mermaid reads the block's text; drop pandoc's <code> wrapper.
    body = re.sub(r'<pre class="mermaid"><code>(.*?)</code></pre>', r'<pre class="mermaid">\1</pre>', body, flags=re.S)
    # Links to files in the repository open them on GitHub.
    body = re.sub(r'href="(?!https?:|#|mailto:)([^"]+)"',
                  lambda m: f'href="{REPO}/{"tree" if m.group(1).endswith("/") else "blob"}/main/{m.group(1)}"', body)
    toc = "".join(f'<a href="#{i}">{re.sub("<[^>]+>", "", t)}</a>'
                  for i, t in re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body))
    # Body last, so placeholder names inside the guide's text (PROTOCOL.md contains "TOC") are left alone.
    page = TEMPLATE.replace("TOC", toc).replace("REPO", REPO).replace("BODY", body)
    (ROOT / "docs" / "index.html").write_text(page)
    print(f"wrote docs/index.html ({len(page):,} bytes)")


if __name__ == "__main__":
    main()
