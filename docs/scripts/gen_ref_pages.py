"""Generate the code reference pages and navigation."""

import pathlib

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()
src = pathlib.Path(__file__).parent.parent.parent / "ties"

mod_symbol = '<code class="doc-symbol doc-symbol-nav doc-symbol-module"></code>'

for path in sorted(src.rglob("*.py")):
    if "testing" in str(path) or "scripts" in str(path) or "examples" in str(path):
        continue

    module_path = path.relative_to(src.parent).with_suffix("")
    doc_path = path.relative_to(src).with_suffix(".md")
    full_doc_path = pathlib.Path("reference", doc_path)

    parts = tuple(module_path.parts)

    if parts[-1] == "__init__":
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
        full_doc_path = full_doc_path.with_name("index.md")
    elif parts[-1].startswith("_"):
        continue

    nav_parts = [f"{mod_symbol} {part}" for part in parts]
    nav[tuple(nav_parts)] = doc_path.as_posix()

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        ident = ".".join(parts)
        fd.write(f"::: {ident}")

    mkdocs_gen_files.set_edit_path(full_doc_path, ".." / path)

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())