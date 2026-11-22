# Inverse Problems submission typesetting

The unmodified publisher class and its optional ORCID icon come from the July 2025 IOP template archive. The class is licensed under LPPL 1.3c or later. The publisher also accepts a standard LaTeX article layout. [Official IOP template archive](https://publishingsupport.iopscience.iop.org/wp-content/uploads/2025/07/ioplatextemplate.zip) · [IOP author instructions](https://publishingsupport.iopscience.iop.org/publishing-support/authors/authoring-for-journals/writing-journal-article/)

Tectonic 0.17.0 on the local Mac compiled both the official example and an anonymous class check. From the submission directory, with the manuscript source and the figures it names available, run

```sh
tectonic -X compile -Z search-path=typeset --keep-logs manuscript.tex
pdfinfo manuscript.pdf
pdftotext -layout manuscript.pdf -
```

Tectonic resolves standard TeX packages on first use. The same source can be built by latexmk or pdfLaTeX if the class and figure files are on the respective search paths. Use `\documentclass[12pt,anonymous]{iopjournal}` for double-anonymous review. The class changes `\title` into a typesetting command inside the document and uses class-provided `\author`, `\affil`, and `\email` commands in place of `\maketitle`. If you omit the template's journal issue placeholders, use `\pagestyle{plain}` and omit `\articletype`. The class loads hyperref and graphicx. Set the final review title and blank PDF author via `\hypersetup{pdftitle={...},pdfauthor={}}`.

Inverse Problems lets authors choose single- or double-anonymous review. Initial submission requires the complete manuscript PDF with the figures and tables in the text, plus relevant supplements. For double-anonymous review, check the PDF properties, filenames, manuscript, justification letter, supplementary files, and any repository or data links for identifying information. Enter author and funding details in the submission system. A direct link to an author-named project discloses identity to reviewers. The anonymous class replaces its special acknowledgments, funding, authorship and data-availability sections with anonymous stubs, but it does not scan arbitrary text, external figures, supplementary files or links. The editor can add authorship details after the review as appropriate. [Inverse Problems peer review](https://publishingsupport.iopscience.iop.org/journals/inverse-problems/about-inverse-problems/) · [IOP anonymisation checklist](https://publishingsupport.iopscience.iop.org/questions/checklist-for-anonymising-your-manuscript/) · [IOP submission checklist](https://publishingsupport.iopscience.iop.org/submission-checklist/)

IOP asks for a readable manuscript, suggests body text of at least 12 pt, and adds submission line numbers automatically. TeX source for later production should compile with class, graphics and bibliography in the same directory, using filenames containing only ASCII letters, digits, and underscores before the extensions. [IOP author instructions](https://publishingsupport.iopscience.iop.org/publishing-support/authors/authoring-for-journals/writing-journal-article/)
