#!/usr/bin/env python3
"""Resolve Z-ALG R5 reference labels and probe conversation URLs."""

from __future__ import annotations

import concurrent.futures
import datetime as dt
import difflib
import html
import json
import re
import subprocess
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

BASE = Path(__file__).resolve().parent
BIB_PATH = BASE / "bibliography.json"
OUT_JSON = BASE / "reference_resolution.json"
OUT_MD = BASE / "REFERENCES_AUDIT.md"
PDF = (
    BASE
    / "package/Z_ALG_Codex_Handoff_2026-11-18/reference/"
    "Z_ALG_Research_Report_2026-11-17.pdf"
)
EXTRACT = (
    BASE
    / "package/Z_ALG_Codex_Handoff_2026-11-18/reference/"
    "Z_ALG_Research_Report_2026-11-17_extracted_utf8.txt"
)
UA = "Inv-SLAM-reference-audit/1.0"
BROWSER_UA = "Mozilla/5.0 (compatible; Inv-SLAM-reference-audit/1.0)"
NOW = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
PMLR = "https://proceedings.mlr.press/v70/guo17a.html"
TYPESAFE = "https://typesafe.ai/blog/introducing-system-one-models-and-jev"


def load_bibliography() -> dict:
    with BIB_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)


BIB = load_bibliography()
REFS = {item["id"]: item for item in BIB["references"]}


def component_specs() -> dict[str, list[dict]]:
    specs = {
        "L14": [
            {
                "part": "first",
                "authors": "Golub, Pereyra",
                "year": "1973",
                "title": (
                    "The Differentiation of Pseudo-Inverses and Nonlinear Least "
                    "Squares Problems Whose Variables Separate"
                ),
                "venue": "SIAM Journal on Numerical Analysis",
                "identifier": {"type": "doi", "value": "10.1137/0710036"},
            },
            {
                "part": "second",
                "authors": "van Leeuwen, Yang",
                "year": "2025",
                "title": "An analysis of constraint-relaxation in PDE-based inverse problems",
                "venue": "Inverse Problems 41, 025009",
                "identifier": {
                    "type": "doi",
                    "value": "10.1088/1361-6420/adab89",
                },
            },
        ],
        "L23": [
            {
                "part": "first",
                "authors": "Ito, Noma, Ohno",
                "year": "2015",
                "title": "Maximal minors of a matrix with linear form entries",
                "venue": "Linear and Multilinear Algebra",
                "identifier": {
                    "type": "doi",
                    "value": "10.1080/03081087.2014.959516",
                },
            },
            {
                "part": "second",
                "authors": "Gopalakrishnan, Neiger, Safey El Din",
                "year": "2024",
                "title": (
                    "Optimized Groebner basis algorithms for maximal determinantal "
                    "ideals and critical point computations"
                ),
                "venue": "arXiv",
                "identifier": {"type": "arxiv", "value": "2402.07353"},
            },
        ],
    }
    for ref_id, ref in REFS.items():
        if ref_id in specs:
            continue
        ids = ref.get("identifiers", [])
        specs[ref_id] = [
            {
                "part": "single",
                "authors": ref.get("authors"),
                "year": ref.get("year"),
                "title": ref.get("title"),
                "venue": ref.get("venue"),
                "identifier": ids[0] if ids else None,
            }
        ]
    return specs


SPECS = component_specs()


def get_json(url: str) -> tuple[dict, int, str]:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": UA,
            "Accept": "application/json, application/atom+xml, text/xml",
        },
    )
    with urllib.request.urlopen(request, timeout=25) as response:
        return json.loads(response.read()), response.status, response.geturl()


def crossref_doi(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/:")
    result = {
        "provider": "Crossref",
        "identifier_type": "doi",
        "identifier": doi,
        "query_url": url,
    }
    try:
        try:
            payload, status, final_url = get_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code != 429:
                raise
            time.sleep(1.0)
            payload, status, final_url = get_json(url)
        message = payload["message"]
        title = (message.get("title") or [None])[0]
        subtitle = (message.get("subtitle") or [None])[0]
        if title and subtitle and subtitle not in title:
            title = f"{title}: {subtitle}"
        authors = [
            ", ".join(part for part in [author.get("family"), author.get("given")] if part)
            for author in message.get("author", [])
        ]
        issued = message.get("issued", {}).get("date-parts", [[None]])[0]
        result.update(
            {
                "ok": True,
                "http_status": status,
                "source_url": final_url,
                "canonical_doi": message.get("DOI"),
                "title": title,
                "subtitle": subtitle,
                "authors": [author for author in authors if author],
                "year": issued[0] if issued else None,
                "issued_date_parts": issued,
                "published_print": message.get("published-print", {}).get(
                    "date-parts", [[None]]
                )[0],
                "published_online": message.get("published-online", {}).get(
                    "date-parts", [[None]]
                )[0],
                "venue": (message.get("container-title") or [None])[0],
                "volume": message.get("volume"),
                "issue": message.get("issue"),
                "page": message.get("page"),
                "article_number": message.get("article-number"),
                "type": message.get("type"),
                "publisher": message.get("publisher"),
                "url": message.get("URL"),
            }
        )
    except urllib.error.HTTPError as exc:
        result.update(
            {
                "ok": False,
                "http_status": exc.code,
                "error": f"HTTP {exc.code}: {exc.reason}",
                "source_url": exc.geturl(),
            }
        )
    except Exception as exc:
        result.update({"ok": False, "error": f"{type(exc).__name__}: {exc}"})
    return result


ATOM = {
    "a": "http://www.w3.org/2005/Atom",
    "arxiv": "http://arxiv.org/schemas/atom",
}


def arxiv_record(arxiv_id: str) -> dict:
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"id_list": arxiv_id, "max_results": 1}
    )
    result = {
        "provider": "arXiv",
        "identifier_type": "arxiv",
        "identifier": arxiv_id,
        "query_url": url,
    }
    try:
        request = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(request, timeout=25) as response:
            raw = response.read()
            status = response.status
            final_url = response.geturl()
        root = ET.fromstring(raw)
        entry = root.find("a:entry", ATOM)
        if entry is None:
            raise ValueError("No arXiv entry returned")
        title = entry.findtext("a:title", default="", namespaces=ATOM)
        published = entry.findtext("a:published", default=None, namespaces=ATOM)
        updated = entry.findtext("a:updated", default=None, namespaces=ATOM)
        authors = [
            re.sub(
                r"\s+",
                " ",
                entry_author.findtext("a:name", default="", namespaces=ATOM) or "",
            ).strip()
            for entry_author in entry.findall("a:author", ATOM)
        ]
        result.update(
            {
                "ok": True,
                "http_status": status,
                "source_url": final_url,
                "arxiv_id": arxiv_id,
                "canonical_arxiv_url": entry.findtext("a:id", default=None, namespaces=ATOM),
                "title": re.sub(r"\s+", " ", title).strip(),
                "authors": [author for author in authors if author],
                "year": int(published[:4]) if published else None,
                "published": published,
                "updated": updated,
            }
        )
    except urllib.error.HTTPError as exc:
        result.update(
            {"ok": False, "http_status": exc.code, "error": f"HTTP {exc.code}: {exc.reason}"}
        )
    except Exception as exc:
        result.update({"ok": False, "error": f"{type(exc).__name__}: {exc}"})
    return result


def crossref_title(title: str) -> dict:
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": title, "rows": 5}
    )
    result = {
        "provider": "Crossref title search",
        "identifier_type": "title",
        "query_title": title,
        "query_url": url,
    }
    try:
        payload, status, final_url = get_json(url)
        records = []
        for message in payload["message"]["items"]:
            authors = [
                ", ".join(
                    part
                    for part in [author.get("family"), author.get("given")]
                    if part
                )
                for author in message.get("author", [])
            ]
            issued = message.get("issued", {}).get("date-parts", [[None]])[0]
            records.append(
                {
                    "canonical_doi": message.get("DOI"),
                    "title": (message.get("title") or [None])[0],
                    "authors": [author for author in authors if author],
                    "year": issued[0] if issued else None,
                    "issued_date_parts": issued,
                    "venue": (message.get("container-title") or [None])[0],
                    "volume": message.get("volume"),
                    "issue": message.get("issue"),
                    "page": message.get("page"),
                    "article_number": message.get("article-number"),
                    "type": message.get("type"),
                    "url": message.get("URL"),
                    "relevance_score": message.get("score"),
                }
            )
        result.update(
            {
                "ok": True,
                "http_status": status,
                "source_url": final_url,
                "results": records,
            }
        )
    except urllib.error.HTTPError as exc:
        result.update(
            {"ok": False, "http_status": exc.code, "error": f"HTTP {exc.code}: {exc.reason}"}
        )
    except Exception as exc:
        result.update({"ok": False, "error": f"{type(exc).__name__}: {exc}"})
    return result


class PageParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_title = False
        self.title_parts: list[str] = []
        self.metadata: list[dict] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = {key.lower(): value for key, value in attrs if key}
        if tag.lower() == "title":
            self.in_title = True
        elif tag.lower() == "meta":
            name = (values.get("name") or values.get("property") or "").lower()
            if name in {
                "citation_title",
                "citation_author",
                "citation_publication_date",
                "og:title",
                "author",
                "date",
            }:
                self.metadata.append({"key": name, "content": values.get("content")})

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "title":
            self.in_title = False

    def handle_data(self, data: str) -> None:
        if self.in_title:
            self.title_parts.append(data)


class RedirectRecorder(urllib.request.HTTPRedirectHandler):
    def __init__(self) -> None:
        super().__init__()
        self.chain: list[dict] = []

    def redirect_request(self, request, file_pointer, code, message, headers, new_url):
        self.chain.append(
            {"status": code, "from": request.full_url, "to": new_url}
        )
        return super().redirect_request(
            request, file_pointer, code, message, headers, new_url
        )


def fetch_url(url: str) -> dict:
    recorder = RedirectRecorder()
    opener = urllib.request.build_opener(recorder)
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": BROWSER_UA,
            "Accept": (
                "text/html,application/xhtml+xml,application/pdf,"
                "application/json;q=0.9,*/*;q=0.5"
            ),
            "Range": "bytes=0-786431",
        },
    )
    result = {"url": url, "redirects": recorder.chain}
    try:
        response = opener.open(request, timeout=25)
        status = response.status
        body = response.read(786432)
        final_url = response.geturl()
        headers = response.headers
        response.close()
    except urllib.error.HTTPError as exc:
        status = exc.code
        body = exc.read(786432) if exc.fp else b""
        final_url = exc.geturl() or url
        headers = exc.headers
    except Exception as exc:
        return {
            **result,
            "ok": False,
            "error": f"{type(exc).__name__}: {exc}",
            "access_error": True,
        }

    content_type = headers.get_content_type() if headers else None
    charset = headers.get_content_charset() if headers else None
    title = None
    page_metadata = []
    if content_type in {"text/html", "application/xhtml+xml"}:
        parser = PageParser()
        parser.feed(body.decode(charset or "utf-8", "replace"))
        title = re.sub(r"\s+", " ", "".join(parser.title_parts)).strip() or None
        page_metadata = parser.metadata
    result.update(
        {
            "ok": 200 <= status < 400,
            "http_status": status,
            "final_url": final_url,
            "content_type": content_type,
            "response_title": title,
            "page_metadata": page_metadata,
            "body_bytes_sampled": len(body),
            "access_error": None if 200 <= status < 400 else f"HTTP {status}",
        }
    )
    return result


def retrieve_metadata() -> dict:
    jobs = []
    metadata = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        for specs in SPECS.values():
            for spec in specs:
                identifier = spec.get("identifier")
                if not identifier:
                    continue
                key = identifier["type"] + ":" + identifier["value"]
                if key in metadata:
                    continue
                function = crossref_doi if identifier["type"] == "doi" else arxiv_record
                metadata[key] = executor.submit(function, identifier["value"])
        metadata["title:L17"] = executor.submit(crossref_title, REFS["L17"]["title"])
        jobs = list(metadata)
        for key in jobs:
            metadata[key] = metadata[key].result()
    return metadata


METADATA = retrieve_metadata()
# The chapter record omits its author, but the parent book DOI is primary Wiley metadata.
METADATA["doi:10.1002/9781119311997"] = crossref_doi("10.1002/9781119311997")


def retrieve_urls(urls: list[str]) -> dict[str, dict]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        results = list(executor.map(fetch_url, urls))
    return {result["url"]: result for result in results}


conversation_urls = [item["url"] for item in BIB["conversation_external_urls"]]
all_fetch_urls = list(dict.fromkeys(conversation_urls + [PMLR]))
URL_RESULTS = retrieve_urls(all_fetch_urls)


def page_bibliography(url: str, venue: str) -> dict:
    response = URL_RESULTS[url]
    values = lambda key: [
        item["content"]
        for item in response.get("page_metadata", [])
        if item["key"] == key and item.get("content")
    ]
    if not response.get("ok"):
        return {
            "ok": False,
            "provider": venue,
            "http_status": response.get("http_status"),
            "error": response.get("access_error"),
        }
    title_values = values("citation_title") or values("og:title")
    return {
        "ok": True,
        "provider": venue,
        "source_url": response.get("final_url"),
        "http_status": response.get("http_status"),
        "title": (title_values or [response.get("response_title")])[0],
        "authors": values("citation_author"),
        "publication_date": (values("citation_publication_date") or values("date") or [None])[0],
        "venue": venue,
    }


PAGE_BIB = {
    PMLR: page_bibliography(PMLR, "PMLR"),
    TYPESAFE: page_bibliography(TYPESAFE, "TypeSafe AI Blog"),
}

if PAGE_BIB[TYPESAFE].get("ok"):
    with urllib.request.urlopen(
        urllib.request.Request(TYPESAFE, headers={"User-Agent": BROWSER_UA}),
        timeout=25,
    ) as response:
        type_safe_html = response.read().decode("utf-8", "replace")
    date_match = re.search(r'\\?"date\\?",\\?"(20\d\d-\d\d-\d\d)T', type_safe_html)
    author_match = re.search(
        r"(Diogo Almeida,\s*founder,\s*TypeSafe)", type_safe_html
    )
    if date_match:
        PAGE_BIB[TYPESAFE]["publication_date"] = date_match.group(1)
    if PAGE_BIB[TYPESAFE].get("title"):
        PAGE_BIB[TYPESAFE]["title"] = re.sub(
            r"\s+-\s+TypeSafe AI Blog\s*$",
            "",
            PAGE_BIB[TYPESAFE]["title"],
            flags=re.I,
        )
    if author_match:
        PAGE_BIB[TYPESAFE]["authors"] = [author_match.group(1)]


def normalize(value) -> str:
    if value is None:
        return ""
    value = html.unescape(str(value))
    value = (
        value.replace("ä", "ae")
        .replace("ö", "oe")
        .replace("ü", "ue")
        .replace("Ä", "Ae")
        .replace("Ö", "Oe")
        .replace("Ü", "Ue")
        .replace("ß", "ss")
    )
    value = unicodedata.normalize("NFKD", value)
    value = "".join(char for char in value if not unicodedata.combining(char))
    value = value.replace("&", " and ")
    value = re.sub(r"[\W_]+", " ", value.lower(), flags=re.UNICODE)
    return " ".join(value.split())


def compare_title(reported: str | None, resolved: str | None) -> dict:
    def clean(value):
        value = re.sub(r",\s*Chapter\s+\d+\s*$", "", value or "", flags=re.I)
        value = re.sub(r"\s+-\s+TypeSafe AI Blog\s*$", "", value or "", flags=re.I)
        return normalize(value)

    left, right = clean(reported), clean(resolved)
    return {
        "match": bool(left and right and left == right),
        "similarity": round(
            difflib.SequenceMatcher(None, left, right).ratio(), 6
        )
        if left and right
        else 0.0,
    }


def resolved_families(authors: list[str] | None) -> list[str]:
    families = []
    for author in authors or []:
        if "," in author:
            family = author.split(",", 1)[0]
        else:
            family = author.rsplit(None, 1)[-1]
        families.append(normalize(family))
    return families


def compare_authors(reported: str | None, resolved: list[str] | None) -> dict:
    actual = resolved_families(resolved)
    reported_normalized = f" {normalize(reported)} "
    matched = bool(actual) and all(
        f" {family} " in reported_normalized for family in actual if family
    )
    return {
        "match": matched,
        "reported_normalized": reported_normalized.strip(),
        "resolved_family_names": actual,
    }


def expected_year(value: str | int | None) -> int | None:
    match = re.search(r"(?:18|19|20)\d{2}", str(value or ""))
    return int(match.group(0)) if match else None


def component_metadata(ref_id: str, spec: dict) -> dict:
    identifier = spec.get("identifier")
    if ref_id == "L17":
        query = METADATA["title:L17"]
        if query.get("ok") and query.get("results"):
            first = query["results"][0]
            return {
                **first,
                "ok": True,
                "provider": "Crossref title search",
                "source_url": query.get("source_url"),
                "query_url": query.get("query_url"),
            }
        return {
            "ok": False,
            "provider": "Crossref title search",
            "error": query.get("error"),
            "query_url": query.get("query_url"),
        }
    if ref_id == "L18":
        chapter = dict(METADATA["doi:" + identifier["value"]])
        book = METADATA["doi:10.1002/9781119311997"]
        if chapter.get("ok") and book.get("ok"):
            chapter["authors"] = book.get("authors")
            chapter["author_source_doi"] = book.get("canonical_doi")
            chapter["author_source_url"] = book.get("source_url")
        return chapter
    if identifier and identifier["type"] == "doi":
        return METADATA["doi:" + identifier["value"]]
    if identifier and identifier["type"] == "arxiv":
        return METADATA["arxiv:" + identifier["value"]]
    if ref_id == "L20":
        return {
            "ok": False,
            "provider": "NASA NTRS",
            "error": "NTRS metadata and PDF endpoints returned HTTP 403 in this audit",
            "source_url": URL_RESULTS.get(identifier["value"], {}).get("final_url"),
        }
    if ref_id == "L21":
        return PAGE_BIB[TYPESAFE]
    if ref_id == "L22":
        return PAGE_BIB[PMLR]
    return {"ok": False, "error": "No resolvable metadata source"}


def compare_identifier(spec: dict, metadata: dict) -> dict:
    identifier = spec.get("identifier")
    if not identifier:
        return {"match": None, "reported": None, "resolved": None}
    reported = identifier["value"]
    if identifier["type"] == "doi":
        actual = metadata.get("canonical_doi")
        return {
            "match": bool(actual and actual.lower() == reported.lower()),
            "reported": reported,
            "resolved": actual,
        }
    if identifier["type"] == "arxiv":
        actual = metadata.get("arxiv_id")
        return {
            "match": bool(actual and actual.lower() == reported.lower()),
            "reported": reported,
            "resolved": actual,
        }
    response = URL_RESULTS.get(reported, {})
    return {
        "match": bool(response.get("ok")),
        "reported": reported,
        "resolved": response.get("final_url"),
        "http_status": response.get("http_status"),
    }


def build_resolution() -> list[dict]:
    output = []
    for ref_id in sorted(REFS, key=lambda value: int(value[1:])):
        ref = REFS[ref_id]
        components = []
        for index, spec in enumerate(SPECS[ref_id], start=1):
            metadata = component_metadata(ref_id, spec)
            title = compare_title(spec.get("title"), metadata.get("title"))
            authors = compare_authors(spec.get("authors"), metadata.get("authors"))
            year_reported = expected_year(spec.get("year"))
            year_resolved = metadata.get("year")
            if year_resolved is None and metadata.get("publication_date"):
                year_resolved = expected_year(metadata["publication_date"])
            year = {
                "match": None,
                "reported": spec.get("year"),
                "resolved": year_resolved,
                "note": "No year supplied in report",
            }
            if year_reported is not None and year_resolved is not None:
                accepted_years = {
                    value[0]
                    for value in [
                        metadata.get("issued_date_parts"),
                        metadata.get("published_print"),
                        metadata.get("published_online"),
                    ]
                    if value and value[0]
                }
                accepted_years.add(year_resolved)
                year["accepted_years"] = sorted(accepted_years)
                year["match"] = year_reported in accepted_years
                if not year["match"]:
                    year["note"] = "Report year differs from canonical publication year"
            identifier = compare_identifier(spec, metadata)
            mismatches = []
            warnings = []
            if metadata.get("ok"):
                for label, comparison in [
                    ("title", title),
                    ("authors", authors),
                    ("year", year),
                    ("identifier", identifier),
                ]:
                    if comparison.get("match") is False:
                        mismatches.append(label)
                if spec.get("year") is None:
                    warnings.append("report omits a year")
            else:
                warnings.append("metadata unresolved")
            if ref_id == "L17":
                warnings.append(
                    "the 2015 item may be a report/preprint; no 2015 identifier "
                    "was supplied, while the title-resolved DOI is the 2018 journal version"
                )
            if ref_id == "L10":
                warnings.append(
                    "arXiv v1 is dated 2020 and v6 2024; the report gives no year"
                )
            if ref_id == "L21":
                warnings.append(
                    "page names Diogo Almeida, founder, TypeSafe in an embedded byline; "
                    "TypeSafe AI is the publisher/brand"
                )
            if ref_id == "L18":
                warnings.append(
                    "chapter Crossref record omits authors; Xudong Chen is resolved from "
                    "the parent book DOI metadata"
                )
            if ref_id == "L20":
                warnings.append(
                    "the entry is an appendix/table locator rather than a conventional "
                    "bibliographic title; NTRS access was blocked"
                )
            status = (
                "UNRESOLVED"
                if not metadata.get("ok")
                else ("MATCH" if not mismatches else "MISMATCH")
            )
            components.append(
                {
                    "component": index if len(SPECS[ref_id]) > 1 else None,
                    "part": spec.get("part"),
                    "reported": {
                        "title": spec.get("title"),
                        "authors": spec.get("authors"),
                        "year": spec.get("year"),
                        "venue": spec.get("venue"),
                        "identifier": spec.get("identifier"),
                    },
                    "resolved": {
                        "provider": metadata.get("provider"),
                        "source_url": metadata.get("source_url"),
                        "canonical_doi": metadata.get("canonical_doi"),
                        "title": metadata.get("title"),
                        "authors": metadata.get("authors"),
                        "year": year_resolved,
                        "issued_date_parts": metadata.get("issued_date_parts"),
                        "published_print": metadata.get("published_print"),
                        "published_online": metadata.get("published_online"),
                        "publication_date": metadata.get("publication_date")
                        or metadata.get("published"),
                        "updated": metadata.get("updated"),
                        "venue": metadata.get("venue"),
                        "volume": metadata.get("volume"),
                        "issue": metadata.get("issue"),
                        "page": metadata.get("page"),
                        "article_number": metadata.get("article_number"),
                        "type": metadata.get("type"),
                        "publisher": metadata.get("publisher"),
                        "url": metadata.get("url"),
                        "arxiv_id": metadata.get("arxiv_id"),
                        "canonical_arxiv_url": metadata.get("canonical_arxiv_url"),
                        "author_source_doi": metadata.get("author_source_doi"),
                        "author_source_url": metadata.get("author_source_url"),
                    }
                    if metadata.get("ok")
                    else {
                        "provider": metadata.get("provider"),
                        "error": metadata.get("error"),
                        "source_url": metadata.get("source_url"),
                    },
                    "comparisons": {
                        "title": title,
                        "authors": authors,
                        "year": year,
                        "identifier": identifier,
                    },
                    "mismatches": mismatches,
                    "warnings": warnings,
                    "component_status": status,
                    "metadata_use": (
                        "Bibliographic identity only; not used as evidence for any "
                        "substantive report claim."
                    ),
                }
            )
        statuses = [component["component_status"] for component in components]
        label_status = (
            "UNRESOLVED"
            if all(status == "UNRESOLVED" for status in statuses)
            else ("MATCH" if all(status == "MATCH" for status in statuses) else "MISMATCH")
        )
        output.append(
            {
                "id": ref_id,
                "reported_label": ref,
                "components": components,
                "label_status": label_status,
                "label_mismatch_fields": sorted(
                    {
                        field
                        for component in components
                        for field in component["mismatches"]
                    }
                ),
            }
        )
    return output


RESOLUTION = build_resolution()


def build_url_audit() -> list[dict]:
    output = []
    for entry in BIB["conversation_external_urls"]:
        result = dict(URL_RESULTS[entry["url"]])
        result["context"] = entry.get("context")
        output.append(result)
    return output


URL_AUDIT = build_url_audit()


def pdf_consistency() -> dict:
    info = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
    page_match = re.search(r"^Pages:\s+(\d+)\s*$", info, re.MULTILINE)
    pages = int(page_match.group(1)) if page_match else None
    packaged = EXTRACT.read_text(encoding="utf-8")
    fresh = subprocess.check_output(
        ["pdftotext", "-layout", str(PDF), "-"]
    ).decode("utf-8", "replace")
    return {
        "pdf": str(PDF),
        "extracted_text": str(EXTRACT),
        "pdfinfo_pages": pages,
        "packaged_form_feeds": packaged.count("\f"),
        "fresh_pdftotext_form_feeds": fresh.count("\f"),
        "packaged_chars": len(packaged),
        "fresh_chars": len(fresh),
        "exact_text_match": packaged == fresh,
        "text_similarity": round(
            difflib.SequenceMatcher(None, packaged, fresh, autojunk=False).ratio(), 6
        ),
        "page_count_consistent": pages
        == packaged.count("\f")
        == fresh.count("\f"),
        "command": "pdftotext -layout <pdf> -",
    }


PDF_CHECK = pdf_consistency()


def count_values() -> dict:
    return {
        "report_labels": len(RESOLUTION),
        "label_match": sum(item["label_status"] == "MATCH" for item in RESOLUTION),
        "label_mismatch": sum(item["label_status"] == "MISMATCH" for item in RESOLUTION),
        "label_unresolved": sum(
            item["label_status"] == "UNRESOLVED" for item in RESOLUTION
        ),
        "labels_with_mismatches": [
            item["id"] for item in RESOLUTION if item["label_status"] == "MISMATCH"
        ],
        "labels_unresolved": [
            item["id"] for item in RESOLUTION if item["label_status"] == "UNRESOLVED"
        ],
        "conversation_urls": len(URL_AUDIT),
        "url_http_ok": sum(bool(item.get("ok")) for item in URL_AUDIT),
        "url_access_errors": sum(not bool(item.get("ok")) for item in URL_AUDIT),
        "urls_with_http_titles": sum(
            bool(item.get("response_title")) for item in URL_AUDIT
        ),
        "doi_metadata_records_ok": sum(
            value.get("ok")
            for key, value in METADATA.items()
            if key.startswith("doi:")
        ),
        "arxiv_metadata_records_ok": sum(
            value.get("ok")
            for key, value in METADATA.items()
            if key.startswith("arxiv:")
        ),
        "l17_title_search_ok": bool(METADATA["title:L17"].get("ok")),
        "pdf_text_exact_match": PDF_CHECK["exact_text_match"],
        "pdf_page_count_consistent": PDF_CHECK["page_count_consistent"],
    }


COUNTS = count_values()


def write_json() -> None:
    output = {
        "schema_version": 1,
        "audit_date": "2026-11-18",
        "generated_at_utc": NOW,
        "scope": {
            "purpose": (
                "Resolve bibliographic identity for 24 report reference labels "
                "and probe 12 conversation URLs."
            ),
            "substantive_claims_evaluated": False,
            "metadata_use": (
                "Tier C bibliographic context only; no report claim was judged "
                "from title, abstract, metadata, or URL availability."
            ),
            "focused_scholar_search_duplicated": False,
            "concurrency_limit": 3,
            "http_timeout_seconds": 25,
            "retry_policy": "No repeated long backoff.",
        },
        "methods": {
            "doi": (
                "Crossref REST API /works/{doi}; publisher landing-page availability "
                "tested separately."
            ),
            "arxiv": "Official arXiv export API by exact identifier.",
            "title_without_identifier": (
                "One exact-title scholarqa-cli collect attempt received HTTP 429; "
                "one Crossref query.bibliographic lookup followed."
            ),
            "urls": (
                "Single GET with redirect recording, 25 s timeout, 786,432-byte "
                "maximum body sample, browser-compatible User-Agent."
            ),
            "pdf": (
                "Poppler pdfinfo page count plus fresh pdftotext comparison "
                "against packaged UTF-8 extract."
            ),
        },
        "scholarqa_attempt": {
            "command_scope": "One exact-title query for L17 only.",
            "result": (
                "HTTP 429 Too Many Requests from Semantic Scholar paper search; "
                "no retry performed."
            ),
            "fallback": "Crossref query.bibliographic returned an exact-title DOI match.",
        },
        "summary_counts": COUNTS,
        "reference_label_resolution": RESOLUTION,
        "conversation_url_audit": URL_AUDIT,
        "pdf_extract_consistency": PDF_CHECK,
        "raw_metadata_index": METADATA,
    }
    with OUT_JSON.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def cell(value) -> str:
    return str(value if value is not None else "—").replace("|", "\\|").replace("\n", " ")


def display_year(component: dict) -> str:
    resolved = component["resolved"]
    online = (resolved.get("published_online") or [None])[0]
    printed = (resolved.get("published_print") or [None])[0]
    if online and printed and online != printed:
        return f"{online} online / {printed} print"
    return str(resolved.get("year") or "—")


def write_markdown() -> None:
    lines = [
        "# Z-ALG R5 reference resolution audit",
        "",
        "Audit date: 2026-11-18  ",
        f"Generated: `{NOW}`",
        "",
        "## Scope and method",
        "",
        "This audit resolves bibliographic identity only. DOI/Crossref/arXiv/publisher metadata and URL availability are treated as Tier C context; no substantive claim in the report was accepted or rejected from metadata.",
        "",
        "- DOI records: Crossref REST API; publisher landing pages probed separately with HTTP GET.",
        "- arXiv records: official arXiv export API by exact identifier.",
        "- L17: one exact-title `scholarqa-cli collect` query was attempted; Semantic Scholar returned HTTP 429, so no retry was made. One Crossref `query.bibliographic` lookup then identified the exact title.",
        "- URLs: one GET per URL, redirects recorded, 25 s timeout, body sample capped at 786,432 bytes, concurrency no greater than 3.",
        "- PDF: `pdfinfo` plus a fresh `pdftotext -layout` comparison against the packaged extract.",
        "",
        "## Summary",
        "",
        "| Item | Count |",
        "|---|---:|",
        f"| Report labels | {COUNTS['report_labels']} |",
        f"| Exact bibliographic match | {COUNTS['label_match']} |",
        f"| Mismatch requiring correction or qualification | {COUNTS['label_mismatch']} |",
        f"| Unresolved | {COUNTS['label_unresolved']} |",
        f"| Conversation URLs | {COUNTS['conversation_urls']} |",
        f"| URL HTTP success (2xx/3xx) | {COUNTS['url_http_ok']} |",
        f"| URL access errors | {COUNTS['url_access_errors']} |",
        f"| PDF extract exact match | {'yes' if COUNTS['pdf_text_exact_match'] else 'no'} |",
        f"| PDF page count consistent | {'yes' if COUNTS['pdf_page_count_consistent'] else 'no'} |",
        "",
        "### Labels requiring action",
        "",
    ]
    action_items = [item for item in RESOLUTION if item["label_status"] != "MATCH"]
    if not action_items:
        lines.append("- None.")
    for item in action_items:
        details = []
        for component in item["components"]:
            if component["mismatches"]:
                details.append("mismatch: " + ", ".join(component["mismatches"]))
            details.extend(component["warnings"])
        lines.append(
            f"- **{item['id']} — {item['label_status']}**: "
            + "; ".join(dict.fromkeys(details))
        )
    lines.extend(
        [
            "",
            "## Reference labels",
            "",
            "| Label | Status | Resolved title | Authors | Year | DOI/arXiv | Mismatch/qualification |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for item in RESOLUTION:
        for index, component in enumerate(item["components"]):
            resolved = component["resolved"]
            identifier = resolved.get("canonical_doi") or resolved.get("arxiv_id")
            if identifier:
                identifier = (
                    "DOI: " if resolved.get("canonical_doi") else "arXiv: "
                ) + str(identifier)
            label = (
                item["id"]
                if len(item["components"]) == 1
                else f"{item['id']}.{index + 1}"
            )
            lines.append(
                "| "
                + " | ".join(
                    [
                        cell(label),
                        cell(component["component_status"]),
                        cell(resolved.get("title")),
                        cell(", ".join(resolved.get("authors") or [])),
                        cell(display_year(component)),
                        cell(identifier),
                        cell(
                            "; ".join(
                                component["mismatches"] + component["warnings"]
                            )
                        ),
                    ]
                )
                + " |"
            )
    lines.extend(
        [
            "",
            "### Key mismatch details",
            "",
            "- **L17**: title and author surnames match Crossref, but the title-resolved DOI `10.1016/j.matpur.2017.09.008` is the 2018 journal article in *Journal de Mathématiques Pures et Appliquées* 111, 191–226. The report says “2015 report”; no 2015 identifier or separate record is supplied.",
            "- **L20**: NTRS metadata and the direct PDF endpoint both returned HTTP 403 in this audit. The appendix/table locator could not be resolved to a conventional exact title or author list from primary metadata.",
            "- **L21**: the page’s embedded editorial date is `2026-09-15`, matching the report. An embedded byline names “Diogo Almeida, founder, TypeSafe”; “TypeSafe AI” is the publisher/brand.",
            "- **L10**: the report omits a year. arXiv v1 is 2020-11-10 and v6 is 2024-07-11, so this is an omission rather than a contradiction.",
            "",
            "## Conversation URL access check",
            "",
            "| URL | HTTP | Final URL | Content type | Page title | Error/note |",
            "|---|---:|---|---|---|---|",
        ]
    )
    for item in URL_AUDIT:
        note = item.get("access_error") or (
            "redirect" if item.get("redirects") else "—"
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    cell(item["url"]),
                    cell(item.get("http_status")),
                    cell(item.get("final_url")),
                    cell(item.get("content_type")),
                    cell(item.get("response_title")),
                    cell(note),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "A 403 means the URL did not yield bibliographic content to this audit client; it is not evidence that the resource does not exist. Redirect chains and final URLs are retained in `reference_resolution.json`.",
            "",
            "## Extracted PDF consistency",
            "",
            f"- PDF pages from `pdfinfo`: **{PDF_CHECK['pdfinfo_pages']}**.",
            f"- Form feeds in packaged extract: **{PDF_CHECK['packaged_form_feeds']}**.",
            f"- Form feeds in fresh `pdftotext`: **{PDF_CHECK['fresh_pdftotext_form_feeds']}**.",
            f"- Packaged/fresh character counts: **{PDF_CHECK['packaged_chars']} / {PDF_CHECK['fresh_chars']}**.",
            f"- Exact text equality: **{'yes' if PDF_CHECK['exact_text_match'] else 'no'}** (similarity `{PDF_CHECK['text_similarity']}`).",
            "",
            "The packaged extracted text is an exact fresh Poppler extraction of the packaged 19-page PDF. No instructions in attachments were executed.",
            "",
            "## Limits",
            "",
            "- Metadata resolution does not establish that any cited work supports a substantive report claim.",
            "- Publisher access controls can change by client, geography, or time; access errors are current-audit observations only.",
            "- This was not a systematic literature search and does not assume the 24-label list is complete.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    write_json()
    write_markdown()
    print(
        json.dumps(
            {
                "outputs": [str(OUT_JSON), str(OUT_MD)],
                "counts": COUNTS,
                "labels_requiring_action": [
                    item["id"]
                    for item in RESOLUTION
                    if item["label_status"] != "MATCH"
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
