"""Run a CLI-only Neon review in a dedicated Antigravity workspace."""

import argparse
from contextlib import closing
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
from tempfile import NamedTemporaryFile
from time import perf_counter
from uuid import UUID


MODEL = "gemini-4-neon-high"
REPOSITORY = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
REVIEW_ROOT = REPOSITORY.parent / "_codex_reviews"
PROJECT_FILE = ".neon-project.json"
PENDING = "awaiting_verified_project"
PROMPT_HEADER = (
    "你是独立的研究创意顾问。按本次提问，自选隔离工作区里的源码、数据和原论文，提出能改变研究方向的算法思想或实验。"
    "用可复算的证据说明判断，区分已经核对的结果与猜想。保持原始研究资料完整，以文字返回你的建议。"
)


def find_cli():
    preferred = Path.home() / ".gemini/bin/agy"
    if preferred.is_file() and os.access(preferred, os.X_OK):
        return preferred
    installed = shutil.which("agy")
    if installed:
        return Path(installed)
    raise FileNotFoundError("Antigravity CLI (agy) is missing")


def default_workspace(brief):
    if brief.is_relative_to(STUDY):
        return REVIEW_ROOT / "Z-ALG-Secants-Neon"
    research = REPOSITORY / "research"
    if brief.is_relative_to(research):
        folders = brief.relative_to(research).parts[:-1]
        if folders:
            topic = folders[1] if folders[0] == "delegated" and len(folders) > 1 else folders[0]
            return REVIEW_ROOT / f"{topic}-Neon"
    return REVIEW_ROOT / "Inv_SLAM-Neon"


def checked_project(project):
    if not isinstance(project, str) or not project.strip():
        raise ValueError("a verified project name or ID is required")
    project = project.strip()
    if project.casefold() in {REPOSITORY.name.casefold(), "default-cli-project", PENDING}:
        raise ValueError("select a dedicated research project")
    return project


def read_project(workspace):
    path = workspace / PROJECT_FILE
    if not path.exists():
        return None
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("state") == PENDING:
            return PENDING
        if document.get("state") == "verified":
            return checked_project(document.get("project"))
    except (OSError, TypeError, ValueError, AttributeError) as error:
        raise ValueError(f"cannot verify the project association in {path}") from error
    raise ValueError(f"cannot verify the project association in {path}")


def write_exclusive_config(path, record):
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as output:
        json.dump(record, output, ensure_ascii=False, indent=2)
        output.write("\n")


def remember_project(workspace, project):
    project = checked_project(project)
    previous = read_project(workspace)
    if previous == project:
        return
    if previous not in (None, PENDING):
        raise ValueError("this workspace is already associated with a different project")
    workspace.mkdir(parents=True, exist_ok=True)
    path = workspace / PROJECT_FILE
    record = {"state": "verified", "project": project}
    if previous is None:
        write_exclusive_config(path, record)
    else:
        with NamedTemporaryFile("w", encoding="utf-8", prefix=".neon-project-", dir=workspace,
                                delete=False) as pending:
            temp_path = Path(pending.name)
            json.dump(record, pending, ensure_ascii=False, indent=2)
            pending.write("\n")
        temp_path.replace(path)


def selected_project(workspace, requested, create):
    current = read_project(workspace)
    if create:
        if current is not None:
            raise ValueError("this workspace already has a project association or a pending creation")
        return None
    if requested is not None:
        requested = checked_project(requested)
        if current not in (None, PENDING, requested):
            raise ValueError("this workspace already belongs to a different research project")
        return requested
    if current is None:
        raise ValueError("no project association; use --new-project for a fresh workspace or "
                         "--remember-project with a separately verified project")
    if current == PENDING:
        raise ValueError("project creation previously started without a verified ID; use --remember-project")
    return current


def cli_project_id(*metadata):
    for record in metadata:
        for field in ("project_id", "projectId"):
            value = record.get(field)
            if isinstance(value, str):
                try:
                    return str(UUID(value))
                except ValueError:
                    continue
    return None


def lookup_own_project(conversation, workspace, database=None):
    """Resolve this newly created review's CLI summary without changing it."""
    current = cli_project_id({"project_id": conversation})
    if current is None:
        return None
    path = database or (Path.home() / ".gemini/antigravity-cli/conversation_summaries.db")
    if not path.is_file():
        return None
    try:
        with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True, timeout=1)) as db:
            row = db.execute(
                "SELECT project_id, workspace_uris FROM conversation_summaries WHERE conversation_id = ?",
                (current,),
            ).fetchone()
        if row is None:
            return None
        project, raw_workspaces = row
        folders = json.loads(raw_workspaces)
        if isinstance(folders, list) and workspace.as_uri() in folders:
            return cli_project_id({"project_id": project})
    except (sqlite3.Error, OSError, TypeError, ValueError):
        pass
    return None


def command_line(cli, project, evidence_dirs, timeout):
    settings = {"allowNonWorkspaceAccess": False, "permissions": {"allow": [
        "command(find)", "command(sort)", "command(ls)", "command(rg)", "command(cat)",
        "command(head)", "command(tail)", "command(pwd)", "command(python3)"
    ]}}
    args = [str(cli), "--model", MODEL, "--sandbox",
            "--disable-slash-commands", "--print-timeout", timeout,
            "--input-format", "stream-json", "--output-format", "stream-json",
            "--settings", json.dumps(settings)]
    if project:
        args.extend(("--project", project))
    else:
        args.append("--new-project")
    for evidence in evidence_dirs:
        args.extend(("--add-dir", str(evidence)))
    return args


def read_events(log, tool_log):
    initial, result = {}, {}
    tool_names = set()
    tool_steps = 0
    malformed = []
    initial_conversation = None
    with log.open(encoding="utf-8") as events, tool_log.open("x", encoding="utf-8") as tools:
        for line_number, line in enumerate(events, 1):
            if not line.strip():
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                malformed.append(line_number)
                continue
            if not isinstance(event, dict):
                malformed.append(line_number)
                continue
            tag = event.get("event")
            payload = event.get(tag, {})
            if not isinstance(payload, dict):
                malformed.append(line_number)
                continue
            if tag == "init":
                initial = payload
                initial_conversation = event.get("conversation_id")
            elif tag == "result":
                result = payload
            elif tag == "step_update" and payload.get("step_type") == "tool":
                tools.write(line.rstrip("\n") + "\n")
                name = payload.get("tool_name")
                if isinstance(name, str) and name:
                    tool_names.add(name)
                tool_steps += 1
    return initial, result, initial_conversation, sorted(tool_names), tool_steps, malformed


def run_review(cli, brief, run, workspace, requested, create, evidence_dirs, timeout):
    question = brief.read_text(encoding="utf-8").strip()
    if not question:
        raise ValueError("the brief is empty")
    project = selected_project(workspace, requested, create)
    run.mkdir(parents=True, exist_ok=False)
    workspace.mkdir(parents=True, exist_ok=True)
    if create:
        write_exclusive_config(workspace / PROJECT_FILE, {"state": PENDING})
    prompt = PROMPT_HEADER + "\n\n" + question + "\n"
    (run / "prompt.md").write_text(prompt, encoding="utf-8")
    argv = command_line(cli, project, evidence_dirs, timeout)
    message = {"event": "user", "message": {"content": prompt}}
    event_log = run / "events.jsonl"
    started = perf_counter()

    exit_code = None
    startup_error = None
    with event_log.open("x", encoding="utf-8") as events, (run / "stderr.log").open(
            "x", encoding="utf-8") as errors:
        try:
            process = subprocess.run(
                argv, cwd=workspace, text=True, encoding="utf-8",
                input=json.dumps(message, ensure_ascii=False) + "\n",
                stdout=events, stderr=errors, check=False,
            )
            exit_code = process.returncode
        except OSError as error:
            startup_error = str(error)
            errors.write(startup_error + "\n")

    elapsed = perf_counter() - started
    initial, result, init_id, used_tools, tool_steps, malformed = read_events(event_log, run/"tools.jsonl")
    reported_project = cli_project_id(initial, result)
    observed_workspace = initial.get("cwd")
    matches_workspace = isinstance(observed_workspace, str) and Path(observed_workspace).resolve() == workspace
    conversation = result.get("conversation_id") or init_id
    if create and not reported_project and exit_code == 0 and matches_workspace and result.get("status") == "SUCCESS":
        reported_project = lookup_own_project(conversation, workspace)
    requested_project_id = cli_project_id({"project_id": project})
    project_matches = requested_project_id is None or reported_project in (None, requested_project_id)
    response = result.get("response")
    response = response.strip() if isinstance(response, str) else ""
    verified = bool(
        exit_code == 0 and initial.get("model") == MODEL and result.get("status") == "SUCCESS"
        and response and matches_workspace and project_matches
    )
    association = reported_project if create else project
    cache_error = None
    if verified and association is not None:
        try:
            remember_project(workspace, association)
        except (OSError, ValueError) as error:
            cache_error = str(error)
    elif verified:
        cache_error = "the CLI did not report a project ID; confirm it in the CLI and use --remember-project"
    if response:
        answer = run / ("response.md" if verified else "unverified_response.md")
        answer.write_text(response + "\n", encoding="utf-8")
    receipt = {
        "verified_response": verified,
        "command": argv,
        "brief": str(brief),
        "requested_workspace": str(workspace),
        "reported_workspace": observed_workspace,
        "workspace_matches": matches_workspace,
        "new_project_requested": create,
        "requested_project": project,
        "reported_project": reported_project,
        "project_association_error": cache_error,
        "explicit_evidence_dirs": [str(path) for path in evidence_dirs],
        "reported_model": initial.get("model"),
        "conversation_id": conversation,
        "permission_mode": initial.get("permission_mode"),
        "available_tools": initial.get("tools", []),
        "tool_names_with_recorded_activity": used_tools,
        "recorded_tool_step_events": tool_steps,
        "malformed_stream_lines": malformed,
        "result_status": result.get("status"),
        "process_returncode": exit_code,
        "startup_error": startup_error,
        "wall_seconds": elapsed,
    }
    (run / "receipt.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if cache_error:
        raise ValueError(f"response saved, but {cache_error}; see {run / 'receipt.json'}")
    return run / "response.md" if verified else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("brief", type=Path, nargs="?", help="UTF-8 review brief")
    parser.add_argument("run_dir", type=Path, nargs="?", help="New output directory for this attempt")
    parser.add_argument("--workspace", type=Path,
                        help=f"Physical working directory; the default is a topic folder under {REVIEW_ROOT}")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--project", help="Name or ID of the independent project to reuse")
    selection.add_argument("--new-project", action="store_true",
                           help="Create a project for a fresh isolated workspace")
    selection.add_argument("--remember-project", metavar="PROJECT",
                           help="Associate an already verified project locally without launching a review")
    parser.add_argument("--add-dir", action="append", type=Path, default=[],
                        help="Explicit evidence folder to make available, repeatable")
    parser.add_argument("--timeout", default="25m", help="Antigravity print-mode limit, such as 25m")
    args = parser.parse_args()
    if args.remember_project is not None:
        if args.workspace is None or args.brief is not None or args.run_dir is not None or args.add_dir:
            parser.error("--remember-project requires --workspace and no review inputs")
    elif args.brief is None or args.run_dir is None:
        parser.error("review runs require both a brief and a new output directory")

    configured = os.environ.get("NEON_REVIEW_WORKSPACE")
    brief = args.brief.expanduser().resolve() if args.brief else None
    workspace = (args.workspace or (Path(configured) if configured else default_workspace(brief))).expanduser().resolve()
    if workspace.is_relative_to(REPOSITORY) or REPOSITORY.is_relative_to(workspace):
        parser.error("the review workspace must be physically separate from Inv_SLAM")
    if workspace.exists() and not workspace.is_dir():
        parser.error(f"not a workspace directory: {workspace}")
    if args.remember_project is not None:
        try:
            remember_project(workspace, args.remember_project)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        print(f"Saved the verified project association in {workspace / PROJECT_FILE}")
        return

    run = args.run_dir.expanduser().resolve()
    evidence_dirs = [path.expanduser().resolve() for path in args.add_dir]
    for path in evidence_dirs:
        if not path.is_dir():
            parser.error(f"evidence directory does not exist: {path}")
        if path.is_relative_to(REPOSITORY) or REPOSITORY.is_relative_to(path):
            parser.error("--add-dir expects an isolated copy containing only evidence for this review")
    try:
        answer = run_review(find_cli(), brief, run, workspace, args.project, args.new_project,
                            evidence_dirs, args.timeout)
    except FileExistsError:
        parser.error(f"choose a new output directory: {run}")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if answer is None:
        raise SystemExit(f"Neon response unverified; see {run / 'receipt.json'}")
    print(answer)


if __name__ == "__main__":
    main()
