"""Clear only explicitly registered files made in a new run's private scratch directory."""
import argparse
import json
import sys
from pathlib import Path


NAME = ".maxwell-scratch"
LEDGER = ".temporary-files.jsonl"


def scratch_for(run):
    directory = run / NAME
    if directory.is_symlink() or not (directory / LEDGER).is_file():
        raise ValueError("没有本工具为本轮新建的临时目录")
    if (directory / LEDGER).is_symlink():
        raise ValueError("清单不能是符号链接")
    return directory


def owned_file(directory, name):
    rel = Path(name)
    if rel.is_absolute() or not rel.parts or ".." in rel.parts or rel.name == LEDGER:
        raise ValueError(f"范围外路径 {name}")
    file = directory / rel
    if file.resolve() != directory.resolve() / rel or file.is_symlink():
        raise ValueError(f"跳过链接或重解析路径 {name}")
    return file


def same_contents(left, right):
    with left.open("rb") as source, right.open("rb") as kept:
        while True:
            block = source.read(1024 * 1024)
            if kept.read(len(block) or 1) != block:
                return False
            if not block:
                return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for name in ("init", "mark", "clean"):
        p = sub.add_parser(name)
        p.add_argument("run", type=Path, help="本轮输出目录")
        if name == "mark":
            p.add_argument("file", help="相对于 .maxwell-scratch 的单个普通文件")
            group = p.add_mutually_exclusive_group(required=True)
            group.add_argument("--copy-at", type=Path, help="当前任务已完成复制、预定保留的文件")
            group.add_argument("--rebuild", help="从保留输入重新生成此工作件的实际命令")
        if name == "clean":
            p.add_argument("--apply", action="store_true", help="作业终止和归档完成后执行")
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    if args.action == "init":
        directory = run / NAME
        directory.mkdir()  # Never take ownership of an existing directory.
        (directory / LEDGER).touch()
        print(directory)
        return 0
    directory = scratch_for(run)
    ledger = directory / LEDGER
    if args.action == "mark":
        file = owned_file(directory, args.file)
        if not file.is_file():
            raise ValueError("只登记已生成的普通临时文件")
        st = file.stat()
        row = dict(name=str(file.relative_to(directory)), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
        if args.copy_at:
            keep = args.copy_at.resolve(strict=True)
            if directory in keep.parents or not keep.is_file() or keep.stat().st_size != st.st_size:
                raise ValueError("保留副本应在临时区以外，且与成功复制记录大小相符")
            row["copy_at"] = str(keep)
        elif args.rebuild.strip():
            row["rebuild"] = args.rebuild.strip()
        else:
            raise ValueError("缺少重建办法")
        with ledger.open("a", encoding="utf8") as stream:
            print(json.dumps(row, ensure_ascii=False), file=stream)
        return 0
    entries = {r["name"]: r for r in (json.loads(s) for s in ledger.read_text(encoding="utf8").splitlines())}
    paid = 0
    count = skipped = 0
    for name, row in entries.items():
        try:
            file = owned_file(directory, name)
        except ValueError as reason:
            print("保留", reason)
            skipped += 1
            continue
        if not file.exists():  # A previous cleanup may already have finished this entry.
            continue
        st = file.stat()
        if not file.is_file() or (st.st_size, st.st_mtime_ns) != (row["bytes"], row["mtime_ns"]):
            print("保留（登记后发生变化）", file)
            skipped += 1
            continue
        if "copy_at" in row:
            keep = Path(row["copy_at"])
            if (not keep.is_file() or keep.resolve() != keep
                    or directory in keep.parents or keep.stat().st_size != st.st_size):
                print("保留（目前找不到登记的完整副本）", file)
                skipped += 1
                continue
            if args.apply and not same_contents(file, keep):
                print("保留（指定副本的内容有差异）", file)
                skipped += 1
                continue
        if args.apply:
            with (run / "scratch-removed.jsonl").open("a", encoding="utf8") as stream:
                file.unlink()
                print(json.dumps(row, ensure_ascii=False), file=stream)
        count += 1
        paid += st.st_size
        print("已清" if args.apply else "可清", file)
    print(f'{"已处理" if args.apply else "预览"} {count} 个文件，逻辑大小 {paid} B；保留变化项 {skipped}。')
    return int(bool(skipped))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as exc:
        sys.exit(str(exc))
